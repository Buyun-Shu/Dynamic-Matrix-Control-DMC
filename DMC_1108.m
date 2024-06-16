%%离散动态矩阵控制（DMC）
clc;        %清除命令窗口
clear all;  %清除工作区变量
close all;  %关闭所有图表
%% 系统模型建立
% num=5;                          %分子多项式
% den=[50 5];                     %分母多项式
num=0.115;
den=[100 0.603];
% num=[10256.86 0];
% den=[1 0 35.2836];
                                %den=[1 0.5 0.1];（二阶系统的分母多项式）
sys=tf(num,den,'InputDelay',1); %创建传递函数系统模型，输入延迟为1个单位
Ts=1;                           %离散化采样时间间隔，1个单位（采样周期）
dsys=c2d(sys,Ts,'foh');         %c2d将连续时间系统sys根据前向欧拉法和采样时间转换为离散时间系统dsys
[num1,den1]=tfdata(dsys,'v');   %提取离散时间系统的分子、分母系数
% %step(dsys);
% grid on
% %T_final=200;
T_final=4000;
% tt=Ts:Ts:T_final*Ts;            %向量tt包含从Ts到T_final*Ts等间隔时间点
% a1=step(dsys,tt);               %存储离散时间系统dsys在时间向量tt上的步响应数据
loaded_data=load('sudu.mat','sudu');
a1 = loaded_data.sudu;
%% 模拟实时采集管径数据
% n = T_final;                    %指定数组长度
% rangeMin = 180;                 %管径下限mm
% rangeMax = 220;                 %管径上限mm
% real_time_data = (rangeMax-rangeMin).*rand(1, n) + rangeMin;
% 
% shape_parameter = 1;            %形状参数
% sensitivity = 10;               %敏感度参数
% sum_rk = real_time_data(1);
%% 设置控制器参数（可以修改之处）
N=20;         %建模时域，默认NT之后阶跃响应区域平稳(原来N=20）
                %建模时域N的作用是确定控制器"使用多少个最近历史输入和输出数据来进行预测"
                %较大N值可提供更长时间范围的预测，但会增加计算复杂度和延迟
                %N越大曲线越光滑

P=6;           %优化时域，被控对象在控制增量作用下未来P个时刻预测输出值***P=6
                %用于确定控制器生成的最优控制序列的长度
                %较大P值可以更好地处理系统的长时间影响和变化，但会增加计算复杂度

M=2;           %控制时域，每一时刻k起的M个控制增量作用***M=2     P大于等于M
                %表示在每个控制时刻向前预测控制信号的时间步数
                %较大的M值使控制器能够预测未来更多时间步的控制信号，但也会增加计算复杂度和延迟

q=1;            %权系数，抑制跟踪误差变化***q=1
                %较大的q值表示更重视减小输出误差，使得控制器更多调节控制信号以追踪期望参考轨迹，但会引入更多的控制操作。
                %较小的q值则意味着更宽容的系统输出误差，控制器可能更趋向于保持稳定而减小控制操作
     
r=0.1;          %权系数，抑制控制量变化***r=0.1
                %较大r值更重视减小控制输入的幅度，使得控制器更谨慎调节控制信号以减小控制活动，但可能响应较慢（稳定、准确）
                %较小r值则意味着更宽松的控制输入限制，控制器可能更积极地调节控制信号以快速响应系统的需求（快速响应）

%w=2.5;         %期望值***
w=1;
a=a1(1:N);      %模型参数，对象单位阶跃响应采样值

alpha=0.8;      %柔化系数***alpha=0.8
                %柔化系数值越大，控制器输出就越平滑；但过大的柔化系数使得控制系统响应变得迟钝
                %0.9时彻底没有超调量，调节时间大约70s
                %0.2前段不光滑，有较大超调量，调节时间也较长
                %0.1前段不光滑，无超调量，调节时间大约60s
%% 计算动态矩阵A
for i=1:P                   %A是一个P行M列的矩阵
    for j=1:M
        if j>i
            A(i,j)=0;
        else
            A(i,j)=a(i-j+1);%根据模型参数采样值计算出A矩阵
        end
    end
end
%% 权矩阵Q、R、反馈系数H、控制向量d
Q=q*eye(P);                 %误差权矩阵，P行P列的对角矩阵
R=r*eye(M);                 %控制权矩阵，M行M列的对角矩阵
H=[1;alpha*ones(N-1,1)];    %校正向量，修正对未来输出的预测，N行1列的列向量，其中首元素是1，其余都是alpha
c=[1,zeros(1,M-1)];         %M维行向量，取首元素计算
d=c*(inv(A'*Q*A+R)*A'*Q);   %控制向量d，一个1行M列的向量
%% 构造移位矩阵S
S=zeros(N);                 %S是一个N行N列的矩阵，用于构建DMC算法中延迟输入的移位操作
for i=1:N
    if i<N
        S(i,i+1)=1;
    end
    if i==N
        S(N,N)=1;
    end
end
%% 初始化DMC
Alpha=[];                   %初始化Alpha向量
for j=1:P
    Alpha=[Alpha,alpha^j];  %存储柔化系数的1~P幂次方
end
Uk=zeros(1,M);              %控制输入（滚动优化后的实际输入）
Yp0k=[];                    %预测从k+1到k+P未来输出，控制量变化0次
for i=1:P
    Yp0k=[Yp0k;0];          
end
Yn0k=[];                    %延迟输入向量
for i=1:N
    Yn0k=[Yn0k;0];          %N行1列的列向量        
end
Yk=zeros(1,N);              %1行N列的零矩阵
uk=0;                       %实际控制输入
yk=0;                       %当前状态
utemp=0;                    %临时变量，用于存储实际控制输入uk的当前更新值
%% 梯度下降法初始化

% 初始化学习率
% learning_rate = 0.01;       % 可根据实际情况进行调整
% max_iterations = 100;       % 最大迭代次数
%% DMC实时控制
for i=1:T_final
    time(i)=i*Ts;    
    y1k=-den1(2)*yk+num1(2)*uk;     %实测输出y1k，当前状态yk 实际控制输入uk？？？
    %y2k=-den1(3)*yk-den1(2)*y1k+num1(3)uk;（可能能用于二阶线性模型）
    
    %加入脉冲响应
%     if(i==500)
%         y1k=y1k+0.3;
%     end
%     if(i==515)
%         y1k=y1k+0.4;
%     end
%     if(i==580)
%         y1k=y1k-0.8;
%     end
%     if(i==600)
%         y1k=y1k-0.3;
%     end

    Yk(i)=y1k;                      %将y1k存储在Yk向量中
    Uk(i)=utemp;                    %将当前实际控制输入uk的更新值存储在Uk向量中
    Yr=(1-Alpha')*w+Alpha'*y1k;     %期望输出预测值（或参考轨迹）
    error=y1k-Yp0k(1);              %输出误差，实际输出y1k与-模型预测输出向量Yp0k首元素   【3-10】
    Ycork=Yn0k+H*error;             %Ycork为校正后的输出预测向量，H为校正向量            【3-11】
    Yn0k=S*Ycork;                   %延迟输入Yn0k，S为移位矩阵                          【3-13】
    Yp0k=Yn0k(1:P);                 %更新 预测输出向量Yp0k，截取Yn0k的前P个元素
    dU=d*(Yr-Yp0k);                 %控制输入增量dU，d为控制向量                         【3-7】

    utemp=uk+dU;                    %更新临时变量utemp，得到 新的实际控制输入uk的值
    Yp1k=Yp0k+a(1:P)*dU;            %计算新的预测输出向量Yp1k
    Yn1k=Yn0k+a(1:N)*dU;            %计算新的延迟输入向量Yn1k
    Yp0k=Yp1k;                      %更新预测输出向量Yp0k，将Yp1k的值赋给Yp0k
    Yn0k=Yn1k;                      %更新延迟输入向量Yn0k，将Yn1k的值赋给Yn0k
    uk=utemp;                       %更新实际控制输入uk的值
    yk=y1k;                         %更新输出yk的值
    
%     if i>=2
%         delta_rk = real_time_data(i) - real_time_data(i-1);
%         sum_rk = sum_rk + real_time_data(i);
%         equal_rk = sum_rk/i;
%         weighting_function = 1/(1+shape_parameter*((abs(delta_rk/equal_rk))^sensitivity));
%         r = r * weighting_function;
%     end
%     %error1=Yk-(Yn0k+([1;alpha*ones(N-1,1)]).*error);        
%     objective = @(alpha) sum((Yk' - (Yn0k + ([1; alpha.*ones(N-1,1)]).*error)).^2)*(1/i);%自定义的目标函数，用于评估参数效果
%    
%     %Yk是1行N列的零矩阵，Yn0k是N行1列的列向量，H是N行1列的列向量
% 
%     % 梯度下降法优化过程
%     for iter = 1:max_iterations
%         % 使用数值微分近似计算梯度
% %         grad_Q = numerical_gradient(q, error1);
% %         grad_R = numerical_gradient(r, error1);
%         grad_alpha = numerical_gradient(alpha, objective);
%         
%         % 更新控制器参数
% %         q = q - learning_rate * grad_Q;
% %         r = r - learning_rate * grad_R;
%         alpha = alpha - learning_rate * grad_alpha;
%     end
%     % 更新权矩阵Q、R、反馈系数H、控制向量d
% %     Q=q*eye(P);                 
%     R=r*eye(M);                 
%     H=[1;alpha*ones(N-1,1)];    
%     d=c*(inv(A'*Q*A+R)*A'*Q);   
end
% disp('Optimized Parameters:')
% % disp(['q: ', num2str(q)]);
% disp(['r: ', num2str(r)]);
% disp(['alpha: ', num2str(alpha)]);

%% 画图
plot(time(1:T_final),Yk(1:T_final));                       %纵坐标是Yk
hold on     
grid on
xlabel('time(s)'); 
ylabel('{\itv}(m/s)');
axis([0,T_final*Ts,0,1.5]);

x_limits = get(gca, 'XLim');                               % 获取当前坐标轴的 x 范围
line(x_limits, [1, 1], 'Color', 'red', 'LineStyle', '--'); % 画 y=1 的直线，设置颜色为红色，线型为虚线
legend('{\ity}','{\itw}', 'Location', 'NorthEast');
title('electromagnetic braking');
% title('Dynamic Matrix Control');
hold off
%% 梯度下降法附加函数 

% % 数值微分函数，用于近似计算参数的梯度
% function grad = numerical_gradient(param, objective)
%     epsilon = 1e-6; % 微小增量
%     %n = numel(param); % 参数个数
%     
%     %grad = zeros(size(param));
%     grad=0;
% %     for i = 1:n
% %         % 微小变化量
% %         delta = zeros(size(param));
% %         delta(i) = epsilon;
%         
%         % 增量方式计算目标函数值
% %         cost_plus = objective(param + delta);
% %         cost_minus = objective(param - delta);
%     cost_plus = objective(param + epsilon);
%     cost_minus = objective(param - epsilon);
%     % 近似计算梯度
%     grad = (cost_plus - cost_minus) / (2 * epsilon);
% %     end
% end