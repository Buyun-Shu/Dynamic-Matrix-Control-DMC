%%��ɢ��̬������ƣ�DMC��
clc;        %��������
clear all;  %�������������
close all;  %�ر�����ͼ��
%% ϵͳģ�ͽ���
% num=5;                          %���Ӷ���ʽ
% den=[50 5];                     %��ĸ����ʽ
num=0.115;
den=[100 0.603];
% num=[10256.86 0];
% den=[1 0 35.2836];
                                %den=[1 0.5 0.1];������ϵͳ�ķ�ĸ����ʽ��
sys=tf(num,den,'InputDelay',1); %�������ݺ���ϵͳģ�ͣ������ӳ�Ϊ1����λ
Ts=1;                           %��ɢ������ʱ������1����λ���������ڣ�
dsys=c2d(sys,Ts,'foh');         %c2d������ʱ��ϵͳsys����ǰ��ŷ�����Ͳ���ʱ��ת��Ϊ��ɢʱ��ϵͳdsys
[num1,den1]=tfdata(dsys,'v');   %��ȡ��ɢʱ��ϵͳ�ķ��ӡ���ĸϵ��
% %step(dsys);
% grid on
% %T_final=200;
T_final=4000;
% tt=Ts:Ts:T_final*Ts;            %����tt������Ts��T_final*Ts�ȼ��ʱ���
% a1=step(dsys,tt);               %�洢��ɢʱ��ϵͳdsys��ʱ������tt�ϵĲ���Ӧ����
loaded_data=load('sudu.mat','sudu');
a1 = loaded_data.sudu;
%% ģ��ʵʱ�ɼ��ܾ�����
% n = T_final;                    %ָ�����鳤��
% rangeMin = 180;                 %�ܾ�����mm
% rangeMax = 220;                 %�ܾ�����mm
% real_time_data = (rangeMax-rangeMin).*rand(1, n) + rangeMin;
% 
% shape_parameter = 1;            %��״����
% sensitivity = 10;               %���жȲ���
% sum_rk = real_time_data(1);
%% ���ÿ����������������޸�֮����
N=20;         %��ģʱ��Ĭ��NT֮���Ծ��Ӧ����ƽ��(ԭ��N=20��
                %��ģʱ��N��������ȷ��������"ʹ�ö��ٸ������ʷ������������������Ԥ��"
                %�ϴ�Nֵ���ṩ����ʱ�䷶Χ��Ԥ�⣬�������Ӽ��㸴�ӶȺ��ӳ�
                %NԽ������Խ�⻬

P=6;           %�Ż�ʱ�򣬱��ض����ڿ�������������δ��P��ʱ��Ԥ�����ֵ***P=6
                %����ȷ�����������ɵ����ſ������еĳ���
                %�ϴ�Pֵ���Ը��õش���ϵͳ�ĳ�ʱ��Ӱ��ͱ仯���������Ӽ��㸴�Ӷ�

M=2;           %����ʱ��ÿһʱ��k���M��������������***M=2     P���ڵ���M
                %��ʾ��ÿ������ʱ����ǰԤ������źŵ�ʱ�䲽��
                %�ϴ��Mֵʹ�������ܹ�Ԥ��δ������ʱ�䲽�Ŀ����źţ���Ҳ�����Ӽ��㸴�ӶȺ��ӳ�

q=1;            %Ȩϵ�������Ƹ������仯***q=1
                %�ϴ��qֵ��ʾ�����Ӽ�С�����ʹ�ÿ�����������ڿ����ź���׷�������ο��켣�������������Ŀ��Ʋ�����
                %��С��qֵ����ζ�Ÿ����ݵ�ϵͳ��������������ܸ������ڱ����ȶ�����С���Ʋ���
     
r=0.1;          %Ȩϵ�������ƿ������仯***r=0.1
                %�ϴ�rֵ�����Ӽ�С��������ķ��ȣ�ʹ�ÿ��������������ڿ����ź��Լ�С���ƻ����������Ӧ�������ȶ���׼ȷ��
                %��Сrֵ����ζ�Ÿ����ɵĿ����������ƣ����������ܸ������ص��ڿ����ź��Կ�����Ӧϵͳ�����󣨿�����Ӧ��

%w=2.5;         %����ֵ***
w=1;
a=a1(1:N);      %ģ�Ͳ���������λ��Ծ��Ӧ����ֵ

alpha=0.8;      %�ữϵ��***alpha=0.8
                %�ữϵ��ֵԽ�󣬿����������Խƽ������������ữϵ��ʹ�ÿ���ϵͳ��Ӧ��óٶ�
                %0.9ʱ����û�г�����������ʱ���Լ70s
                %0.2ǰ�β��⻬���нϴ󳬵���������ʱ��Ҳ�ϳ�
                %0.1ǰ�β��⻬���޳�����������ʱ���Լ60s
%% ���㶯̬����A
for i=1:P                   %A��һ��P��M�еľ���
    for j=1:M
        if j>i
            A(i,j)=0;
        else
            A(i,j)=a(i-j+1);%����ģ�Ͳ�������ֵ�����A����
        end
    end
end
%% Ȩ����Q��R������ϵ��H����������d
Q=q*eye(P);                 %���Ȩ����P��P�еĶԽǾ���
R=r*eye(M);                 %����Ȩ����M��M�еĶԽǾ���
H=[1;alpha*ones(N-1,1)];    %У��������������δ�������Ԥ�⣬N��1�е���������������Ԫ����1�����඼��alpha
c=[1,zeros(1,M-1)];         %Mά��������ȡ��Ԫ�ؼ���
d=c*(inv(A'*Q*A+R)*A'*Q);   %��������d��һ��1��M�е�����
%% ������λ����S
S=zeros(N);                 %S��һ��N��N�еľ������ڹ���DMC�㷨���ӳ��������λ����
for i=1:N
    if i<N
        S(i,i+1)=1;
    end
    if i==N
        S(N,N)=1;
    end
end
%% ��ʼ��DMC
Alpha=[];                   %��ʼ��Alpha����
for j=1:P
    Alpha=[Alpha,alpha^j];  %�洢�ữϵ����1~P�ݴη�
end
Uk=zeros(1,M);              %�������루�����Ż����ʵ�����룩
Yp0k=[];                    %Ԥ���k+1��k+Pδ��������������仯0��
for i=1:P
    Yp0k=[Yp0k;0];          
end
Yn0k=[];                    %�ӳ���������
for i=1:N
    Yn0k=[Yn0k;0];          %N��1�е�������        
end
Yk=zeros(1,N);              %1��N�е������
uk=0;                       %ʵ�ʿ�������
yk=0;                       %��ǰ״̬
utemp=0;                    %��ʱ���������ڴ洢ʵ�ʿ�������uk�ĵ�ǰ����ֵ
%% �ݶ��½�����ʼ��

% ��ʼ��ѧϰ��
% learning_rate = 0.01;       % �ɸ���ʵ��������е���
% max_iterations = 100;       % ����������
%% DMCʵʱ����
for i=1:T_final
    time(i)=i*Ts;    
    y1k=-den1(2)*yk+num1(2)*uk;     %ʵ�����y1k����ǰ״̬yk ʵ�ʿ�������uk������
    %y2k=-den1(3)*yk-den1(2)*y1k+num1(3)uk;�����������ڶ�������ģ�ͣ�
    
    %����������Ӧ
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

    Yk(i)=y1k;                      %��y1k�洢��Yk������
    Uk(i)=utemp;                    %����ǰʵ�ʿ�������uk�ĸ���ֵ�洢��Uk������
    Yr=(1-Alpha')*w+Alpha'*y1k;     %�������Ԥ��ֵ����ο��켣��
    error=y1k-Yp0k(1);              %�����ʵ�����y1k��-ģ��Ԥ���������Yp0k��Ԫ��   ��3-10��
    Ycork=Yn0k+H*error;             %YcorkΪУ��������Ԥ��������HΪУ������            ��3-11��
    Yn0k=S*Ycork;                   %�ӳ�����Yn0k��SΪ��λ����                          ��3-13��
    Yp0k=Yn0k(1:P);                 %���� Ԥ���������Yp0k����ȡYn0k��ǰP��Ԫ��
    dU=d*(Yr-Yp0k);                 %������������dU��dΪ��������                         ��3-7��

    utemp=uk+dU;                    %������ʱ����utemp���õ� �µ�ʵ�ʿ�������uk��ֵ
    Yp1k=Yp0k+a(1:P)*dU;            %�����µ�Ԥ���������Yp1k
    Yn1k=Yn0k+a(1:N)*dU;            %�����µ��ӳ���������Yn1k
    Yp0k=Yp1k;                      %����Ԥ���������Yp0k����Yp1k��ֵ����Yp0k
    Yn0k=Yn1k;                      %�����ӳ���������Yn0k����Yn1k��ֵ����Yn0k
    uk=utemp;                       %����ʵ�ʿ�������uk��ֵ
    yk=y1k;                         %�������yk��ֵ
    
%     if i>=2
%         delta_rk = real_time_data(i) - real_time_data(i-1);
%         sum_rk = sum_rk + real_time_data(i);
%         equal_rk = sum_rk/i;
%         weighting_function = 1/(1+shape_parameter*((abs(delta_rk/equal_rk))^sensitivity));
%         r = r * weighting_function;
%     end
%     %error1=Yk-(Yn0k+([1;alpha*ones(N-1,1)]).*error);        
%     objective = @(alpha) sum((Yk' - (Yn0k + ([1; alpha.*ones(N-1,1)]).*error)).^2)*(1/i);%�Զ����Ŀ�꺯����������������Ч��
%    
%     %Yk��1��N�е������Yn0k��N��1�е���������H��N��1�е�������
% 
%     % �ݶ��½����Ż�����
%     for iter = 1:max_iterations
%         % ʹ����ֵ΢�ֽ��Ƽ����ݶ�
% %         grad_Q = numerical_gradient(q, error1);
% %         grad_R = numerical_gradient(r, error1);
%         grad_alpha = numerical_gradient(alpha, objective);
%         
%         % ���¿���������
% %         q = q - learning_rate * grad_Q;
% %         r = r - learning_rate * grad_R;
%         alpha = alpha - learning_rate * grad_alpha;
%     end
%     % ����Ȩ����Q��R������ϵ��H����������d
% %     Q=q*eye(P);                 
%     R=r*eye(M);                 
%     H=[1;alpha*ones(N-1,1)];    
%     d=c*(inv(A'*Q*A+R)*A'*Q);   
end
% disp('Optimized Parameters:')
% % disp(['q: ', num2str(q)]);
% disp(['r: ', num2str(r)]);
% disp(['alpha: ', num2str(alpha)]);

%% ��ͼ
plot(time(1:T_final),Yk(1:T_final));                       %��������Yk
hold on     
grid on
xlabel('time(s)'); 
ylabel('{\itv}(m/s)');
axis([0,T_final*Ts,0,1.5]);

x_limits = get(gca, 'XLim');                               % ��ȡ��ǰ������� x ��Χ
line(x_limits, [1, 1], 'Color', 'red', 'LineStyle', '--'); % �� y=1 ��ֱ�ߣ�������ɫΪ��ɫ������Ϊ����
legend('{\ity}','{\itw}', 'Location', 'NorthEast');
title('electromagnetic braking');
% title('Dynamic Matrix Control');
hold off
%% �ݶ��½������Ӻ��� 

% % ��ֵ΢�ֺ��������ڽ��Ƽ���������ݶ�
% function grad = numerical_gradient(param, objective)
%     epsilon = 1e-6; % ΢С����
%     %n = numel(param); % ��������
%     
%     %grad = zeros(size(param));
%     grad=0;
% %     for i = 1:n
% %         % ΢С�仯��
% %         delta = zeros(size(param));
% %         delta(i) = epsilon;
%         
%         % ������ʽ����Ŀ�꺯��ֵ
% %         cost_plus = objective(param + delta);
% %         cost_minus = objective(param - delta);
%     cost_plus = objective(param + epsilon);
%     cost_minus = objective(param - epsilon);
%     % ���Ƽ����ݶ�
%     grad = (cost_plus - cost_minus) / (2 * epsilon);
% %     end
% end