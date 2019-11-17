% ��Ⱥ�㷨�Ļ���ԭ��:
% 1��������·�����ͷ���Ϣ�ء�
% 2��������û�߹���·�ڣ��������ѡһ��·�ߡ�ͬʱ���ͷ���·�������йص���Ϣ�ء�
% 3����Ϣ��Ũ����·�����ȳɷ��ȡ������������ٴ�������·��ʱ����ѡ����Ϣ��Ũ�Ƚϸ�·����
% 4������·���ϵ���Ϣ��Ũ��Խ��Խ��
% 5��������Ⱥ�ҵ�����Ѱʳ·����

clear;
close all;
clc;

tic;

% ��Ⱥ����
Ant = 100;
% ��������
Ger = 50;
% ����������Χ������
xmin = -1;
xmax = 1;
ymin = -1;
ymax = 1;
tcl = 0.01;

% Ŀ�꺯��
f = 'cos(2*pi*x).*cos(2*pi*y).*exp(-(x.^2+y.^2)/10)';
[x,y] = meshgrid(xmin:tcl:xmax, ymin:tcl:ymax);
vxp = x;
vyp = y;
vzp = eval(f);

figure(1)
mesh(vxp,vyp,vzp);
hold on;
% ��ʼ������λ��
X = zeros(Ant,2);
T0 = zeros(Ant,1);
fitness = @(x)cos(2*pi*x(:,1)).*cos(2*pi*x(:,2)).*exp(-(x(:,1).^2+x(:,2).^2)/10);
for i = 1:Ant
    % ������������ֲ�����
    X(i,1) = (xmin + (xmax-xmin)*rand(1));
    X(i,2) = (ymin + (ymax-ymin)*rand(1));    
    % T0����Ϣ�أ�����ֵԽ����Ϣ��Ũ��Խ��
    % T0(i) = cos(2*pi*X(i,1)).*cos(2*pi*X(i,2)).*exp(-(X(i,1).^2+X(i,2).^2)/10);
    % ������Ϣ��Ũ��
    T0(i) = fitness(X(i,:));
end
plot3(X(:,1),X(:,2),T0,'k*');
hold on;grid on;
title('���ϵĳ�ʼ�ֲ�λ��');
xlabel('x');ylabel('y');zlabel('f(x,y)');

% ��ʼ�Ż�
T_Best = zeros(Ger,1);
Prob = zeros(Ger,Ant);
max_local = zeros(Ger,1);
max_global = zeros(Ger,1);
for i_ger = 1:Ger
    % ȫ��ת��ѡ������
    P0 = 0.2;
    % ��Ϣ������ϵ��
    P = 0.8;
    % ת�Ʋ�������
    lamda = 1/i_ger;   
    % ��ȡ�����Ϣ�أ���������λ��
    [T_Best(i_ger), BestIndex] = max(T0);   
    % ��ȡȫ��ת�Ƹ���
    for j_g = 1:Ant
        r = T0(BestIndex) - T0(j_g);
        Prob(i_ger, j_g) = r/T0(BestIndex);
    end
   
    for j_g_tr = 1:Ant
        if Prob(i_ger,j_g_tr)<P0
            % �ֲ�����
            temp(:,1) = X(j_g_tr,1) + (2*rand(1)-1)*lamda;
            temp(:,2) = X(j_g_tr,2) + (2*rand(1)-1)*lamda;
        else
            % ȫ������
            temp(:,1) = X(j_g_tr,1) + (xmax-xmin)*(rand(1)-0.5);
            temp(:,2) = X(j_g_tr,2) + (xmax-xmin)*(rand(1)-0.5);            
        end
        % ��ֹ����Խ�磬�ڿ������ڽ�������
        if temp(:,1)<xmin
            temp(:,1) = xmin;
        end
        if temp(:,1)>xmax
            temp(:,1) = xmax;
        end
        if temp(:,2)<ymin
            temp(:,2) = ymin;
        end
        if temp(:,2)>ymax
            temp(:,2) = ymax;
        end
        % ����λ�ú���ֵ��ԭλ�ú���ֵ���бȽϣ��ж��Ƿ����λ�ý��и���
        % if cos(2*pi*temp1).*cos(2*pi*temp2)*exp(-(temp1.^2+temp2.^2)/10) > cos(2*pi*X(j_g_tr,1)).*cos(2*pi*X(j_g_tr,2))*exp(-(X(j_g_tr,1).^2+X(j_g_tr,2).^2)/10)
        if fitness(temp) > fitness(X(j_g_tr,:)) 
            X(j_g_tr,1) = temp(:,1);
            X(j_g_tr,2) = temp(:,2);
        end
    end
    
    % ��Ϣ�ظ��£�������������Ϣ��+��ǰ��Ϣ��
    % ��Ϣ�ظ���ģ�ͣ�����ģ�ͣ�Ant-Cycleģ�ͣ�
    for t_t = 1:Ant
        % T0(t_t) = (1-P)*T0(t_t) + cos(2*pi*X(t_t,1)).*cos(2*pi*X(t_t,2)).*exp(-(X(t_t,1).^2+X(t_t,2).^2)/10);
        T0(t_t) = (1-P)*T0(t_t) + fitness(X(t_t,:));
    end
    
    % �ҳ������Ϣ��Ũ�ȣ���������Ӧ�ĸ���
    [c_iter, i_iter] = max(T0);
    maxpoint_iter = [X(i_iter,1), X(i_iter,2)];
    % max_local(i_ger) = cos(2*pi*X(i_iter,1)).*cos(2*pi*X(i_iter,2)).*exp(-(X(i_iter,1).^2+X(i_iter,2).^2)/10);
    max_local(i_ger) = fitness(X(i_iter,:));
    
    % ��ÿ��ȫ�����Ž�洢��max_global������
    if i_ger>=2
        if max_local(i_ger) > max_global(i_ger-1)
            max_global(i_ger) = max_local(i_ger);
        else
            max_global(i_ger) = max_global(i_ger-1);
        end
    else
        max_global(i_ger) = max_local(i_ger);
    end
end

%%
figure(2)
mesh(vxp,vyp,vzp);
hold on;
x = X(:,1);
y = X(:,2);
plot3(x,y,eval(f),'b*')
hold on;
grid on;
title('���ϵ����շֲ�λ��')
xlabel('x');
ylabel('y');
zlabel('f(x,y)');

figure(3)
plot(1:Ger, max_global, 'b-')
hold on;
title('���ź����仯����');
xlabel('iteration');
ylabel('f(x)');
grid on;
[c_max, i_max] = max(T0);
Maxpoint = [X(i_max,1),X(i_max,2)]
% maxvalue = cos(2*pi*X(i_max,1)).*cos(2*pi*X(i_max,2)).*exp(-(X(i_max,1).^2+X(i_max,2).^2)/10)
maxvalue = fitness(X(i_max,:))









