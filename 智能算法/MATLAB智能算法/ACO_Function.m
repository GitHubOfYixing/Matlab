% 蚁群算法的基本原理:
% 1、蚂蚁在路径上释放信息素。
% 2、碰到还没走过的路口，就随机挑选一条路走。同时，释放与路径长度有关的信息素。
% 3、信息素浓度与路径长度成反比。后来的蚂蚁再次碰到该路口时，就选择信息素浓度较高路径。
% 4、最优路径上的信息素浓度越来越大。
% 5、最终蚁群找到最优寻食路径。

clear;
close all;
clc;

tic;

% 种群数量
Ant = 100;
% 迭代次数
Ger = 50;
% 设置搜索范围及步长
xmin = -1;
xmax = 1;
ymin = -1;
ymax = 1;
tcl = 0.01;

% 目标函数
f = 'cos(2*pi*x).*cos(2*pi*y).*exp(-(x.^2+y.^2)/10)';
[x,y] = meshgrid(xmin:tcl:xmax, ymin:tcl:ymax);
vxp = x;
vyp = y;
vzp = eval(f);

figure(1)
mesh(vxp,vyp,vzp);
hold on;
% 初始化蚂蚁位置
X = zeros(Ant,2);
T0 = zeros(Ant,1);
fitness = @(x)cos(2*pi*x(:,1)).*cos(2*pi*x(:,2)).*exp(-(x(:,1).^2+x(:,2).^2)/10);
for i = 1:Ant
    % 在区间内随机分布蚂蚁
    X(i,1) = (xmin + (xmax-xmin)*rand(1));
    X(i,2) = (ymin + (ymax-ymin)*rand(1));    
    % T0：信息素，函数值越大，信息素浓度越大
    % T0(i) = cos(2*pi*X(i,1)).*cos(2*pi*X(i,2)).*exp(-(X(i,1).^2+X(i,2).^2)/10);
    % 计算信息素浓度
    T0(i) = fitness(X(i,:));
end
plot3(X(:,1),X(:,2),T0,'k*');
hold on;grid on;
title('蚂蚁的初始分布位置');
xlabel('x');ylabel('y');zlabel('f(x,y)');

% 开始优化
T_Best = zeros(Ger,1);
Prob = zeros(Ger,Ant);
max_local = zeros(Ger,1);
max_global = zeros(Ger,1);
for i_ger = 1:Ger
    % 全局转移选择因子
    P0 = 0.2;
    % 信息素蒸发系数
    P = 0.8;
    % 转移步长参数
    lamda = 1/i_ger;   
    % 获取最大信息素，并保留其位置
    [T_Best(i_ger), BestIndex] = max(T0);   
    % 求取全局转移概率
    for j_g = 1:Ant
        r = T0(BestIndex) - T0(j_g);
        Prob(i_ger, j_g) = r/T0(BestIndex);
    end
   
    for j_g_tr = 1:Ant
        if Prob(i_ger,j_g_tr)<P0
            % 局部搜索
            temp(:,1) = X(j_g_tr,1) + (2*rand(1)-1)*lamda;
            temp(:,2) = X(j_g_tr,2) + (2*rand(1)-1)*lamda;
        else
            % 全局搜索
            temp(:,1) = X(j_g_tr,1) + (xmax-xmin)*(rand(1)-0.5);
            temp(:,2) = X(j_g_tr,2) + (xmax-xmin)*(rand(1)-0.5);            
        end
        % 防止搜索越界，在可行域内进行搜索
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
        % 将新位置函数值与原位置函数值进行比较，判断是否对其位置进行更新
        % if cos(2*pi*temp1).*cos(2*pi*temp2)*exp(-(temp1.^2+temp2.^2)/10) > cos(2*pi*X(j_g_tr,1)).*cos(2*pi*X(j_g_tr,2))*exp(-(X(j_g_tr,1).^2+X(j_g_tr,2).^2)/10)
        if fitness(temp) > fitness(X(j_g_tr,:)) 
            X(j_g_tr,1) = temp(:,1);
            X(j_g_tr,2) = temp(:,2);
        end
    end
    
    % 信息素更新：保留下来的信息素+当前信息素
    % 信息素更新模型：蚁周模型（Ant-Cycle模型）
    for t_t = 1:Ant
        % T0(t_t) = (1-P)*T0(t_t) + cos(2*pi*X(t_t,1)).*cos(2*pi*X(t_t,2)).*exp(-(X(t_t,1).^2+X(t_t,2).^2)/10);
        T0(t_t) = (1-P)*T0(t_t) + fitness(X(t_t,:));
    end
    
    % 找出最大信息素浓度，及其所对应的个体
    [c_iter, i_iter] = max(T0);
    maxpoint_iter = [X(i_iter,1), X(i_iter,2)];
    % max_local(i_ger) = cos(2*pi*X(i_iter,1)).*cos(2*pi*X(i_iter,2)).*exp(-(X(i_iter,1).^2+X(i_iter,2).^2)/10);
    max_local(i_ger) = fitness(X(i_iter,:));
    
    % 将每代全局最优解存储到max_global矩阵中
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
title('蚂蚁的最终分布位置')
xlabel('x');
ylabel('y');
zlabel('f(x,y)');

figure(3)
plot(1:Ger, max_global, 'b-')
hold on;
title('最优函数变化趋势');
xlabel('iteration');
ylabel('f(x)');
grid on;
[c_max, i_max] = max(T0);
Maxpoint = [X(i_max,1),X(i_max,2)]
% maxvalue = cos(2*pi*X(i_max,1)).*cos(2*pi*X(i_max,2)).*exp(-(X(i_max,1).^2+X(i_max,2).^2)/10)
maxvalue = fitness(X(i_max,:))









