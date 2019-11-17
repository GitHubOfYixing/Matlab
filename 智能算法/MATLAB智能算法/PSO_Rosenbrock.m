% 求解无约束优化问题（Rosenbrock函数，维数为5）
% 
clc;
clear;
close all;

%**************************************************%
[x,y] = meshgrid(-3:0.01:3,-3:0.01:3);
z = 100*x.*y.*exp(-x.^2-y.^2-x.*y/5);
% z = (1-x).^2+100*(y-x.^2).^2;

% z1 = (1-x1).^2+100*(x2-x1.^2).^2;
% z = z1 + (1-x2).^2+100*(x3-x2.^2).^2;
% z3 = z2 + (1-x3).^2+100*(x4-x3.^2).^2;
% z4 = z3 + (1-x4).^2+100*(x5-x4.^2).^2;

figure(1);
mesh(x,y,z);%三维网格图
% meshc(x,y,z);%带等高线的三维网格图
% meshz(x,y,z);%带底座的三维网格图
% surf(x,y,z);%三维曲面图
title('z = 100*x*y*exp(-x^2-y^2-x*y/5)');
% title('z = (1-x).^2+100*(y-x.^2).^2');
xlabel('x');
ylabel('y');
zlabel('z');

% 迭代次数
Index = 1;
% 最大迭代次数
MaxIndex = 100;
% 认知系数
c1 = 0.6;
% 社会学习系数
c2 = 0.7;
% 惯性系数，随迭代次数增加而递减
wmax = 0.9;
wmin = 0.4;
w = wmax - (wmax-wmin)*(Index/MaxIndex)^2;
% 粒子最大飞行速度（-3,3）
vmax = 6;
% 适应度函数
fitness = @(x)100*x(:,1).*x(:,2).*exp(-x(:,1).^2-x(:,2).^2-x(:,1).*x(:,2)/5);
% fitness = @(x)(1-x(:,1)).^2+100*(x(:,2)-x(:,1).^2).^2;
% fitness = @(x)(1-x(:,1)).^2+100*(x(:,2)-x(:,1).^2).^2 + (1-x(:,2)).^2+100*(x(:,3)-x(:,2).^2).^2;
% 随机初始化粒子种群速度V与位置X
% 设置种群大小为20，每个粒子的维数为2维
N = 20; Dim = 2;
% 初始化每个粒子位置x
x = -3 + 6*rand(N,Dim);
% 初始化每个粒子速度v
v = vmax*rand(N,Dim);
% 计算每个个体的适应度
f = fitness(x);

% 设置初始个体历史最佳适应度
fpbest = f;
% 设置初始个体历史最佳位置
pbest = x;
% 设置初始种群历史最佳适应度，i为最佳适应度所对应的最佳位置
[fgbest,i] = min(fpbest);
% 设置初始种群历史最佳位置
gbest = pbest(i,:);
                    
% 开始迭代
figure(2)
t = zeros(1,MaxIndex);
y = zeros(1,MaxIndex);
while Index <= MaxIndex
    scatter(x(:,1),x(:,2),'yo'); % 二维动态显示粒子的更新过程
    title(sprintf('fmin(%.4f)',fgbest))
    axis([-3 3 -3 3]);

%     scatter3(x(:,1),x(:,2),f,'yo'); % 三维动态显示粒子的更新过程
%     axis([-3 3 -3 3 -20 20]);

    % 绘制收敛曲线
%     t(Index) = Index;
%     y(Index) = fgbest;
%     plot(t,y,'b','MarkerSize',5);
%     p = plot(t,y,'b','MarkerSize',5);
%     set(p,'XData',t,'YData',y)
%     drawnow
%     axis([0 MaxIndex -30 -10]);
    
    % 更新惯性权值（较大权值有利于全局搜索，较小权值有利于局部搜索）
    w = wmax - (wmax-wmin)*(Index/MaxIndex)^2;
    % 更新粒子
    % 由个体和种群的最佳信息数据进行更新移动速度（注意：这里需要将全局最优值进行升维处理，使其与x维度一致）
    v = w*v + c1*rand()*(pbest - x) + c2*rand()*(ones(N,1)*gbest - x); 
    % 通过速度更新位置
    x = x + v;
    
    % 评估每个粒子适应度
    f = fitness(x);
    % 查找是否存在更小的适应度值
    k = find(f(:) < fpbest(:));
    % 如果存在，则更新粒子个体
    if size(k)
        % 更新个体最佳适应度
        fpbest(k) = f(k);
        % 更新最优个体位置
        pbest(k,:) = x(k,:);
    end
    % 更新种群最佳适应度
    [fgbest,i] = min(fpbest); 
    % 更新种群最佳位置
    gbest = pbest(i,:);
    
    Index = Index + 1;
    pause(0.1)
end
text(gbest(:,1)-0.7,gbest(:,2)+0.3,sprintf('(%.4f, %.4f)',gbest(:,1),gbest(:,2)))
gbest
fgbest

% fmincon
fun = @(x)100*x(:,1).*x(:,2).*exp(-x(:,1).^2-x(:,2).^2-x(:,1).*x(:,2)/5);
% fun = @(x)(1-x(1)).^2+100*(x(2)-x(1).^2).^2;
% fun = @(x)(1-x(:,1)).^2+100*(x(:,2)-x(:,1).^2).^2 + (1-x(:,2)).^2+100*(x(:,3)-x(:,2).^2).^2;
x0 = [-2,2];
% A = [2,2];
% b = 1;
% x = fmincon(fun,x0,A,b)
xstd = fmincon(fun,x0)
fstd = fun(xstd)


