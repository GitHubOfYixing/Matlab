% 遗传算法
% clc;clear;close all;

% umin = -3;
% umax = 3;
% delta = 0.0001;
% 确定变量的取值范围:(-3,3); 编码精度: 0.0001
% 设置交配概率Pc：通常取0.4-0.9之间
% Pc = 0.7;
% 设置变异概率Pm：通常取0.001-0.1之间
% Pm = 0.1;
% 迭代次数
% N_r = 100;

% x = -2.5:0.0001:2.5;
% y = (x - 4*x.^2 + x.^4);
% plot(x,y,'-r');

% f = 'x - 8*x.^2 + 3*x.^4';
% fitness = '20 - (x - 8*x.^2 + 3*x.^4)';
f = '5*x.*exp(x) + cos(x.^2./6)';
fitness = '20 - (5*x.*exp(x) + cos(x.^2./6))';
N = 20;
M = zeros(1,N);
for i=1:N
    [g_min,Sdecu,J_min] = Genetic_Algorithm(f,fitness,-100,0,0.0001,0.8,0.01,100);
    M(i) = J_min;
end
% 查找出现次数最多的那个数，均值，中位数
sprintf('真值: %.3f\n众数: %.3f\n均值: %.3f\n中数: %.3f\n大值: %.3f\n小值: %.3f\n',g_min, mode(M), mean(M), median(M), max(M) ,min(M))


