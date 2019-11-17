% 曲线与曲面绘制
clc;
clear;
close all;

% 向量曲线绘制
syms t
r = {cos(t), sin(t), t};
% class(r(1))
figure(1)
subplot(1,3,1); fplot(r{1},[0 5*pi],'r--');axis equal
title('余弦函数'),xlabel('t'),ylabel(sprintf('%s',r{1})); %普通函数
subplot(1,3,2); fplot(r{1},r{2},[0 5*pi],'b.-');axis equal
title('圆函数'),xlabel(sprintf('%s',r{1})),ylabel(sprintf('%s',r{2})); %参数函数
subplot(1,3,3); fplot3(r{1},r{2},r{3},[0 5*pi],'k-.');axis equal
title('圆柱螺旋函数'),xlabel(sprintf('%s',r{1})),ylabel(sprintf('%s',r{2})),zlabel(sprintf('%s',r{3}));%空间参数函数

% 向量曲面绘制（meshgrid、mesh、surf）
syms u v
% r = {cos(u)*cos(v), cos(u)*sin(v), sin(u)};
figure(2)
u = -pi/2:0.1:pi/2;
v = 0:0.1:2*pi;
[u,v] = meshgrid(u,v); %生成格点矩阵
x = cos(u).*cos(v);
y = cos(u).*sin(v);
z = sin(u);
subplot(2,2,1)
mesh(x,y,z) %绘制曲面网格图（球面）
axis equal

u = 1:0.1:5;
v = 0:0.1:2*pi;
[u,v] = meshgrid(u,v); %生成格点矩阵
x = sin(u).*cos(v);
y = sin(u).*sin(v);
z = v;
subplot(2,2,2)
mesh(x,y,z) %绘制曲面网格图（正螺旋面）
axis equal

u = 1:0.1:10;
v = 0:0.1:2*pi;
[u,v] = meshgrid(u,v); %生成格点矩阵
x = sin(u).*cos(v);
y = sin(u).*sin(v);
z = log(abs(tan(u/2))) + cos(u);
subplot(2,2,3)
mesh(x,y,z) %绘制曲面网格图（伪球面）
axis equal

t = 0:0.1:2.01*pi;
v = -1:0.1:1;
[v,t] = meshgrid(v,t); %生成格点矩阵
r = 2+0.5*v.*cos(t/2);
x = r.*cos(t);
y = r.*sin(t);
z = 0.5*v.*sin(t/2);
subplot(2,2,4)
mesh(x,y,z) %绘制曲面网格图（麦比乌斯曲面）
axis equal

% close Figure 1 2

% 曲率和挠率计算
syms t real
r = [sin(t), cos(t), t];
r1 = diff(r);
r2 = diff(r1); % r2 = diff(r,2);
m = cross(r1,r2); % 计算外积
% 计算曲率
k = simplify((sqrt(m(1)^2 + m(2)^2 + m(3)^2))) / simplify(sqrt(r1(1)^2 + r1(2)^2 + r1(3)^2))^3;
sprintf('曲率：%s',k);

r3 = diff(r2);
C = [r1;r2;r3];
% 计算挠率
tt = simplify(det(C) / sqrt(m(1)^2 + m(2)^2 + m(3)^2)^2);
sprintf('挠率：%s',tt);

% 计算弧长
s = int(norm(r1),t,0,2*pi);
sprintf('弧长：%s',s);

% 向量曲面计算
% 限定符号变量类型 real | positive | integer | rational
syms u v R real
r = {R*cos(u)*cos(v), R*cos(u)*sin(v), R*sin(u)}; %u = theta; v = phi
% r = {cos(u),sin(u),v};
ru = diff(r,u);
rv = diff(r,v);
simplify(cross(ru,rv));

% 曲面的第一基本形式：dr = ru*du + rv*dv
% conj(Z) = real(Z) - i*imag(Z)
E = simplify(dot(ru,ru));
F = simplify(dot(ru,rv));
G = simplify(dot(rv,rv));
% class(G)
G2 = string(G);
% strfind：查找字符串中的特定子串，并返回子串第一个字符的位置
% strrep：将字符串中特定子串替换为另一子串
% strtok：以空格为标志提取字符串中的子串
s_c = {"sin(u)^2 - 1", "1 - sin(u)^2", "cos(u)^2 - 1", "1 - cos(u)^2"};
s_c_ch = {"-cos(u)^2", "cos(u)^2", "-sin(u)^2", "sin(u)^2"};
class(s_c{1});
% 判断字符串数组中的字符串是否存在于被查找的字符串中，并对其做相应的替换
for i=1:length(s_c)
    if contains(G2, s_c{i})
        G = strrep(G2, s_c{i}, s_c_ch{i});
        G = simplify(str2sym(G));
    end
end

syms du2 du dv dv2
I = E*du2 + 2*F*du*dv + G*dv2;

% 曲面上u曲线与v曲线夹角
theta = acos(F/sqrt(E*G));

% 曲面域D的面积（二重积分计算）
area = int(int(sqrt(E*G-F),u,-pi/2,pi/2),v,0,2*pi);

% 曲面第二基本形式
ruu = diff(ru,u);
ruv = diff(ru,v);
rvv = diff(rv,v);

% 消去绝对值符号
P = simplify(sqrt(E*G-F^2));
P2 = string(P);
if contains(P2, 'abs')
    P2 = strrep(P2,'abs','');
end
P = str2sym(P2);

L = simplify(dot(ruu, cross(ru,rv)))/P;
M = simplify(dot(ruv, cross(ru,rv)))/P;
N = simplify(dot(rvv, cross(ru,rv)))/P;
II = L*du2 + 2*M*du*dv + N*dv2;

% 法截线曲率
k = II/I;

% 欧拉公式：kn = k1*cos(theta)^2 + k2*sin(theta)^2
% 高斯曲率（全曲率、总曲率）：K = k1*k2
K = (L*N-M^2)/(E*G-F^2);
% 平均曲率：H = (k1 + k2)/2
H = (L*G-2*M*F+N*E)/(2*(E*G-F^2));










