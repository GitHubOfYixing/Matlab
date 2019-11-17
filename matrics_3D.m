% 矩阵还原
clc;
clear;
close all;

A = [1 2;3 4];
% A = [1 2 3;4 5 6;7 8 9];  % 保持中心位置不变，边块元素始终在边块位置，角块元素始终在角块位置
% A = [1 2 3 4;5 6 7 8;9 10 11 12;13 14 15 16];

% 1. 查找中心元素
n = length(A);
% 取余运算：rem(a,b)；取模运算：mod(a,b)；向下取整：floor(a/b)；向上取整：ceil(a/b)；最接近取整：round(a/b)
if ~rem(n,2)
    P0 = [A(n/2,n/2), A(n/2,n/2+1), A(n/2+1,n/2), A(n/2+1,n/2+1)]; % 提取中心元素
    r = randperm(length(P0)); % 生成关于列数的随机排列列数序列
    P0 = P0(r); % 随机打乱数组元素顺序
    A(n/2,n/2) = P0(1);
    A(n/2,n/2+1) = P0(2);
    A(n/2+1,n/2) = P0(3);
    A(n/2+1,n/2+1) = P0(4);
else
    A0c = A(ceil(n/2),ceil(n/2));
end

% 2. 查找边界元素（由内向外依次查找）
if n >= 2 % 3阶及以上才存在边块
    P1 = [A(1,2:n-1) A(n,2:n-1) A(2:n-1,1)' A(2:n-1,n)'];
    r = randperm(length(P1)); % 生成关于列数的随机排列列数序列
    P1 = P1(r); % 随机打乱数组元素顺序
    A(1,2:n-1) = P1(1:n-2);
    A(n,2:n-1) = P1(n-1:2*n-4);
    A(2:n-1,1) = P1(2*n-3:3*n-6)';
    A(2:n-1,n) = P1(3*n-5:4*n-8)';
end

% 3. 查找角点元素
if n >= 2
    P2 = [A(1,1) A(1,n) A(n,1) A(n,n)];
    r = randperm(length(P2)); % 生成关于列数的随机排列列数序列
    P2 = P2(r); % 随机打乱数组元素顺序
    A(1,1) = P2(1);
    A(1,n) = P2(2);
    A(n,1) = P2(3);
    A(n,n) = P2(4);
end



