% 广义逆矩阵求解
clc;clear;close all;

% A = [2 4 1 1;1 2 -1 2;-1 -2 -2 1];
% b = [3 0 3]';
% A = [0 1 -1;-1 0 1;1 -1 0;0 1 -1];
% b = [1 -1 0 1]';
A = [1 1 0;1 1 2;0 0 2];
b = [1 1 1]';

% 对矩阵A进行最大秩分解
R = rref(A);    %计算行最简式
[m, n] = size(A);
B(1:m,:) = 0;
D(:,1:n) = 0;
for i = 1:m
    for j = 1:n
        if(R(i,j)==1.0 && sum(R(1:i-1,j))==0 && sum(R(i+1:m,j))==0)
            B = [B,A(:,j)];
            D = [D;R(i,:)];
        end
    end
end
B = B(:,2:end)
D = D(2:end,:)

% B = [1 0;1 2;0 2];
% D = [1 1 0;0 0 1];

% 计算A+
% 显示分数解
format rat
Aplus = D'*inv(D*D')*inv(B'*B)*B'

% 判断A*Aplus*b = b
b2 = A*Aplus*b;
if isequal(b2, b)
    disp('Ax=b是相容方程组，存在最小范数解')
else
    disp('Ax=b是不相容方程组，存在最小二乘解（最佳逼近解）')
end

% 求最小二乘解的通解及最佳逼近解
syms L E u
x = (Aplus*b)'

% A = [1 2 -1;0 -1 2];
% Aplus = (1/14)*[5 4;6 2;3 8];
% 
% R1 = Aplus*A
% R2 = Aplus*A*Aplus*A
% 
% format short
% R1-R2
% norm(R1-R2)







