% 矩阵行列变换
clc;
clear;
close all;

format short % 小数显示
format rat % 分数显示
A = [1 2;4 6];

A1 = [1 6;4 2];
A2 = [1 2;6 4];
A3 = [4 2;1 6];
A4 = [2 1;4 6];
Ac1 = A1\A;
Ac2 = A2\A;
Ac3 = A3\A;
Ac4 = A4\A;

B1 = [6 2;4 1];
B2 = [1 4;2 6];
Bc1 = B1\A;
Bc2 = B2\A;

C1 = [4 6;1 2];
C2 = [2 1;6 4];
Cc1 = C1\A;
Cc2 = C2\A;

D1 = [2 6;1 4];
D2 = [6 4;2 1];
D3 = [4 1;6 2];
D4 = [1 2;4 6];
Dc1 = D1\A;
Dc2 = D2\A;
Dc3 = D3\A;
Dc4 = D4\A;



