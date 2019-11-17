% �������������
clc;
clear;
close all;

% �������߻���
syms t
r = {cos(t), sin(t), t};
% class(r(1))
figure(1)
subplot(1,3,1); fplot(r{1},[0 5*pi],'r--');axis equal
title('���Һ���'),xlabel('t'),ylabel(sprintf('%s',r{1})); %��ͨ����
subplot(1,3,2); fplot(r{1},r{2},[0 5*pi],'b.-');axis equal
title('Բ����'),xlabel(sprintf('%s',r{1})),ylabel(sprintf('%s',r{2})); %��������
subplot(1,3,3); fplot3(r{1},r{2},r{3},[0 5*pi],'k-.');axis equal
title('Բ����������'),xlabel(sprintf('%s',r{1})),ylabel(sprintf('%s',r{2})),zlabel(sprintf('%s',r{3}));%�ռ��������

% ����������ƣ�meshgrid��mesh��surf��
syms u v
% r = {cos(u)*cos(v), cos(u)*sin(v), sin(u)};
figure(2)
u = -pi/2:0.1:pi/2;
v = 0:0.1:2*pi;
[u,v] = meshgrid(u,v); %���ɸ�����
x = cos(u).*cos(v);
y = cos(u).*sin(v);
z = sin(u);
subplot(2,2,1)
mesh(x,y,z) %������������ͼ�����棩
axis equal

u = 1:0.1:5;
v = 0:0.1:2*pi;
[u,v] = meshgrid(u,v); %���ɸ�����
x = sin(u).*cos(v);
y = sin(u).*sin(v);
z = v;
subplot(2,2,2)
mesh(x,y,z) %������������ͼ���������棩
axis equal

u = 1:0.1:10;
v = 0:0.1:2*pi;
[u,v] = meshgrid(u,v); %���ɸ�����
x = sin(u).*cos(v);
y = sin(u).*sin(v);
z = log(abs(tan(u/2))) + cos(u);
subplot(2,2,3)
mesh(x,y,z) %������������ͼ��α���棩
axis equal

t = 0:0.1:2.01*pi;
v = -1:0.1:1;
[v,t] = meshgrid(v,t); %���ɸ�����
r = 2+0.5*v.*cos(t/2);
x = r.*cos(t);
y = r.*sin(t);
z = 0.5*v.*sin(t/2);
subplot(2,2,4)
mesh(x,y,z) %������������ͼ�������˹���棩
axis equal

% close Figure 1 2

% ���ʺ����ʼ���
syms t real
r = [sin(t), cos(t), t];
r1 = diff(r);
r2 = diff(r1); % r2 = diff(r,2);
m = cross(r1,r2); % �������
% ��������
k = simplify((sqrt(m(1)^2 + m(2)^2 + m(3)^2))) / simplify(sqrt(r1(1)^2 + r1(2)^2 + r1(3)^2))^3;
sprintf('���ʣ�%s',k);

r3 = diff(r2);
C = [r1;r2;r3];
% ��������
tt = simplify(det(C) / sqrt(m(1)^2 + m(2)^2 + m(3)^2)^2);
sprintf('���ʣ�%s',tt);

% ���㻡��
s = int(norm(r1),t,0,2*pi);
sprintf('������%s',s);

% �����������
% �޶����ű������� real | positive | integer | rational
syms u v R real
r = {R*cos(u)*cos(v), R*cos(u)*sin(v), R*sin(u)}; %u = theta; v = phi
% r = {cos(u),sin(u),v};
ru = diff(r,u);
rv = diff(r,v);
simplify(cross(ru,rv));

% ����ĵ�һ������ʽ��dr = ru*du + rv*dv
% conj(Z) = real(Z) - i*imag(Z)
E = simplify(dot(ru,ru));
F = simplify(dot(ru,rv));
G = simplify(dot(rv,rv));
% class(G)
G2 = string(G);
% strfind�������ַ����е��ض��Ӵ����������Ӵ���һ���ַ���λ��
% strrep�����ַ������ض��Ӵ��滻Ϊ��һ�Ӵ�
% strtok���Կո�Ϊ��־��ȡ�ַ����е��Ӵ�
s_c = {"sin(u)^2 - 1", "1 - sin(u)^2", "cos(u)^2 - 1", "1 - cos(u)^2"};
s_c_ch = {"-cos(u)^2", "cos(u)^2", "-sin(u)^2", "sin(u)^2"};
class(s_c{1});
% �ж��ַ��������е��ַ����Ƿ�����ڱ����ҵ��ַ����У�����������Ӧ���滻
for i=1:length(s_c)
    if contains(G2, s_c{i})
        G = strrep(G2, s_c{i}, s_c_ch{i});
        G = simplify(str2sym(G));
    end
end

syms du2 du dv dv2
I = E*du2 + 2*F*du*dv + G*dv2;

% ������u������v���߼н�
theta = acos(F/sqrt(E*G));

% ������D����������ػ��ּ��㣩
area = int(int(sqrt(E*G-F),u,-pi/2,pi/2),v,0,2*pi);

% ����ڶ�������ʽ
ruu = diff(ru,u);
ruv = diff(ru,v);
rvv = diff(rv,v);

% ��ȥ����ֵ����
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

% ����������
k = II/I;

% ŷ����ʽ��kn = k1*cos(theta)^2 + k2*sin(theta)^2
% ��˹���ʣ�ȫ���ʡ������ʣ���K = k1*k2
K = (L*N-M^2)/(E*G-F^2);
% ƽ�����ʣ�H = (k1 + k2)/2
H = (L*G-2*M*F+N*E)/(2*(E*G-F^2));










