% �����Լ���Ż����⣨Rosenbrock������ά��Ϊ5��
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
mesh(x,y,z);%��ά����ͼ
% meshc(x,y,z);%���ȸ��ߵ���ά����ͼ
% meshz(x,y,z);%����������ά����ͼ
% surf(x,y,z);%��ά����ͼ
title('z = 100*x*y*exp(-x^2-y^2-x*y/5)');
% title('z = (1-x).^2+100*(y-x.^2).^2');
xlabel('x');
ylabel('y');
zlabel('z');

% ��������
Index = 1;
% ����������
MaxIndex = 100;
% ��֪ϵ��
c1 = 0.6;
% ���ѧϰϵ��
c2 = 0.7;
% ����ϵ����������������Ӷ��ݼ�
wmax = 0.9;
wmin = 0.4;
w = wmax - (wmax-wmin)*(Index/MaxIndex)^2;
% �����������ٶȣ�-3,3��
vmax = 6;
% ��Ӧ�Ⱥ���
fitness = @(x)100*x(:,1).*x(:,2).*exp(-x(:,1).^2-x(:,2).^2-x(:,1).*x(:,2)/5);
% fitness = @(x)(1-x(:,1)).^2+100*(x(:,2)-x(:,1).^2).^2;
% fitness = @(x)(1-x(:,1)).^2+100*(x(:,2)-x(:,1).^2).^2 + (1-x(:,2)).^2+100*(x(:,3)-x(:,2).^2).^2;
% �����ʼ��������Ⱥ�ٶ�V��λ��X
% ������Ⱥ��СΪ20��ÿ�����ӵ�ά��Ϊ2ά
N = 20; Dim = 2;
% ��ʼ��ÿ������λ��x
x = -3 + 6*rand(N,Dim);
% ��ʼ��ÿ�������ٶ�v
v = vmax*rand(N,Dim);
% ����ÿ���������Ӧ��
f = fitness(x);

% ���ó�ʼ������ʷ�����Ӧ��
fpbest = f;
% ���ó�ʼ������ʷ���λ��
pbest = x;
% ���ó�ʼ��Ⱥ��ʷ�����Ӧ�ȣ�iΪ�����Ӧ������Ӧ�����λ��
[fgbest,i] = min(fpbest);
% ���ó�ʼ��Ⱥ��ʷ���λ��
gbest = pbest(i,:);
                    
% ��ʼ����
figure(2)
t = zeros(1,MaxIndex);
y = zeros(1,MaxIndex);
while Index <= MaxIndex
    scatter(x(:,1),x(:,2),'yo'); % ��ά��̬��ʾ���ӵĸ��¹���
    title(sprintf('fmin(%.4f)',fgbest))
    axis([-3 3 -3 3]);

%     scatter3(x(:,1),x(:,2),f,'yo'); % ��ά��̬��ʾ���ӵĸ��¹���
%     axis([-3 3 -3 3 -20 20]);

    % ������������
%     t(Index) = Index;
%     y(Index) = fgbest;
%     plot(t,y,'b','MarkerSize',5);
%     p = plot(t,y,'b','MarkerSize',5);
%     set(p,'XData',t,'YData',y)
%     drawnow
%     axis([0 MaxIndex -30 -10]);
    
    % ���¹���Ȩֵ���ϴ�Ȩֵ������ȫ����������СȨֵ�����ھֲ�������
    w = wmax - (wmax-wmin)*(Index/MaxIndex)^2;
    % ��������
    % �ɸ������Ⱥ�������Ϣ���ݽ��и����ƶ��ٶȣ�ע�⣺������Ҫ��ȫ������ֵ������ά����ʹ����xά��һ�£�
    v = w*v + c1*rand()*(pbest - x) + c2*rand()*(ones(N,1)*gbest - x); 
    % ͨ���ٶȸ���λ��
    x = x + v;
    
    % ����ÿ��������Ӧ��
    f = fitness(x);
    % �����Ƿ���ڸ�С����Ӧ��ֵ
    k = find(f(:) < fpbest(:));
    % ������ڣ���������Ӹ���
    if size(k)
        % ���¸��������Ӧ��
        fpbest(k) = f(k);
        % �������Ÿ���λ��
        pbest(k,:) = x(k,:);
    end
    % ������Ⱥ�����Ӧ��
    [fgbest,i] = min(fpbest); 
    % ������Ⱥ���λ��
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


