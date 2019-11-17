% һ����������Ҫ�ݷ�ȫ��31��ʡ����У���Ҫѡ����̵�·��
% ��Ⱥ�㷨���TSP����
% C      n�����е����꣬n��2�ľ���
% NC_max ����������
% m      ���ϸ���
% Alpha  ������Ϣ����Ҫ�̶ȵĲ���
% Beta   ��������ʽ������Ҫ�̶ȵĲ���
% Rho    ��Ϣ������ϵ��
% Q      ��Ϣ������ǿ��ϵ��
% R_best �������·��
% L_best �������·�ߵĳ���

clear; close all; clc;

% ���ϸ���
m=50;
% ������Ϣ����Ҫ�̶ȵĲ���
Alpha=1;
% ��������ʽ������Ҫ�̶ȵĲ���
Beta=5;
% ��Ϣ������ϵ��
Rho=0.1;
% ����������
NC_max=200;
% ��Ϣ������ǿ��ϵ��
Q=100;

% 31��ʡ������
C=[
1304 2312;
3639 1315;
4177 2244;
3712 1399;
3488 1535;
3326 1556;
3238 1229;
4196 1004;
4312 790;
4386 570;
3007 1970;
2562 1756;
2788 1491;
2381 1676;
1332 695;
3715 1678;
3918 2179; 
4061 2370;
3780 2212;
3676 2578;
4029 2838;
4263 2931;
3429 1908;
3507 2367;
3394 2643;
3439 3201;
2935 3240;
3140 3550;
2545 2357;
2778 2826;
2370 2975];

% ��һ����������ʼ��
% n��ʾ����Ĺ�ģ�����и�����
n=size(C,1);
% D��ʾ��ȫͼ�ĸ�Ȩ�ڽӾ���
D=zeros(n,n);
% ֻ�����϶Խ�Ԫ�ص�ֵ
for i=1:n
    for j=i:n
        % ����������������֮��ľ���
        D(i,j)=sqrt((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2);
    end
end
% �Խ�Ԫ����Ϊ0�����������������Ҫȡ��������eps��������Ծ��ȣ���ʾ
D = D + D' + diag(eps);

% EtaΪ�������ӣ�������Ϊ����ĵ���
Eta=1./D;
% TauΪ��Ϣ�ؾ���
Tau=ones(n,n);
% �洢����¼·��������
Tabu=zeros(m,n);
% ��������������¼��������
NC=1; 
% �������·��
R_best=zeros(NC_max,n);
% �������·�ߵĳ���
L_best=inf.*ones(NC_max,1);
% ����·�ߵ�ƽ������
L_ave=zeros(NC_max,1);

% ֹͣ����֮һ���ﵽ������������ֹͣ
while NC<=NC_max 
    % �ڶ�������mֻ���Ϸŵ�n��������
    % �漴��ȡ
    Randpos=[];
    for i=1:(ceil(m/n))
        Randpos=[Randpos,randperm(n)];
    end
    Tabu(:,1)=(Randpos(1,1:m))';
    
    % ��������mֻ���ϰ����ʺ���ѡ����һ�����У���ɸ��Ե�����
    % ���ڳ��в�����
    for j=2:n 
        for i=1:m
            % ��¼�ѷ��ʵĳ��У������ظ�����
            visited=Tabu(i,1:(j-1)); 
            % �����ʵĳ���
            J=zeros(1,(n-j+1));
            % �����ʳ��е�ѡ����ʷֲ�
            P=J;
            Jc=1;
            for k=1:n
                % ��ʼʱ��0
                if isempty(find(visited==k, 1)) 
                    J(Jc)=k;
                    Jc=Jc+1;
                end
            end
            
            % �����ѡ���еĸ��ʷֲ�
            for k=1:length(J)
                P(k)=(Tau(visited(end),J(k))^Alpha)*(Eta(visited(end),J(k))^Beta);
            end
            P=P/(sum(P));
            
            % ������ԭ��ѡȡ��һ������
            % cumsumԪ���ۼӼ����
            Pcum=cumsum(P); 
            % ������ĸ��ʴ���ԭ���ľ�ѡ������·��
            Select=find(Pcum>=rand);
            to_visit=J(Select(1));
            Tabu(i,j)=to_visit;
        end
    end
    
    if NC>=2
        Tabu(1,:)=R_best(NC-1,:);
    end
    
    % ���Ĳ�����¼���ε������·��
    L=zeros(m,1); 
    for i=1:m
        R=Tabu(i,:);
        for j=1:(n-1)
            % ԭ������ϵ�j�����е���j+1�����еľ���
            L(i)=L(i)+D(R(j),R(j+1));
        end
        % һ���������߹��ľ���
        L(i)=L(i)+D(R(1),R(n)); 
    end
    
    % ��Ѿ���ȡ��С 
    L_best(NC)=min(L); 
    pos=find(L==L_best(NC));
    % ���ֵ���������·��
    R_best(NC,:)=Tabu(pos(1),:);
    % ���ֵ������ƽ������
    L_ave(NC)=mean(L);
    NC=NC+1;

    % ���岽��������Ϣ��
    % ��ʼʱ��Ϣ��Ϊn*n��0����
    Delta_Tau=zeros(n,n); 
    for i=1:m
        for j=1:(n-1)
            % �˴�ѭ����·����i,j���ϵ���Ϣ������
            Delta_Tau(Tabu(i,j),Tabu(i,j+1))=Delta_Tau(Tabu(i,j),Tabu(i,j+1))+Q/L(i);
        end
        % �˴�ѭ��������·���ϵ���Ϣ������
        Delta_Tau(Tabu(i,n),Tabu(i,1))=Delta_Tau(Tabu(i,n),Tabu(i,1))+Q/L(i);
    end
    
    % ������Ϣ�ػӷ������º����Ϣ��
    Tau=(1-Rho).*Tau+Delta_Tau; 
    
    % �����������ɱ�����
    % ֱ������������
    Tabu=zeros(m,n);
end

%% ���߲���������
% �ҵ����·��
Pos=find(L_best==min(L_best));
% ���������������·��
Shortest_Route=R_best(Pos(1),:)
% ��������������̾���
Shortest_Length=L_best(Pos(1))

figure(1)
plot(L_best)
xlabel('��������')
ylabel('Ŀ�꺯��ֵ')
title('��Ӧ�Ƚ�������')

figure(2)
% ���Ƶ�һ����ͼ��
subplot(1,2,1)
% ��·��ͼ
% C Coordinate�ڵ����꣬��һ��N��2�ľ���洢
% R Route ·��
N=length(R);
scatter(C(:,1),C(:,2));
hold on;
plot([C(R(1),1),C(R(N),1)],[C(R(1),2),C(R(N),2)],'g')
hold on;
for ii=2:N
    plot([C(R(ii-1),1),C(R(ii),1)],[C(R(ii-1),2),C(R(ii),2)],'g')
    hold on;
end
title('�����������Ż����')

% ���Ƶڶ�����ͼ��
subplot(1,2,2)
plot(L_best)
hold on
plot(L_ave,'r')
title('ƽ���������̾���')






