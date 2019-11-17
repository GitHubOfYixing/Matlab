% һ����������Ҫ�ݷ�ȫ��31��ʡ����У���Ҫѡ����̵�·��
% ��Ⱥ�㷨���TSP����

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
NC_max=100;
% ��Ϣ������ǿ��ϵ��
Q=100;

% ��һ����������ʼ��
% n��ʾ����Ĺ�ģ�����и�����
n=8;
% D��ʾ��ȫͼ�ĸ�Ȩ�ڽӾ���
D=zeros(n,n);
% ֻ�����϶Խ�Ԫ�ص�ֵ
D(1,2)=3; D(1,3)=1; D(1,4)=2; D(1,5)=4; D(1,6)=3; D(1,7)=2; D(1,8)=1;
D(2,3)=3; D(2,4)=1; D(2,5)=2; D(2,6)=4; D(2,7)=1; D(2,8)=3;
D(3,4)=3; D(3,5)=3; D(3,6)=2; D(3,7)=4; D(3,8)=1;
D(4,5)=1; D(4,6)=2; D(4,7)=3; D(4,8)=2;
D(5,6)=1; D(5,7)=1; D(5,8)=3;
D(6,7)=2; D(6,8)=1;
D(7,8)=4;
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

plot(L_best)
hold on
plot(L_ave,'r')
title('ƽ���������̾���')






