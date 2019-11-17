% 一个旅行商人要拜访全国31个省会城市，需要选择最短的路径
% 蚁群算法解决TSP问题

clear; close all; clc;

% 蚂蚁个数
m=50;
% 表征信息素重要程度的参数
Alpha=1;
% 表征启发式因子重要程度的参数
Beta=5;
% 信息素蒸发系数
Rho=0.1;
% 最大迭代次数
NC_max=100;
% 信息素增加强度系数
Q=100;

% 第一步：变量初始化
% n表示问题的规模（城市个数）
n=8;
% D表示完全图的赋权邻接矩阵
D=zeros(n,n);
% 只计算上对角元素的值
D(1,2)=3; D(1,3)=1; D(1,4)=2; D(1,5)=4; D(1,6)=3; D(1,7)=2; D(1,8)=1;
D(2,3)=3; D(2,4)=1; D(2,5)=2; D(2,6)=4; D(2,7)=1; D(2,8)=3;
D(3,4)=3; D(3,5)=3; D(3,6)=2; D(3,7)=4; D(3,8)=1;
D(4,5)=1; D(4,6)=2; D(4,7)=3; D(4,8)=2;
D(5,6)=1; D(5,7)=1; D(5,8)=3;
D(6,7)=2; D(6,8)=1;
D(7,8)=4;
% 对角元不能为0，但后面的启发因子要取倒数，用eps（浮点相对精度）表示
D = D + D' + diag(eps);

% Eta为启发因子，这里设为距离的倒数
Eta=1./D;
% Tau为信息素矩阵
Tau=ones(n,n);
% 存储并记录路径的生成
Tabu=zeros(m,n);
% 迭代计数器，记录迭代次数
NC=1; 
% 各代最佳路线
R_best=zeros(NC_max,n);
% 各代最佳路线的长度
L_best=inf.*ones(NC_max,1);
% 各代路线的平均长度
L_ave=zeros(NC_max,1);

% 停止条件之一：达到最大迭代次数，停止
while NC<=NC_max 
    % 第二步：将m只蚂蚁放到n个城市上
    % 随即存取
    Randpos=[];
    for i=1:(ceil(m/n))
        Randpos=[Randpos,randperm(n)];
    end
    Tabu(:,1)=(Randpos(1,1:m))';
    
    % 第三步：m只蚂蚁按概率函数选择下一座城市，完成各自的周游
    % 所在城市不计算
    for j=2:n 
        for i=1:m
            % 记录已访问的城市，避免重复访问
            visited=Tabu(i,1:(j-1)); 
            % 待访问的城市
            J=zeros(1,(n-j+1));
            % 待访问城市的选择概率分布
            P=J;
            Jc=1;
            for k=1:n
                % 开始时置0
                if isempty(find(visited==k, 1)) 
                    J(Jc)=k;
                    Jc=Jc+1;
                end
            end
            
            % 计算待选城市的概率分布
            for k=1:length(J)
                P(k)=(Tau(visited(end),J(k))^Alpha)*(Eta(visited(end),J(k))^Beta);
            end
            P=P/(sum(P));
            
            % 按概率原则选取下一个城市
            % cumsum元素累加即求和
            Pcum=cumsum(P); 
            % 若计算的概率大于原来的就选择这条路线
            Select=find(Pcum>=rand);
            to_visit=J(Select(1));
            Tabu(i,j)=to_visit;
        end
    end
    
    if NC>=2
        Tabu(1,:)=R_best(NC-1,:);
    end
    
    % 第四步：记录本次迭代最佳路线
    L=zeros(m,1); 
    for i=1:m
        R=Tabu(i,:);
        for j=1:(n-1)
            % 原距离加上第j个城市到第j+1个城市的距离
            L(i)=L(i)+D(R(j),R(j+1));
        end
        % 一轮下来后走过的距离
        L(i)=L(i)+D(R(1),R(n)); 
    end
    
    % 最佳距离取最小
    L_best(NC)=min(L); 
    pos=find(L==L_best(NC));
    % 此轮迭代后的最佳路线
    R_best(NC,:)=Tabu(pos(1),:);
    % 此轮迭代后的平均距离
    L_ave(NC)=mean(L);
    NC=NC+1;

    % 第五步：更新信息素
    % 开始时信息素为n*n的0矩阵
    Delta_Tau=zeros(n,n); 
    for i=1:m
        for j=1:(n-1)
            % 此次循环在路径（i,j）上的信息素增量
            Delta_Tau(Tabu(i,j),Tabu(i,j+1))=Delta_Tau(Tabu(i,j),Tabu(i,j+1))+Q/L(i);
        end
        % 此次循环在整个路径上的信息素增量
        Delta_Tau(Tabu(i,n),Tabu(i,1))=Delta_Tau(Tabu(i,n),Tabu(i,1))+Q/L(i);
    end
    
    % 考虑信息素挥发，更新后的信息素
    Tau=(1-Rho).*Tau+Delta_Tau; 
    
    % 第六步：禁忌表清零
    % 直到最大迭代次数
    Tabu=zeros(m,n);
end

%% 第七步：输出结果
% 找到最佳路径
Pos=find(L_best==min(L_best));
% 最大迭代次数后最佳路径
Shortest_Route=R_best(Pos(1),:)
% 最大迭代次数后最短距离
Shortest_Length=L_best(Pos(1))

figure(1)
plot(L_best)
xlabel('迭代次数')
ylabel('目标函数值')
title('适应度进化曲线')

plot(L_best)
hold on
plot(L_ave,'r')
title('平均距离和最短距离')






