% 一个旅行商人要拜访8个省会城市，需要选择最短的路径
% 遗传算法解决TSP问题

clear; close all; clc;

% 种群大小（m个染色体）
m=100;
% 交配概率Pc：通常取0.4-0.9之间
Pc = 0.7;
% 变异概率Pm：通常取0.001-0.1之间
Pm = 0.1;
% 最大迭代次数
N_max=200;
% 城市代号（采用实数编码规则）
city = [1 2 3 4 5 6 7 8];
city_num = length(city);
% 各城市之间的距离,完全图的赋权邻接矩阵
city_dst=zeros(city_num,city_num);
% 只计算上对角元素的值
city_dst(1,2)=3; city_dst(1,3)=1; city_dst(1,4)=2; city_dst(1,5)=4; city_dst(1,6)=3; city_dst(1,7)=2; city_dst(1,8)=1;
city_dst(2,3)=3; city_dst(2,4)=1; city_dst(2,5)=2; city_dst(2,6)=4; city_dst(2,7)=1; city_dst(2,8)=3;
city_dst(3,4)=3; city_dst(3,5)=3; city_dst(3,6)=2; city_dst(3,7)=4; city_dst(3,8)=1;
city_dst(4,5)=1; city_dst(4,6)=2; city_dst(4,7)=3; city_dst(4,8)=2;
city_dst(5,6)=1; city_dst(5,7)=1; city_dst(5,8)=3;
city_dst(6,7)=2; city_dst(6,8)=1;
city_dst(7,8)=4;
% 对角元不能为0，但后面的启发因子要取倒数，用eps（浮点相对精度）表示
city_dst = city_dst + city_dst';

%% 第一步：初始化种群
S = zeros(m,city_num);
city_fitness = zeros(m,1);
for i = 1:m
    S(i,:) = randperm(city_num);
end

for ki=1:N_max
    for i = 1:m
        % 计算每个个体的适应度函数
        city_fit = 0;
        for j = 1:city_num-1
            city_fit = city_fit + city_dst(S(i,j),S(i,j+1));
        end
        city_fitness(i) = city_fit;
    end
    % 计算每一次迭代过程中的最小值
    [f_min(ki),ri] = min(city_fitness);
    min_route(ki,:) = S(ri,:);

    %% 第二步：计算累积概率
    % 计算每个个体被选中的概率
    Rw = city_fitness/sum(city_fitness);
    % 计算每个染色体的累积概率
    Rp = zeros(1,m);
    for k=1:m
        Rp(k) = sum(Rw(1:k));
    end

    %% 第三步：选择染色体 
    % 模拟轮盘转动，产生m个随机数，筛选染色体
    % 判断该随机数位于哪个累积概率区间，就选择哪个染色体，总共选择m个染色体
    % 这些染色体中有可能被选中一次或多次，也有可能一次都没有被选中
    K = zeros(1,m);
    for i=1:m
        Pr = rand();
        K(i) = find(Pr<Rp,1);
    end
    % 最终被选上的染色体
    Kf = S(K,:);

    %% 第四步：进行单点交叉
    % 将染色体进行配对
    % 群体个体两两配对，共m/2对（配对方式：(1,m/2+1),(2,m/2+2),(3,m/2+3),...）
    % 将前一半的染色体作为父染色体，后一半的染色体作为母染色体，并按顺序进行配对
    Nr = rand(1,m/2);
    % 根据交配概率确定实际交配对数：m/2 =（A+B）（A：实际发生交配的染色体对；B：不发生交配的染色体对）
    % 随机产生m/2个随机数，并将其与交配概率Pc比较
    % 若小于交配概率，则表示进行交叉
    flag = find(Nr<Pc);
    % 确定需要发生交叉的染色体对（将前一半的染色体作为父染色体，后一半的染色体作为母染色体，并按顺序进行配对）
    Kfy = Kf(flag,:);
    Kfx = Kf(flag+m/2,:);
    % sprintf('成功配对染色体个数：%d',2*size(Kfy,1))

    % 保存未发生交叉的染色体
    flag_ = find(Nr>=Pc);
    % 未发生交配的父染色体
    Kfy_ = Kf(flag_,:);
    % 未发生交配的母染色体
    Kfx_ = Kf(flag_+m/2,:);
    % 未发生交配的总染色体
    Kf_ = [Kfy_;Kfx_];

    % 对配对好的染色体进行单点交叉
    % 随机在city_num-1个交叉点上进行单点交叉
    % 对不同对的染色体进行随机选择交叉点位置进行交叉
    Ac = '';
    Bc = '';
    for i=1:size(Kfy,1)
        A = Kfy(i,:);
        B = Kfx(i,:);
        % 随机确定交叉位置
        index = randperm(city_num-1,1);
        % 交替位置交叉法
        % 左边为父代的两个染色体，右边为子代染色体。 
        % (1)将左上的数组第一个元素放入右上数组的第一位置中
        % (2)再转移到左下数组第一个元素
        % (3)查看右上数组是否已经包含了该元素
        % (4)如果未包含将其插入右上数组中，否则插入右下数组中。
        % (5)接着从左上数组的第二个元素开始，到左下第二个元素，和前次同意的判断操作。
        % (6)如此类推直到右边两个数组都被填满了为止。
        Ac = A(1:index);
        Bc = B(1:index);
        flagA = index+1;
        flagB = index+1;
        while index <= city_num
            if flagA <= city_num 
                % 如果A数组中元素没有在Ac中出现，则追加在Ac中
                if isempty(find(A(index)==Ac, 1))
                    Ac(flagA) = A(index); 
                    flagA = flagA + 1;
                end
                % 如果B数组中元素没有在Ac中出现，则追加在Ac中
                if isempty(find(B(index)==Ac, 1))
                    Ac(flagA) = B(index);
                    flagA = flagA + 1;
                end
            end
            if flagB <= city_num 
                % 如果B数组中元素没有在Bc中出现，则追加在Bc中
                if isempty(find(B(index)==Bc, 1))
                    Bc(flagB) = B(index); 
                    flagB = flagB + 1;
                end
                % 如果A数组中元素没有在Bc中出现，则追加在Bc中
                if isempty(find(A(index)==Bc, 1))
                    Bc(flagB) = A(index);
                    flagB = flagB + 1;
                end
            end
            index = index+1;
        end
        Mc(i,:) = Ac;
        Nc(i,:) = Bc;
    end
    C = [Mc;Nc];

    %% 第五步：变异
    % 位置倒换法，即染色体上随机的产生两个位置上数值互换
    for i=1:size(C,1)
        % 产生city_num个随机数
        n = rand(1,city_num);
        % 将city_num个随机数与变异概率进行比较，判断是否存在变异基因位（是否需要进行变异）
        k = find(n < Pm);
        % 若存在，则对相应基因位进行变异处理
        if ~isempty(k)
            % 进行变异（有可能只有一个位置发生交换，也有可能多个位置发生交换）
            % 将第k个变异位置与第city_num-k+1变异位置进行互换
            C(i,[k city_num-k+1]) = C(i,[city_num-k+1 k]);
        end
    end

    %% 更新种群（新种群 = 被选中但未发生交配的父群 + 被选中且发生交配（变异）后的子群）
    S = [Kf_;C];
end
[min_distance,k] = min(f_min)
route_min = min_route(k,:) 

plot(f_min)
