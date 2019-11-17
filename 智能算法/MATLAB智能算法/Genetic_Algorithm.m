% 遗传算法求解函数的最小值：f(x) = x - 4*x^2 + x^4;
% clc;clear;close all;
function [g_min,Sdecu,J_min] = Genetic_Algorithm(f,fitness,umin,umax,delta,Pc,Pm,N_r)

global S;

syms x
% f = x - 8*x^2 + 3*x^4;
% f = "x - 8*x.^2 + 3*x.^4";
F = eval(sprintf('@(x) %s',f));
Fitness = eval(sprintf('@(x) %s',fitness));

% fx = vpa(vpasolve(diff(F,x)==0,x),5);
% g = matlabFunction(f);
% 所有极值点
% g_all = roundn(double(F(fx)),-4);
% 实际最小值
% g_min =  roundn(double(min(g_all)),-4);
g_min = min(F(umin:delta:umax));

% umin = -3;
% umax = 3;
% delta = 0.0001;
% (1) 确定变量的取值范围:(-3,3); 编码精度: 0.0001
% 进行二进制编码
% 计算需要的二进制位数(delta = (umax - umin)/(2^lamda - 1); lamda = log((umax - umin)/delta + 1))
lamda = ceil(log2((umax - umin)/delta + 1)); % 向上取整：lamda=16
% (2)确定适应度函数fitness(x): 单值、连续、非负、最大化
% 目标函数：f(x) = x - 4*x^2 + x^4：开口向上，存在小于0的部分，是最小值问题
% 适应度函数：f = 20 - (x - 4*x^2 + x^4)：转换为极大值问题

% (3)设置群体规模N：通常介于(n,2n)=(16,32)
N = ceil(lamda*3/2);
% 将种群设置为偶数，方便进行后续的配对
if mod(N,2)==1
    N = N+1;
end
% sprintf("种群大小：%d",N)

% 产生初始种群（随机产生）（编码）
S = dec2bin(randi(2^lamda,N,1)-1,lamda);
% (4)设置交配概率Pc：通常取0.4-0.9之间
% Pc = 0.7;
% (5)设置变异概率Pm：通常取0.001-0.1之间
% Pm = 0.01;

%% 开始进行迭代
% 迭代次数
% N_r = 300;
f_max = zeros(1,N_r);
delta = (umax - umin)/(2^lamda - 1);
for ki=1:N_r
    
    % (6)进行复制：比例选择(轮盘赌选择法)
    % 根据适应度函数计算个体的适应度
    % 将二进制数转换成十进制数（解码）
    Sdec = bin2dec(S);
    % 将十进制数转换成[umin,umax]之间的数（解码）
    Sdecu = umin + Sdec*delta;
    % 计算每个个体的适应度值
    % fitness = 20 - F(Sdecu);
    Fit = Fitness(Sdecu);
    % 计算每一次迭代过程中的最大值
    f_max(ki) = max(Fit);
    % 计算每个染色体被选中的概率
    Rw = Fit/sum(Fit);

    % 计算每个染色体的累积概率
    Rp = zeros(1,N);
    for k=1:N
        Rp(k) = sum(Rw(1:k));
    end

    % 模拟轮盘转动，产生N个随机数，筛选染色体
    % 判断该随机数位于哪个累积概率区间，就选择哪个染色体，总共选择N个染色体
    % 这些染色体中有可能被选中一次或多次，也有可能一次都没有被选中
    K = zeros(1,N);
    for i=1:N
        Pr = rand();
        K(i) = find(Pr<Rp,1);
    end
    % 最终被选上的染色体
    Kf = S(K,:);
    
    % (7)进行交叉：单点交叉
    % 群体个体两两配对，共N/2对（配对方式：(1,N/2+1),(2,N/2+2),(3,N/2+3),...）
    % 将前一半的染色体作为父染色体，后一半的染色体作为母染色体，并按顺序进行配对
    
    % 根据交配概率确定实际交配对数：N/2 =（A+B）（A：实际发生交配的染色体对；B：不发生交配的染色体对）
    % 随机产生N/2个随机数，并将其与交配概率Pc比较
    Nr = rand(1,N/2);
    % 若小于交配概率，则表示进行交叉
    flag = find(Nr<Pc);
    % 确定需要发生交叉的染色体对（将前一半的染色体作为父染色体，后一半的染色体作为母染色体，并按顺序进行配对）
    Kfy = Kf(flag,:);
    Kfx = Kf(flag+N/2,:);
    % sprintf("成功配对染色体数量：%d",2*size(Kfy,1))
    
    % 保存未发生交叉的染色体
    flag_ = find(Nr>=Pc);
    % 未发生交配的父染色体
    Kfy_ = Kf(flag_,:);
    % 未发生交配的母染色体
    Kfx_ = Kf(flag_+N/2,:);
    % 未发生交配的总染色体
    Kf_ = [Kfy_;Kfx_];

    % 对配对好的染色体进行单点交叉
    % 随机在lamda-1个交叉点上进行单点交叉
    % 对不同对的染色体进行随机选择交叉点位置进行交叉
    Ac = '';
    Bc = '';
    for i=1:size(Kfy,1)
        A = Kfy(i,:);
        B = Kfx(i,:);
        % 随机确定交叉位置
        index = randperm(lamda-1,1);
        % [A;B]
        % 交叉操作
        temp = A(index+1:end);
        Ac(i,:) = [A(1:index),B(index+1:end)];
        Bc(i,:) = [B(1:index),temp];
    end
    C = [Ac;Bc];

    % (8)变异：按基因位变异（变异基因数量：B = Pm*M*L = 0.01*size(C,1)*lamda）
    Kc1 = 0;
    Kc2 = 0;
    for i=1:size(C,1)
        % 产生lamda个随机数
        n = rand(1,lamda);
        % 将lamda个随机数与变异概率进行比较，判断是否存在变异基因位（是否需要进行变异）
        ps = find(n < Pm);
        % 若存在，则对相应基因位进行变异处理
        if ~isempty(ps)
            % 进行基因位变异（有可能只有一个基因位发生改变，也有可能多个基因位发生改变）
            Kc1 = Kc1 + 1; % 统计变异染色体数量
            for k=1:length(ps)
                C(i,:) = [C(i,1:ps(k)-1),dec2bin(~bin2dec(C(i,ps(k)))),C(i,ps(k)+1:end)];
                Kc2 = Kc2 + 1; % 统计变异基因数量
            end
        end
    end
%     sprintf("变异染色体数量：%d",Kc1)
%     sprintf("变异基因数量：%d",Kc2)
%     sprintf("变异率：%.3f",Kc2/(size(C,1)*lamda))

    % (9)更新种群（新种群 = 被选中但未发生交配的父群 + 被选中且发生交配（变异）后的子群）
    S = [Kf_;C];
end

% 遗传算法求解的最小值
f_min = 20 - f_max;
% 统计数组中各元素出现的次数与频数
J = tabulate(f_min);
% 计算最小值出现的位置（可能有多个最小值）
J_cnt = J(:,1)==min(J(:,1));
% 计算最小值
J_min = roundn(double(min(J(:,1))),-4);
% 计算得到最小值时的变量X的值
Sdecu = Sdecu(J_cnt);


