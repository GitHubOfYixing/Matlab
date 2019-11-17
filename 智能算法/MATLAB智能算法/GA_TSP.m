% һ����������Ҫ�ݷ�8��ʡ����У���Ҫѡ����̵�·��
% �Ŵ��㷨���TSP����

clear; close all; clc;

% ��Ⱥ��С��m��Ⱦɫ�壩
m=100;
% �������Pc��ͨ��ȡ0.4-0.9֮��
Pc = 0.7;
% �������Pm��ͨ��ȡ0.001-0.1֮��
Pm = 0.1;
% ����������
N_max=200;
% ���д��ţ�����ʵ���������
city = [1 2 3 4 5 6 7 8];
city_num = length(city);
% ������֮��ľ���,��ȫͼ�ĸ�Ȩ�ڽӾ���
city_dst=zeros(city_num,city_num);
% ֻ�����϶Խ�Ԫ�ص�ֵ
city_dst(1,2)=3; city_dst(1,3)=1; city_dst(1,4)=2; city_dst(1,5)=4; city_dst(1,6)=3; city_dst(1,7)=2; city_dst(1,8)=1;
city_dst(2,3)=3; city_dst(2,4)=1; city_dst(2,5)=2; city_dst(2,6)=4; city_dst(2,7)=1; city_dst(2,8)=3;
city_dst(3,4)=3; city_dst(3,5)=3; city_dst(3,6)=2; city_dst(3,7)=4; city_dst(3,8)=1;
city_dst(4,5)=1; city_dst(4,6)=2; city_dst(4,7)=3; city_dst(4,8)=2;
city_dst(5,6)=1; city_dst(5,7)=1; city_dst(5,8)=3;
city_dst(6,7)=2; city_dst(6,8)=1;
city_dst(7,8)=4;
% �Խ�Ԫ����Ϊ0�����������������Ҫȡ��������eps��������Ծ��ȣ���ʾ
city_dst = city_dst + city_dst';

%% ��һ������ʼ����Ⱥ
S = zeros(m,city_num);
city_fitness = zeros(m,1);
for i = 1:m
    S(i,:) = randperm(city_num);
end

for ki=1:N_max
    for i = 1:m
        % ����ÿ���������Ӧ�Ⱥ���
        city_fit = 0;
        for j = 1:city_num-1
            city_fit = city_fit + city_dst(S(i,j),S(i,j+1));
        end
        city_fitness(i) = city_fit;
    end
    % ����ÿһ�ε��������е���Сֵ
    [f_min(ki),ri] = min(city_fitness);
    min_route(ki,:) = S(ri,:);

    %% �ڶ����������ۻ�����
    % ����ÿ�����屻ѡ�еĸ���
    Rw = city_fitness/sum(city_fitness);
    % ����ÿ��Ⱦɫ����ۻ�����
    Rp = zeros(1,m);
    for k=1:m
        Rp(k) = sum(Rw(1:k));
    end

    %% ��������ѡ��Ⱦɫ�� 
    % ģ������ת��������m���������ɸѡȾɫ��
    % �жϸ������λ���ĸ��ۻ��������䣬��ѡ���ĸ�Ⱦɫ�壬�ܹ�ѡ��m��Ⱦɫ��
    % ��ЩȾɫ�����п��ܱ�ѡ��һ�λ��Σ�Ҳ�п���һ�ζ�û�б�ѡ��
    K = zeros(1,m);
    for i=1:m
        Pr = rand();
        K(i) = find(Pr<Rp,1);
    end
    % ���ձ�ѡ�ϵ�Ⱦɫ��
    Kf = S(K,:);

    %% ���Ĳ������е��㽻��
    % ��Ⱦɫ��������
    % Ⱥ�����������ԣ���m/2�ԣ���Է�ʽ��(1,m/2+1),(2,m/2+2),(3,m/2+3),...��
    % ��ǰһ���Ⱦɫ����Ϊ��Ⱦɫ�壬��һ���Ⱦɫ����ΪĸȾɫ�壬����˳��������
    Nr = rand(1,m/2);
    % ���ݽ������ȷ��ʵ�ʽ��������m/2 =��A+B����A��ʵ�ʷ��������Ⱦɫ��ԣ�B�������������Ⱦɫ��ԣ�
    % �������m/2����������������뽻�����Pc�Ƚ�
    % ��С�ڽ�����ʣ����ʾ���н���
    flag = find(Nr<Pc);
    % ȷ����Ҫ���������Ⱦɫ��ԣ���ǰһ���Ⱦɫ����Ϊ��Ⱦɫ�壬��һ���Ⱦɫ����ΪĸȾɫ�壬����˳�������ԣ�
    Kfy = Kf(flag,:);
    Kfx = Kf(flag+m/2,:);
    % sprintf('�ɹ����Ⱦɫ�������%d',2*size(Kfy,1))

    % ����δ���������Ⱦɫ��
    flag_ = find(Nr>=Pc);
    % δ��������ĸ�Ⱦɫ��
    Kfy_ = Kf(flag_,:);
    % δ���������ĸȾɫ��
    Kfx_ = Kf(flag_+m/2,:);
    % δ�����������Ⱦɫ��
    Kf_ = [Kfy_;Kfx_];

    % ����Ժõ�Ⱦɫ����е��㽻��
    % �����city_num-1��������Ͻ��е��㽻��
    % �Բ�ͬ�Ե�Ⱦɫ��������ѡ�񽻲��λ�ý��н���
    Ac = '';
    Bc = '';
    for i=1:size(Kfy,1)
        A = Kfy(i,:);
        B = Kfx(i,:);
        % ���ȷ������λ��
        index = randperm(city_num-1,1);
        % ����λ�ý��淨
        % ���Ϊ����������Ⱦɫ�壬�ұ�Ϊ�Ӵ�Ⱦɫ�塣 
        % (1)�����ϵ������һ��Ԫ�ط�����������ĵ�һλ����
        % (2)��ת�Ƶ����������һ��Ԫ��
        % (3)�鿴���������Ƿ��Ѿ������˸�Ԫ��
        % (4)���δ��������������������У�����������������С�
        % (5)���Ŵ���������ĵڶ���Ԫ�ؿ�ʼ�������µڶ���Ԫ�أ���ǰ��ͬ����жϲ�����
        % (6)�������ֱ���ұ��������鶼��������Ϊֹ��
        Ac = A(1:index);
        Bc = B(1:index);
        flagA = index+1;
        flagB = index+1;
        while index <= city_num
            if flagA <= city_num 
                % ���A������Ԫ��û����Ac�г��֣���׷����Ac��
                if isempty(find(A(index)==Ac, 1))
                    Ac(flagA) = A(index); 
                    flagA = flagA + 1;
                end
                % ���B������Ԫ��û����Ac�г��֣���׷����Ac��
                if isempty(find(B(index)==Ac, 1))
                    Ac(flagA) = B(index);
                    flagA = flagA + 1;
                end
            end
            if flagB <= city_num 
                % ���B������Ԫ��û����Bc�г��֣���׷����Bc��
                if isempty(find(B(index)==Bc, 1))
                    Bc(flagB) = B(index); 
                    flagB = flagB + 1;
                end
                % ���A������Ԫ��û����Bc�г��֣���׷����Bc��
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

    %% ���岽������
    % λ�õ���������Ⱦɫ��������Ĳ�������λ������ֵ����
    for i=1:size(C,1)
        % ����city_num�������
        n = rand(1,city_num);
        % ��city_num��������������ʽ��бȽϣ��ж��Ƿ���ڱ������λ���Ƿ���Ҫ���б��죩
        k = find(n < Pm);
        % �����ڣ������Ӧ����λ���б��촦��
        if ~isempty(k)
            % ���б��죨�п���ֻ��һ��λ�÷���������Ҳ�п��ܶ��λ�÷���������
            % ����k������λ�����city_num-k+1����λ�ý��л���
            C(i,[k city_num-k+1]) = C(i,[city_num-k+1 k]);
        end
    end

    %% ������Ⱥ������Ⱥ = ��ѡ�е�δ��������ĸ�Ⱥ + ��ѡ���ҷ������䣨���죩�����Ⱥ��
    S = [Kf_;C];
end
[min_distance,k] = min(f_min)
route_min = min_route(k,:) 

plot(f_min)
