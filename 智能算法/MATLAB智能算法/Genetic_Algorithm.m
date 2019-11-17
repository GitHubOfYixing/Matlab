% �Ŵ��㷨��⺯������Сֵ��f(x) = x - 4*x^2 + x^4;
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
% ���м�ֵ��
% g_all = roundn(double(F(fx)),-4);
% ʵ����Сֵ
% g_min =  roundn(double(min(g_all)),-4);
g_min = min(F(umin:delta:umax));

% umin = -3;
% umax = 3;
% delta = 0.0001;
% (1) ȷ��������ȡֵ��Χ:(-3,3); ���뾫��: 0.0001
% ���ж����Ʊ���
% ������Ҫ�Ķ�����λ��(delta = (umax - umin)/(2^lamda - 1); lamda = log((umax - umin)/delta + 1))
lamda = ceil(log2((umax - umin)/delta + 1)); % ����ȡ����lamda=16
% (2)ȷ����Ӧ�Ⱥ���fitness(x): ��ֵ���������Ǹ������
% Ŀ�꺯����f(x) = x - 4*x^2 + x^4���������ϣ�����С��0�Ĳ��֣�����Сֵ����
% ��Ӧ�Ⱥ�����f = 20 - (x - 4*x^2 + x^4)��ת��Ϊ����ֵ����

% (3)����Ⱥ���ģN��ͨ������(n,2n)=(16,32)
N = ceil(lamda*3/2);
% ����Ⱥ����Ϊż����������к��������
if mod(N,2)==1
    N = N+1;
end
% sprintf("��Ⱥ��С��%d",N)

% ������ʼ��Ⱥ����������������룩
S = dec2bin(randi(2^lamda,N,1)-1,lamda);
% (4)���ý������Pc��ͨ��ȡ0.4-0.9֮��
% Pc = 0.7;
% (5)���ñ������Pm��ͨ��ȡ0.001-0.1֮��
% Pm = 0.01;

%% ��ʼ���е���
% ��������
% N_r = 300;
f_max = zeros(1,N_r);
delta = (umax - umin)/(2^lamda - 1);
for ki=1:N_r
    
    % (6)���и��ƣ�����ѡ��(���̶�ѡ��)
    % ������Ӧ�Ⱥ�������������Ӧ��
    % ����������ת����ʮ�����������룩
    Sdec = bin2dec(S);
    % ��ʮ������ת����[umin,umax]֮����������룩
    Sdecu = umin + Sdec*delta;
    % ����ÿ���������Ӧ��ֵ
    % fitness = 20 - F(Sdecu);
    Fit = Fitness(Sdecu);
    % ����ÿһ�ε��������е����ֵ
    f_max(ki) = max(Fit);
    % ����ÿ��Ⱦɫ�屻ѡ�еĸ���
    Rw = Fit/sum(Fit);

    % ����ÿ��Ⱦɫ����ۻ�����
    Rp = zeros(1,N);
    for k=1:N
        Rp(k) = sum(Rw(1:k));
    end

    % ģ������ת��������N���������ɸѡȾɫ��
    % �жϸ������λ���ĸ��ۻ��������䣬��ѡ���ĸ�Ⱦɫ�壬�ܹ�ѡ��N��Ⱦɫ��
    % ��ЩȾɫ�����п��ܱ�ѡ��һ�λ��Σ�Ҳ�п���һ�ζ�û�б�ѡ��
    K = zeros(1,N);
    for i=1:N
        Pr = rand();
        K(i) = find(Pr<Rp,1);
    end
    % ���ձ�ѡ�ϵ�Ⱦɫ��
    Kf = S(K,:);
    
    % (7)���н��棺���㽻��
    % Ⱥ�����������ԣ���N/2�ԣ���Է�ʽ��(1,N/2+1),(2,N/2+2),(3,N/2+3),...��
    % ��ǰһ���Ⱦɫ����Ϊ��Ⱦɫ�壬��һ���Ⱦɫ����ΪĸȾɫ�壬����˳��������
    
    % ���ݽ������ȷ��ʵ�ʽ��������N/2 =��A+B����A��ʵ�ʷ��������Ⱦɫ��ԣ�B�������������Ⱦɫ��ԣ�
    % �������N/2����������������뽻�����Pc�Ƚ�
    Nr = rand(1,N/2);
    % ��С�ڽ�����ʣ����ʾ���н���
    flag = find(Nr<Pc);
    % ȷ����Ҫ���������Ⱦɫ��ԣ���ǰһ���Ⱦɫ����Ϊ��Ⱦɫ�壬��һ���Ⱦɫ����ΪĸȾɫ�壬����˳�������ԣ�
    Kfy = Kf(flag,:);
    Kfx = Kf(flag+N/2,:);
    % sprintf("�ɹ����Ⱦɫ��������%d",2*size(Kfy,1))
    
    % ����δ���������Ⱦɫ��
    flag_ = find(Nr>=Pc);
    % δ��������ĸ�Ⱦɫ��
    Kfy_ = Kf(flag_,:);
    % δ���������ĸȾɫ��
    Kfx_ = Kf(flag_+N/2,:);
    % δ�����������Ⱦɫ��
    Kf_ = [Kfy_;Kfx_];

    % ����Ժõ�Ⱦɫ����е��㽻��
    % �����lamda-1��������Ͻ��е��㽻��
    % �Բ�ͬ�Ե�Ⱦɫ��������ѡ�񽻲��λ�ý��н���
    Ac = '';
    Bc = '';
    for i=1:size(Kfy,1)
        A = Kfy(i,:);
        B = Kfx(i,:);
        % ���ȷ������λ��
        index = randperm(lamda-1,1);
        % [A;B]
        % �������
        temp = A(index+1:end);
        Ac(i,:) = [A(1:index),B(index+1:end)];
        Bc(i,:) = [B(1:index),temp];
    end
    C = [Ac;Bc];

    % (8)���죺������λ���죨�������������B = Pm*M*L = 0.01*size(C,1)*lamda��
    Kc1 = 0;
    Kc2 = 0;
    for i=1:size(C,1)
        % ����lamda�������
        n = rand(1,lamda);
        % ��lamda��������������ʽ��бȽϣ��ж��Ƿ���ڱ������λ���Ƿ���Ҫ���б��죩
        ps = find(n < Pm);
        % �����ڣ������Ӧ����λ���б��촦��
        if ~isempty(ps)
            % ���л���λ���죨�п���ֻ��һ������λ�����ı䣬Ҳ�п��ܶ������λ�����ı䣩
            Kc1 = Kc1 + 1; % ͳ�Ʊ���Ⱦɫ������
            for k=1:length(ps)
                C(i,:) = [C(i,1:ps(k)-1),dec2bin(~bin2dec(C(i,ps(k)))),C(i,ps(k)+1:end)];
                Kc2 = Kc2 + 1; % ͳ�Ʊ����������
            end
        end
    end
%     sprintf("����Ⱦɫ��������%d",Kc1)
%     sprintf("�������������%d",Kc2)
%     sprintf("�����ʣ�%.3f",Kc2/(size(C,1)*lamda))

    % (9)������Ⱥ������Ⱥ = ��ѡ�е�δ��������ĸ�Ⱥ + ��ѡ���ҷ������䣨���죩�����Ⱥ��
    S = [Kf_;C];
end

% �Ŵ��㷨������Сֵ
f_min = 20 - f_max;
% ͳ�������и�Ԫ�س��ֵĴ�����Ƶ��
J = tabulate(f_min);
% ������Сֵ���ֵ�λ�ã������ж����Сֵ��
J_cnt = J(:,1)==min(J(:,1));
% ������Сֵ
J_min = roundn(double(min(J(:,1))),-4);
% ����õ���Сֵʱ�ı���X��ֵ
Sdecu = Sdecu(J_cnt);


