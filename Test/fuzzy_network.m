% ģ��������ANFIS��anfisedit
% �����ƽ����������⾫�ȱƽ�����һ�����������ܼ��ϵ���������
% ������ƣ�һ������ֵ������Ӧ��һ�������һ�������ֻ�ܵ�һ�����Ƶ�����
clc;
clear;
close all;
%%
% x = (0:0.1:5)';
% y = sin(x)./exp(x/5);
% trnData = [x y];
% numMFs = 6;
% mfType = 'gbellmf';
% epoch_n = 30;
% in_fis = genfis1(trnData, numMFs, mfType);
% out_fis = anfis(trnData, in_fis, 30);
% plot(x,y,x,evalfis(x, out_fis));
% legend('Training Data', 'ANFIS Output');
%%
% fismat = readfis('tipper');
% out = evalfis([2 1;4 9], fismat);
%%
Data = [[34 11 11.9 792 3277, 0];
        [35 15 12   777 3101, 0];
        [34 16 15   793 3206, 0];
        [27 10 14   801 3299, 0];
        [33 9  13.5 820 3743, 0];
        [39 13 11   793 3742, 1];
        [45 19 13   811 3955, 0];
        [36 7  16   806 3701, 1];
        [41 12 11   789 3933, 1];
        [40 13 12   790 3935, 1];
        [41 13 11   789 3935, 1]];
% ��Data�������ݱ��浽originData.mat��
save('Data','Data')
% (1)��������
load('Data.mat')
% (2)����ANFIS�Զ�����һ��FIS�ṹ��Ϊ��ʼFIS
% ��originData�е�ǰ10��������Ϊѵ������(ǰ5��������Ϊx�����һ������Ϊy)
X_train = Data(1:10,1:5);
Y_train = Data(1:10,6);
trainData = [X_train Y_train];
in_format = genfis1(trainData);
% (3)�Գ�ʼFIS(in_format)����ѵ��������������ѵ��200�κ�õ�һ��ѵ���õ�ANFISϵͳ��
[format1, error1, stepsize] = anfis(trainData, in_format, 200);
% (4)�����������ݶ�ѵ���õ�ģ����ϵͳ������֤��
y1 = evalfis(X_train, format1); % ��ѵ����������
X_test = Data(11,1:5);
Y_test = Data(11,6);
y2 = evalfis(X_test, format1); % �ò�����������
plot(1:10,Y_train,'b--',1:10,y1,'k.')
legend('Y:ѵ������', 'y1:��������')
xlabel('�������')
ylabel('ϵͳ���')









