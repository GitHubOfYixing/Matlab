% 模糊工具箱ANFIS：anfisedit
% 函数逼近：能以任意精度逼近任意一个定义在致密集上的连续函数
% 解耦控制：一个输入值控制相应的一个输出，一个输出又只受到一个控制的作用
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
% 将Data变量数据保存到originData.mat中
save('Data','Data')
% (1)加载数据
load('Data.mat')
% (2)利用ANFIS自动生成一个FIS结构作为初始FIS
% 将originData中的前10行数据作为训练样本(前5个数据作为x，最后一个数据为y)
X_train = Data(1:10,1:5);
Y_train = Data(1:10,6);
trainData = [X_train Y_train];
in_format = genfis1(trainData);
% (3)对初始FIS(in_format)进行训练。对样本数据训练200次后得到一个训练好的ANFIS系统。
[format1, error1, stepsize] = anfis(trainData, in_format, 200);
% (4)运用评价数据对训练好的模糊神经系统进行验证。
y1 = evalfis(X_train, format1); % 用训练样本仿真
X_test = Data(11,1:5);
Y_test = Data(11,6);
y2 = evalfis(X_test, format1); % 用测试样本仿真
plot(1:10,Y_train,'b--',1:10,y1,'k.')
legend('Y:训练样本', 'y1:测试样本')
xlabel('样本序号')
ylabel('系统输出')









