%% 小波分析的一般步骤
% （1）根据问题需求，选择或设计小波母函数，及其重构小波函数
% （2）确定小波变换的类型和维度（选择离散栅格小波变换还是序列小波变换，问题是一维的还是多维的）
% （3）选择恰当的MATLAB函数对信号进行小波变换，对结果进行显示、分析和处理
% （4）如果对变换结果进行了处理，可选择恰当的MATLAB函数对信号进行重建（重构）

snr = 4;
% 设置信噪比
init = 2055615866;
% 设置随机数初值
[si,xi] = wnoise(1,11,snr,init);
% 产生矩形波信号和含白噪声信号
lev = 5;
xd = wden(xi,'heursure','s','one',lev,'sym8');
figure
subplot(231);
plot(si);
axis([1 2048 -15 15]);
title('原始信号');

subplot(232);
plot(xi);
axis([1 2048 -15 15]);
title('含噪声信号');

ssi = fft(si);
ssi = abs(ssi);
xxi =  fft(xi);
absx = abs(xxi);

subplot(233);
plot(ssi);
title('原始信号的频谱');

subplot(234);
plot(absx);
axis([1 2048 -15 15]);
title('含噪信号的频谱'); % 进行低通滤波

indd2 = 200:1800;
xxi(indd2) = zeros(size(indd2));
xden = ifft(xxi); % 进行傅里叶反变换
xden = real(xden);
xden = abs(xden);

subplot(235);
plot(xd);
axis([1 2048 -15 15]);
title('小波去燥后的信号');

subplot(236);
plot(xden);
axis([1 2048 -15 15]);
title('傅里叶去燥后的信号');

%% 用Mallat算法进行小波谱分析
% （1）定义频率分别为5和10的正弦波
clc;
clear;
close all;
f1 = 5; % 频率1
f2 = 10;% 频率2
fs = 2*(f1+f2); % 采样频率
Ts = 1/fs; % 采样间隔
N = 12; % 采样点数
n = 1:N;
y = sin(2*pi*f1*n*Ts) + sin(2*pi*f2*n*Ts); % 正弦波混合

figure(1)
plot(y);
title('两个正弦波信号');

figure(2)
stem(abs(fft(y)));
title('两信号频谱');

% （2）对所定义的波形进行小波滤波器谱分析
h = wfilters('db6','l'); % 低通
g = wfilters('db6','h'); % 高通
h = [h, zeros(1, N-length(h))]; % 补零（圆周卷积，且增大分辨率便于观察）
g = [g, zeros(1, N-length(g))]; % 补零（圆周卷积，且增大分辨率便于观察）

figure(3)
stem(abs(fft(h)));
title('低通滤波器图');

figure(4)
stem(abs(fft(g)));
title('高通滤波器图');

% （3）选择Mallat分解算法（圆周卷积的快速傅里叶变换实现）进行波形处理
sig1 = ifft(fft(y).*fft(h)); % 低通（低频分量）
sig2 = ifft(fft(y).*fft(g)); % 高通（高频分量）

figure(5); % 信号图
subplot(2,1,1);
plot(real(sig1));
title('分解信号1');

subplot(2,1,2);
plot(real(sig2));
title('分解信号2');

figure(6); % 频谱图
subplot(2,1,1);
stem(abs(fft(sig1)));
title('分解信号1频谱');

subplot(2,1,2);
stem(abs(fft(sig2)));
title('分解信号2频谱');

% （4）利用Mallat重构算法对变换结构进行处理，并对重构后的图形进行比较
sig1 = dyaddown(sig1); % 2抽取
sig2 = dyaddown(sig2); % 2抽取
sig1 = dyadup(sig1); % 2插值
sig2 = dyadup(sig2); % 2插值
sig1 = sig1(1,1:N); % 去掉最后一个0
sig2 = sig2(1,1:N); % 去掉最后一个0
hr = h(end:-1:1); % 重构低通
gr = g(end:-1:1); % 重构高通
hr = circshift(hr', 1)'; % 位置调整圆周右移一位
gr = circshift(gr', 1)'; % 位置调整圆周右移一位
sig1 = ifft(fft(hr).*fft(sig1)); % 低频
sig2 = ifft(fft(gr).*fft(sig2)); % 高频
sig = sig1 + sig2; % 源信号

figure(7);
subplot(2,1,1)
plot(real(sig1));
title('重构低频信号');

subplot(2,1,2)
plot(real(sig2));
title('重构高频信号');

figure(8);
subplot(2,1,1)
stem(abs(fft(sig1)));
title('重构低频信号频谱');

subplot(2,1,2)
stem(abs(fft(sig2)));
title('重构高频信号频谱');

figure(9);
plot(real(sig),'r','linewidth',2);
hold on;
plot(y);
legend('重构信号','原始信号');
title('重构信号与原始信号比较');
























