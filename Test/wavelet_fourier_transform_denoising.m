%% С��������һ�㲽��
% ��1��������������ѡ������С��ĸ�����������ع�С������
% ��2��ȷ��С���任�����ͺ�ά�ȣ�ѡ����ɢդ��С���任��������С���任��������һά�Ļ��Ƕ�ά�ģ�
% ��3��ѡ��ǡ����MATLAB�������źŽ���С���任���Խ��������ʾ�������ʹ���
% ��4������Ա任��������˴�����ѡ��ǡ����MATLAB�������źŽ����ؽ����ع���

snr = 4;
% ���������
init = 2055615866;
% �����������ֵ
[si,xi] = wnoise(1,11,snr,init);
% �������β��źźͺ��������ź�
lev = 5;
xd = wden(xi,'heursure','s','one',lev,'sym8');
figure
subplot(231);
plot(si);
axis([1 2048 -15 15]);
title('ԭʼ�ź�');

subplot(232);
plot(xi);
axis([1 2048 -15 15]);
title('�������ź�');

ssi = fft(si);
ssi = abs(ssi);
xxi =  fft(xi);
absx = abs(xxi);

subplot(233);
plot(ssi);
title('ԭʼ�źŵ�Ƶ��');

subplot(234);
plot(absx);
axis([1 2048 -15 15]);
title('�����źŵ�Ƶ��'); % ���е�ͨ�˲�

indd2 = 200:1800;
xxi(indd2) = zeros(size(indd2));
xden = ifft(xxi); % ���и���Ҷ���任
xden = real(xden);
xden = abs(xden);

subplot(235);
plot(xd);
axis([1 2048 -15 15]);
title('С��ȥ�����ź�');

subplot(236);
plot(xden);
axis([1 2048 -15 15]);
title('����Ҷȥ�����ź�');

%% ��Mallat�㷨����С���׷���
% ��1������Ƶ�ʷֱ�Ϊ5��10�����Ҳ�
clc;
clear;
close all;
f1 = 5; % Ƶ��1
f2 = 10;% Ƶ��2
fs = 2*(f1+f2); % ����Ƶ��
Ts = 1/fs; % �������
N = 12; % ��������
n = 1:N;
y = sin(2*pi*f1*n*Ts) + sin(2*pi*f2*n*Ts); % ���Ҳ����

figure(1)
plot(y);
title('�������Ҳ��ź�');

figure(2)
stem(abs(fft(y)));
title('���ź�Ƶ��');

% ��2����������Ĳ��ν���С���˲����׷���
h = wfilters('db6','l'); % ��ͨ
g = wfilters('db6','h'); % ��ͨ
h = [h, zeros(1, N-length(h))]; % ���㣨Բ�ܾ����������ֱ��ʱ��ڹ۲죩
g = [g, zeros(1, N-length(g))]; % ���㣨Բ�ܾ����������ֱ��ʱ��ڹ۲죩

figure(3)
stem(abs(fft(h)));
title('��ͨ�˲���ͼ');

figure(4)
stem(abs(fft(g)));
title('��ͨ�˲���ͼ');

% ��3��ѡ��Mallat�ֽ��㷨��Բ�ܾ���Ŀ��ٸ���Ҷ�任ʵ�֣����в��δ���
sig1 = ifft(fft(y).*fft(h)); % ��ͨ����Ƶ������
sig2 = ifft(fft(y).*fft(g)); % ��ͨ����Ƶ������

figure(5); % �ź�ͼ
subplot(2,1,1);
plot(real(sig1));
title('�ֽ��ź�1');

subplot(2,1,2);
plot(real(sig2));
title('�ֽ��ź�2');

figure(6); % Ƶ��ͼ
subplot(2,1,1);
stem(abs(fft(sig1)));
title('�ֽ��ź�1Ƶ��');

subplot(2,1,2);
stem(abs(fft(sig2)));
title('�ֽ��ź�2Ƶ��');

% ��4������Mallat�ع��㷨�Ա任�ṹ���д��������ع����ͼ�ν��бȽ�
sig1 = dyaddown(sig1); % 2��ȡ
sig2 = dyaddown(sig2); % 2��ȡ
sig1 = dyadup(sig1); % 2��ֵ
sig2 = dyadup(sig2); % 2��ֵ
sig1 = sig1(1,1:N); % ȥ�����һ��0
sig2 = sig2(1,1:N); % ȥ�����һ��0
hr = h(end:-1:1); % �ع���ͨ
gr = g(end:-1:1); % �ع���ͨ
hr = circshift(hr', 1)'; % λ�õ���Բ������һλ
gr = circshift(gr', 1)'; % λ�õ���Բ������һλ
sig1 = ifft(fft(hr).*fft(sig1)); % ��Ƶ
sig2 = ifft(fft(gr).*fft(sig2)); % ��Ƶ
sig = sig1 + sig2; % Դ�ź�

figure(7);
subplot(2,1,1)
plot(real(sig1));
title('�ع���Ƶ�ź�');

subplot(2,1,2)
plot(real(sig2));
title('�ع���Ƶ�ź�');

figure(8);
subplot(2,1,1)
stem(abs(fft(sig1)));
title('�ع���Ƶ�ź�Ƶ��');

subplot(2,1,2)
stem(abs(fft(sig2)));
title('�ع���Ƶ�ź�Ƶ��');

figure(9);
plot(real(sig),'r','linewidth',2);
hold on;
plot(y);
legend('�ع��ź�','ԭʼ�ź�');
title('�ع��ź���ԭʼ�źűȽ�');
























