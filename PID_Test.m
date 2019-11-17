% PID控制算法
clc;clear;close all;

% 用户设定值Sv
% 传感器测量值Pv
% 根据其差异进行PWM调节（快速性，低超调，稳态时间端）
% 考虑器件的惯性作用：控制信号改变时，执行器件的状态不会立即发生改变
% 进行模糊控制：不对某一个值进行判断，而是对一个区间进行判断
% PID控制系数（表示对这三个部分的信任程度）：比例系数P（当前值）；积分系数I（过去值：需设定存储容量大小）；微分系数D（预测值）
% 基于偏差（Sv-Pv）的控制策略
% 积分部分：开机以来所有时刻采样到的数据？
% 微分部分：最近几次的偏差？而不是对将来状态的预测？
% 采样间隔：每隔多久读取一次数据
% Er = Sv - Pv：（>0控制未达标；=0控制正好合适；<0控制超标）
% 比例控制：pout=Kp*Er；偏差为0时，不输出控制信号，执行器件处于失控状态
% 比例控制添加常数：pout=Kp*Er + pout0
% 积分控制：（偏差序列：Er(1),Er(2),...,Er(k-1),Er(k)）
% 偏差序列求和：Sk = Er(1) + Er(2) + ... + Er(k-1) + Er(k)（等权重求和；不等权重求和（降序权重））
% Sk>0：过去大多数时间未达标（加强控制信号，控制幅度与Sk大小有关）；Sk<0：过去大多数时间超标（减弱控制信号）；
% 积分控制：Iout = Ki*Sk（Ki为积分系数），当Sk为0时，历史情况会干扰当前状态：Iout = Ki*Sk + iout0
% 微分控制：偏差的变化率：Er(2)-Er(1),Er(3)-Er(2),...,Er(k-1)-Er(k-2),Er(k)-Er(k-1)
% 最近两次的偏差变化率：Dk = Er(k)-Er(k-1)
% Dk>0：偏差变化率出现增大趋势；Dk=0：偏差变化率维持不变；Dk<0：偏差变化率出现减小趋势
% 误差变化趋势具有继承不变性
% 微分控制：Dout = Kd*Dk + dout0（当Dk为0时，dout0不为0，Dout不为0） 

num = 1;
denom = [1 1 1];
t = 0:0.1:30;

Gp = tf(num, denom);
H = 1;

M = feedback(Gp, H);
step(M,t,'r')
P = step(M,t,'r');

% 无阻尼固有频率：wn 
wn = 1;
% 阻尼比：kesai（[0,1]:欠阻尼；0:无阻尼；1;临界阻尼；>1:过阻尼）
kesai = 0.5;
% 有阻尼固有频率：wd 
wd = wn*sqrt(1 - kesai^2);
% （1）上升时间（第一次达到稳态所需要的时间）：tr
% kesai<0有效
atan(sqrt(1-kesai^2)/kesai)
tr = (pi - atan(sqrt(1-kesai^2)/kesai))/wd;
% （2）峰值时间：tp
tp = pi/wd
% （3）最大超调量：Mp
Mp = exp(-kesai*pi/sqrt(1-kesai^2))
% （4）调整时间：ts
ts = 4/(kesai*wn)
% （5）振荡次数：N
N = ceil(2*sqrt(1-kesai^2)/(pi*kesai))

% axis([0 50 0 2])
% Temp1 = step(M,'r');
% T = find(abs(Temp1-Temp1(end))<1e-1,1);
% T1 = 30*T/length(Temp1)
hold on

%%
Kp = 0.5;
Ki = 0.3;
Kd = 0.2;

Gc = pid(Kp, Ki, Kd);

Mc = feedback(Gc*Gp, H);
step(Mc,t,'b')
% axis([0 50 0 2])
% Temp2 = step(Mc,'b');
% % 超调量
% P2 = max(Temp2) - Temp2(end);
% % 首次达到稳态的时间
% T = find(abs(Temp2-Temp2(end))<1e-2,1);
% T2 = 30*T/length(Temp2)
grid on

% 遗传算法求解最优PID控制参数
% 二阶欠阻尼系统单位阶跃响应性能指标
% （1）上升时间：tr
% （2）峰值时间：tp
% （3）最大超调量：Mp
% （4）调整时间：ts
% （5）振荡次数：N


















