fs = 48e3;%采样频率
T = 0.05;%信号时长
SNR = -15;%信噪比
f = 15e3;%信号频率
A = 1;%信号振幅
t = 0:1/fs:T-1/fs;%时间向量
sig = A*sin(2*pi*f*t);%单频信号
n = sqrt(A^2/(2*10^(SNR/10)))*randn(size(sig));%高斯白噪声
% 生成1秒的0信号
zeros_1s = zeros(1, fs*T); 
% 将1秒的0信号添加到信号末尾
sig_1 = [sig, zeros_1s];
t = 0:1/fs:T+T-1/fs;%时间向量
figure
plot(t,sig_1);
xlabel('time/s')
ylabel('幅度')
title('fs =200k f=15k')
ylim([-1.5 1.5])
n = sqrt(A^2/(2*10^(SNR/10)))*randn(size(sig_1));
sig_2 = sig_1 + n;
figure
t = 0:1/fs:T+T-1/fs;%时间向量
plot(t,sig_2);
title('信号加噪声');
writeRAF(sig_1,'50',1,fs,0);%将发射信号改为.raf格式文件并保存在存储设备根目录下,.raf格式文件名为数字或字母
% 绘制频谱图
N = length(sig_1);
Y=fft(sig_1)/N;
frequencies = fs*(0:(N/2))/N;% 频率向量
P2 = abs(Y); % 幅度向量
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure;
plot(frequencies, P1);
xlabel('频率 (Hz)');
ylabel('幅度');
title('50K抽样信号频谱图');%可以改变抽样频率绘制