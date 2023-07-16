
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
n = sqrt(A^2/(2*10^(SNR/10)))*randn(size(sig_1));
sig_2 = sig_1 + n;%信号+噪声
mf_1 = (fliplr(sig_1));
figure
plot(mf_1);
delay_sec = T+T;
% 使用 circshift 函数向左时移 y 向量
delay_samples = round(delay_sec * fs); 
mf= circshift(mf_1 , [-1, -delay_samples]);
% 绘制时移后的 y 向量
t_delay = t-delay_samples/fs;
% 匹配滤波器的冲激响应
figure(3);
subplot(2,1,1);
plot(t,sig_1)
title('原信号')
xlabel('time/s')
ylabel('幅度')
subplot(2,1,2);
plot(t_delay,mf);
title('MF时域图')
xlabel('time/s')
ylabel('幅度')
% subplot(2,1,2);
% plot(t,sig_2);
% title('加了噪声的信号')
N = length(mf);
MF = fft(mf, N);
%滤波器特性
SNR = -15;%信噪比
f = 5e3;%信号频率
A = 1;%信号振幅
t = 0:1/fs:T-1/fs;%时间向量
sig_w = A*sin(2*pi*f*t);%单频信号
% 生成1秒的0信号
zeros_1s = zeros(1, fs*T);
sig_w_cha=[sig_w, zeros_1s];
figure
plot(sig_w_cha)
wres=conv(sig_w_cha,mf);
t = 0:1/fs:T+T-1/fs;
Tw=T+T;
m=length(wres)*Tw/(length(sig_w_cha ));
m=0:1/fs:m-1/fs;
figure
plot(m,wres)
title('频率通过匹配滤波器')
% 将1秒的0信号添加到信号末尾
%延时
delay_times=0.05;
zeros_delays = zeros(1, delay_times*fs); 
sig_1_delay = [zeros_delays ,sig, zeros_1s];
n = sqrt(A^2/(2*10^(SNR/10)))*randn(size(sig_1_delay));
sig_2_delay = sig_1_delay;
t_dealy = 0:1/fs:T+T+delay_times-1/fs;
T_delay=T+T+delay_times;
delay_final_res=conv(sig_2_delay,mf);
m=length(delay_final_res)*T_delay/(length(sig_2_delay ));
m=0:1/fs:m-1/fs;
figure
subplot(3,1,1)
plot(t,sig_1);
title('没有延时的信号')
subplot(3,1,2)
plot(t_dealy,sig_2_delay);
title('延时信号')
xlabel('time/s')
ylabel('幅度')
subplot(3,1,3)
plot(m,delay_final_res)
title('延时信号通过匹配滤波器')
xlabel('time/s')
ylabel('幅度')
%正常信号
final_res=conv(sig_2,mf);
m_ori=length(final_res)*T/length(sig);
m_ori=0:1/fs:m_ori-1/fs;
frequencies = -fs/2:fs/length(mf):fs/2-fs/length(mf); % 频率向量
amplitudes = abs(MF); % 幅度向量
figure;
subplot(2,1,1)
plot(frequencies, amplitudes);
xlim([0 fs/2-fs/length(mf)])
xlabel('频率 (Hz)');
ylabel('幅度');
title('匹配滤波器频率响应');
subplot(2,1,2);
plot(m_ori,final_res);
title('信噪比-20DB信号通过匹配滤波器')
xlabel('time/s')
ylabel('幅度')