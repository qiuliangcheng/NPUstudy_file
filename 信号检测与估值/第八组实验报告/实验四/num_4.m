%%
%--------------------第四题(1)-----------------------------------
% 参数设置
A = 1; % 信号幅度
T = 1; % 信号周期
fs = 1000; % 采样频率
t = (0:1/fs:T-1/fs)'; % 时间向量
wo = 2 * pi * 50; % 初始频率
dw = 2 * pi * 0.5; % 频率步长
N = 100; % 匹配滤波器数量
SNR=-10;
var = A^2 / (2 * 10^(SNR/10));
n = sqrt(A^2/(2*10^(SNR/10)))*randn(size(t));
% 生成信号
rand_phase = 2 * pi * rand(); % 随机相位
w = wo ;
r=0.1:0.05:10;
o=zeros(length(r),1);
l=zeros(length(r),1);
y=1;
for j=0.1:0.05:10
    A=j;
    o(y)=10 * log10(A^2 / (2*var));
    x = A * sin(w * t+ rand_phase) + n; % 添加噪声
    matched_filter_output = zeros(N, 1);
    for i = 0:N-1
        filtered_signal = sin((wo + i * dw) * (T-t));
        matched_filter_output(i + 1) = sum(abs(hilbert(conv(x , filtered_signal))));
    end
% 频率估计
    [max_val, max_idx] = max(matched_filter_output);
    estimated_freq = (wo + (max_idx - 1) * dw) / (2 * pi);
    l(y)=estimated_freq;
    y=y+1;
end
figure
plot(o,l);
title('信噪比与估计频率')
xlabel('信噪比')
ylabel('频率估计')
%figure
%subplot(311)
%plot(t,x)
%title('信噪比为8 频率为50hz')
%xlabel('时间')
%ylabel('幅度')
%subplot(312)
%A=0.01;
%w=2 * pi * 10;
%x = A * sin(w * t+ rand_phase) + n; % 添加噪声
%plot(t,x)
%title('信噪比为-50 频率为10hz')
%xlabel('时间')
%ylabel('幅度')
%subplot(313)
%A=0.1;
%w=2 * pi * 20;
%x = A * sin(w * t+ rand_phase) + n; % 添加噪声
%plot(t,x)
%title('信噪比为-30 频率为20hz')
%xlabel('时间')
%ylabel('幅度')
%db=10 * log10(A^2 / (2*var))
% 匹配滤波器和最大值选择器

% 显示结果
disp(['True fs: ', num2str(wo/(2*pi)), ' Hz']);
fprintf('Estimated frequency: %.2f Hz\n', estimated_freq);
%%
%-------------------------第四题（2）--------------------------------
% 参数设置   到达时间估计

A = 1; % 幅度
w0 = 2 * pi * 10; % 角频率
t1 = 0.5; % 到达时间
t_min=0.01;
tm=100*t_min;
fs=1000;
dt=1/fs;
T = 5; % 信号持续时间
SNR = 10; % 信噪比 (dB)
var = A^2 / (2 * 10^(SNR/10));

t = 0:dt:T-dt;
% 生成随机相位
m = 2 * pi * rand(1);
% 生成信号
zeros_1s = zeros(1, fs*t1); 
% 将1秒的0信号添加到信号末尾
x = cos(w0 * (t - t1) + m);
t=0:dt:T+t1-dt;
n = sqrt(A^2/(2*10^(SNR/10)))*randn(size(t));
x = [zeros_1s ,x];
x_noisy = A*x+n;
% figure
% subplot(311)
% plot(t,x_noisy);
% title('幅度为1 信噪比为 4 时间延时为0.5s')
% xlabel('时间')
% ylabel('幅度')
% t1=0.1;
% A=0.01;
% SNR=-50;
% zeros_1s = zeros(1, fs*t1); 
% % 将1秒的0信号添加到信号末尾
% t = 0:dt:T-dt;
% x = cos(w0 * (t - t1) + m);
% t=0:dt:T+t1-dt;
% x = [zeros_1s ,x];
% n = sqrt(A^2/(2*10^(SNR/10)))*randn(size(t));
% x_noisy = A*x+n;
% subplot(312)
% plot(t,x_noisy);
% title('幅度为0.01 信噪比为 -50 时间延时为0.1s')
% xlabel('时间')
% ylabel('幅度')
% 
% t1=1;
% A=8;
% SNR=8;
% zeros_1s = zeros(1, fs*t1); 
% % 将1秒的0信号添加到信号末尾
% t = 0:dt:T-dt;
% x = cos(w0 * (t - t1) + m);
% t=0:dt:T+t1-dt;
% n = sqrt(A^2/(2*10^(SNR/10)))*randn(size(t));
% x = [zeros_1s ,x];
% x_noisy = A*x+n;
% subplot(313)
% plot(t,x_noisy);
% title('幅度为8 信噪比为 8 时间延时为1s')
% xlabel('时间')
% ylabel('幅度')
% 添加高斯白噪声
r=0.1:0.05:10;
o=zeros(length(r),1);
l=zeros(length(r),1);
e=1;
for j=0.1:0.05:10
    A=j;
    x_noisy = A*x+n;
    o(e)=10 * log10(A^2 / (2*var));
    % 匹配滤波器
    t2=0:dt:T+tm;
    h = sin(w0 * (T+tm-t2));
    % zeros_2s = zeros(1, fs*(T+tm)); 
    % h = [h,zeros_2s];
    %figure
    %subplot(211)
    %plot(t,x_noisy)
    % t2=0:dt:2*(T+tm);
    %subplot(212)
    %plot(t2,h);
    %title('匹配滤波器信号')
    y = conv(x_noisy, h);

    %figure
    %plot(y,'b')
    %figure
    % 包络检测
    envelope = abs(hilbert(y));
    m_guji=length(envelope)*(T+t1)/length(x);
    m_guji=0:1/fs:m_guji-1/fs;
    plot(m_guji,envelope,'r')
    % 估计到达时间
    if(SNR>=-10)
        window_size = (tm-t1)*900;  % 窗口大小，假设为100
        tolerance = 0.0005;  % 误差容限，假设为2%
        start_idx = 1;  % 起始索引，假设为1
        window1=0;
        while start_idx + window_size <= length(envelope)
            if(envelope(start_idx)>100)
                window1 = envelope(start_idx : start_idx + window_size - 1);  % 获取当前窗口
                if abs(window1(end) - window1(1)) / window1(1) <= tolerance  % 判断窗口最后一个值与第一个值的误差是否小于误差容限
                first_val = window1(1);  % 取出第一个值
                break;
                end
            end
            start_idx = start_idx + 1;  % 窗口向右滑动1个样本
        end
    else
        [max1 ,start_idx]=max(envelope);
    end
    res=m_guji(start_idx);
    t1_estimated = res-T;
    l(e)=t1_estimated;
    e=e+1;
end
% 结果显示
plot(o,l);
title('时延估计与信噪比')
xlabel('信噪比')
ylabel('时间估计')

%%
figure;
subplot(2, 1, 1);
plot(t, x_noisy);
title('Noisy Signal');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2, 1, 2);
plot(m_guji, envelope);
title('Envelope of Matched Filter Output');
xlabel('Time (s)');
ylabel('Amplitude');
disp(['True arrival time: ', num2str(t1), ' s']);
disp(['Estimated arrival time: ', num2str(t1_estimated), ' s']);


%%
%频率与到达时间共同估计
clear;
clc;
% 参数设置
A = 1; % 信号幅度
T = 5; % 信号周期
fs = 1000; % 采样频率
t1=0.5;%真实到达时间
dt=1/fs;
SNR=-10;%信噪比
wo = 2 * pi * 10; % 初始频率
dw = 2 * pi * 0.5; % 频率步长
rand_phase = 2 * pi * rand(); % 随机相位
zeros_1s = zeros(1, fs*t1); 
tm=2;
t=0:dt:T-dt;
sig_x = sin(wo * (t-t1)+ rand_phase);
var = A^2 / (2 * 10^(SNR/10));
t=0:dt:T+t1-dt;
n = sqrt(A^2/(2*10^(SNR/10)))*randn(size(t));
N = 100; % 匹配滤波器数量

% 生成信号
w = wo ;
sig2_x = [zeros_1s ,sig_x];
x = A*sig2_x + n; % 添加噪声
% figure
% subplot(311)
% plot(t,x_noisy);
% title('幅度为1 信噪比为 -10 时间延时为0.5s 频率为 10Hz')
% xlabel('时间')
% ylabel('幅度')
% t1=0.1;
% A=0.01;
% SNR=-50;
% wo=2*pi*20;
% zeros_1s = zeros(1, fs*t1); 
% % 将1秒的0信号添加到信号末尾
% t = 0:dt:T-dt;
% x=sin(wo * (t-t1)+ rand_phase);
% t=0:dt:T+t1-dt;
% x = [zeros_1s ,x];
% n = sqrt(A^2/(2*10^(SNR/10)))*randn(size(t));
% x_noisy = A*x+n;
% subplot(312)
% plot(t,x_noisy);
% title('幅度为0.01 信噪比为-50 时间延时为0.1s 频率为20Hz')
% xlabel('时间')
% ylabel('幅度')
% t1=1;
% A=8;
% SNR=8;
% wo=50
% zeros_1s = zeros(1, fs*t1); 
% % 将1秒的0信号添加到信号末尾
% t = 0:dt:T-dt;
% x=sin(2*pi*wo * (t-t1)+ rand_phase);
% t=0:dt:T+t1-dt;
% n = sqrt(A^2/(2*10^(SNR/10)))*randn(size(t));
% x = [zeros_1s ,x];
% x_noisy = A*x+n;
% subplot(313)
% plot(t,x_noisy);
% title('幅度为8 信噪比为 8 时间延时为1s 频率为50Hz')
% xlabel('时间')
% ylabel('幅度')
A=8;
db=10 * log10(A^2 / (2*var))
x = A*sig2_x + n; % 添加噪声
t_esi = (0:1/fs:T+tm-1/fs')'; % 时间向量
% 匹配滤波器和最大值选择器
matched_filter_output = zeros(N, 1);
for i = 0:N-1
    filtered_signal = sin((wo + i * dw) * (T+tm-t_esi));
    c((i+1),:)=abs(hilbert(conv(x , filtered_signal)));
    matched_filter_output(i + 1) = sum(abs(hilbert(conv(x , filtered_signal))));
end
% 频率估计
[max_val, max_idx] = max(matched_filter_output);
time_esta=c(max_idx,:);
time_guji=length(time_esta)*(T+t1)/length(sig2_x );
time_guji=0:1/fs:time_guji-1/fs;
plot(time_guji,time_esta,'r')
estimated_freq = (wo + (max_idx - 1) * dw) / (2 * pi);
if(SNR>=-10)
    window_size = (tm-t1)*900;  % 窗口大小，假设为100
    tolerance = 0.0005;  % 误差容限，假设为2%
    start_idx = 1;  % 起始索引，假设为1
    window1=0;
    while start_idx + window_size <= length(time_esta)
        if(time_esta(start_idx)>100)
            window1 = time_esta(start_idx : start_idx + window_size - 1);  % 获取当前窗口
            if abs(window1(end) - window1(1)) / window1(1) <= tolerance  % 判断窗口最后一个值与第一个值的误差是否小于误差容限
            first_val = window1(1);  % 取出第一个值
            break;
            end
        end
        start_idx = start_idx + 1;  % 窗口向右滑动1个样本
    end
else
    [max1 ,start_idx]=max(time_esta);
end
res=time_guji(start_idx);
t1_estimated = res-T;
% 显示结果
disp(['True fs: ', num2str(wo/(2*pi)), ' Hz']);
fprintf('Estimated frequency: %.2f Hz\n', estimated_freq);
disp(['True arrival time: ', num2str(t1), ' s']);
disp(['Estimated arrival time: ', num2str(t1_estimated), ' s']);