%实验五--第一题
close all;clc;
%FIR数字滤波器及软件实现 
%调用 xtg 产生信号 xt, xt 长度 N=1000,并显示 xt 及其频谱
N=1000;xt=xtg(N);
fp=120; fs=150;Rp=0.1;As=60;Fs=1000;    %输入给定指标
%(1)用窗函数法设计滤波器
wc=(fp+fs)/Fs;      %理想低通滤波器截止频率(关于pi归一化）
B=2*pi*(fs-fp)/Fs;  %过渡带宽度指标
Nb=ceil(11*pi/B);   %blackman窗的长度N
hn=fir1(Nb-1,wc,blackman(Nb));
Hw=abs(fft(hn,1024));   % 求设计的滤波器频率特性
ywt=fftfilt(hn,xt,N);   %调用函数fftfilt对xt滤波
%窗函数法设计法的绘图部分（滤波器损耗函数，滤波器输出信号波形）
f=[0:1023]*Fs/1024;
figure(2)
subplot(2,1,1)
plot(f,20*log10(Hw/max(Hw)));grid;
title('(3) 窗函数低通滤波器幅频特性')
xlabel('f/HZ');ylabel('幅度');
axis([0,Fs/2,-120,20]);
t=[0:N-1]/Fs;Tp=N/Fs; 
subplot(2,1,2)
plot(t,ywt);grid;
axis([0,Tp/2,-1,1]);
xlabel('t/s');ylabel('y_w(t)');
title('(4) 滤波噪声后的信号波形') 
%(2)用等波纹最佳逼近法设计滤波器
fb=[fp,fs];m=[1,0];         % 确定remezord函数所需参数f,m,dev
dev=[(10^(Rp/20)-1)/(10^(Rp/20)+1),10^(-As/20)];
[Ne,fo,mo,W]=remezord(fb,m,dev,Fs); % 确定remez函数所需参数
hn=remez(Ne,fo,mo,W);       % 调用remez函数进行设计
Hw=abs(fft(hn,1024));       % 求设计的滤波器频率特性
yet=fftfilt(hn,xt,N);       % 调用函数fftfilt对xt滤波
%等波纹设计法的绘图部分（滤波器损耗函数，滤波器输出信号yw(nT)波形）
f=[0:1023]*Fs/1024;
figure(3)
subplot(2,1,1)
plot(f,20*log10(Hw/max(Hw)));grid;
title('(5) 等波纹最佳逼近低通滤波器幅频特性')
axis([0,Fs/2,-80,10]);
xlabel('f/HZ');ylabel('幅度')
t=[0:N-1]/Fs;Tp=N/Fs; 
subplot(2,1,2);plot(t,yet);grid;
axis([0,Tp/2,-1,1]);
xlabel('t/s');ylabel('y_e(t)');
title('(6) 滤波噪声后的信号波形') 
