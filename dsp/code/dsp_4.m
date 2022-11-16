
%实验四--第一题
%%
clear all;clc;

N=1600;   %N为信号st的长度。
Fs=10000;T=1/Fs;Tp=N*T; %采样频率Fs=10kHz，Tp为采样时间
t=0:T:(N-1)*T;k=0:N-1;f=k/Tp;
fc1=Fs/40;      %第1路调幅信号的载波频率fc1=250Hz,
fm1=fc1/5;      %第1路调幅信号的调制信号频率fm1=50Hz
fc2=Fs/20;      %第2路调幅信号的载波频率fc2=500Hz
fm2=fc2/5;      %第2路调幅信号的调制信号频率fm2=50Hz
fc3=Fs/10;      %第3路调幅信号的载波频率fc3=1000Hz,
fm3=fc3/5;      %第3路调幅信号的调制信号频率fm3=25Hz
xt1=cos(2*pi*fm1*t).*cos(2*pi*fc1*t); %产生第1路调幅信号
xt2=cos(2*pi*fm2*t).*cos(2*pi*fc2*t); %产生第2路调幅信号
xt3=cos(2*pi*fm3*t).*cos(2*pi*fc3*t); %产生第3路调幅信号
st=xt1+xt2+xt3;         %三路调幅信号相加
fxt=fft(st,N);          %计算信号st的频谱
%====以下为绘图部分，绘制st的时域波形和幅频特性曲线====================
subplot(2,1,1)
plot(t,st);grid;xlabel('t/s');ylabel('s(t)');
axis([0,Tp/8,min(st),max(st)]);title('(a) s(t)的波形')
subplot(2,1,2)
stem(f,abs(fxt)/max(abs(fxt)),'.');grid;title('(b) s(t)的频谱')
axis([0,Fs/5,0,1.2]);
xlabel('f/Hz');ylabel('幅度')


%%
%实验四--第二题
clear all;clc; 
%IIR数字滤波器的设计
Fs=10000;T=1/Fs;   %采样频率
%调用信号产生函数mstg产生由三路抑制载波调幅信号相加构成的复合信号st 
st=mstg;
%低通滤波器设计与实现
fp=280;fs=450;
wp=2*fp/Fs;ws=2*fs/Fs;rp=0.1;rs=60;   %DF指标（低通滤波器的通、阻带边界频）
[N,wp]=ellipord(wp,ws,rp,rs);         %调用ellipord计算椭圆DF阶数N和通带截止频率wp
[B,A]=ellip(N,rp,rs,wp);              %调用ellip计算椭圆带通DF系统函数系数向量B和A
y1t=filter(B,A,st);                   %滤波器的软件实现 
%低通滤波器设计与实现的绘图部分
figure(2);
subplot(3,2,1);myplot(B,A);     %调用绘图函数myplot绘制损耗函数曲线
title('低通滤波器损耗函数曲线');
yt='y_1(t)';
subplot(3,2,2);tplot(y1t,T,yt); %调用绘图函数tplot绘制滤波器输出波形
title('低通滤波器的输出y1t波形');
%带通滤波器设计与实现
fpl=440;fpu=560;fsl=275;fsu=900;
wp=[2*fpl/Fs,2*fpu/Fs];ws=[2*fsl/Fs,2*fsu/Fs];rp=0.1;rs=60; %DF指标（带通滤波器的通、阻带边界频）,wp为3db的截止频率
[N,wp]=ellipord(wp,ws,rp,rs);    %调用ellipord计算椭圆DF阶数N和通带3db截止频率wp
[B,A]=ellip(N,rp,rs,wp);         %调用ellip计算椭圆带通DF系统函数系数向量B和A
y2t=filter(B,A,st);              %滤波器软件实现
%带通滤波器设计与实现绘图部分
figure(2);
subplot(3,2,3);myplot(B,A);     %调用绘图函数myplot绘制损耗函数曲线
title('带通滤波器损耗函数曲线');
yt='y_2(t)';
subplot(3,2,4);tplot(y2t,T,yt); %调用绘图函数tplot绘制滤波器输出波形
title('带通滤波器的输出y2t波形');
%高通滤波器设计与实现
fp=890;fs=600;
wp=2*fp/Fs;ws=2*fs/Fs;rp=0.1;rs=60;   %DF指标（高通滤波器的通、阻带边界频）
[N,wp]=ellipord(wp,ws,rp,rs);         %调用ellipord计算椭圆DF阶数N和通带截止频率wp
[B,A]=ellip(N,rp,rs,wp,'high');       %调用ellip计算椭圆带通DF系统函数系数向量B和A
y3t=filter(B,A,st);                   %滤波器软件实现
%高低通滤波器设计与实现绘图部分
figure(2);
subplot(3,2,5);myplot(B,A);     %调用绘图函数myplot绘制损耗函数曲线
title('高通滤波器损耗函数曲线');
yt='y_3(t)';
subplot(3,2,6);tplot(y3t,T,yt); %调用绘图函数tplot绘制滤波器输出波形
title('高通滤波器的输出y3t波形');


