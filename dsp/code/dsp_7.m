clc;
clear;
%------------产生三个频率的信号------------
N=1000;FS=2000;
T=1/FS;Tp=N*T;t=0:T:(N-1)*T;k=0:N-1;f=k/Tp;
f1=100;f2=500;f3=800;
y1=cos(2*pi*f1*t);
y2=cos(2*pi*f2*t);
y3=cos(2*pi*f3*t);
y=y1+y2+y3;
F=fft(y,N);
subplot(2,1,1);
plot(t,y);
title('产生信号Y的波形');grid;
xlabel('t'); ylabel('y'); 
xlim([0,Tp/10])
subplot(2,1,2);
stem(f,abs(F)/max(abs(F)),'.');
title('y的频谱');grid;
xlabel('f/Hz'); ylabel('幅度'); 
xlim([0,1000]);
%---------设计一个低通滤波器-----------
fp=120;fs=460;wp=2*fp*T;ws=2*fs*T;
wc=(ws+wp)/2;
ND=50;
f=FS/N*k;
hn=fir1(ND-1,wc,blackman(ND));
figure(2);
subplot(2,1,1)
plot(hn);
title('低通滤波器h(n)的波形');grid;
xlabel('n'); ylabel('hn'); 
subplot(2,1,2)
HW=fft(hn,N);
plot(f,20*log10(abs(HW)/max(abs(HW))));
title('低通滤波器h(n)幅频特性曲线');grid;
xlabel('f/Hz'); ylabel('|H(jw)|'); 
axis([0,1000,-100,5]);
%-------------滤除800HZ的频率------------
fp=520;fs=740;wp=2*fp/FS;ws=2*fs/FS;wc=(wp+ws)/2;
Nd=50;
hn2=fir1(Nd-1,wc,blackman(Nd));
k=0:N-1;f=FS/N*k;
Hw=fft(hn2,N);%此低通滤波器的频谱
figure(3);
subplot(2,1,1); plot(hn2);title('低通滤波器h2(n)的波形');grid;
xlabel('n'); ylabel('h2(n)'); axis auto normal;
subplot(2,1,2); plot(f,20*log10(abs(Hw)/max(abs(Hw))));
title('低通滤波器h2(n)幅频特性曲线');grid;
xlabel('f/Hz'); ylabel('|H(jw)|'); axis([0,1000,-100,5]);
%-----------利用重叠保留法和重叠相加法进行卷积运算------------
y1=add(y,hn);
y2=over_save(y,hn2);
N=length(y1);
YK1=fft(y1,N);
YK2=fft(y2,N);
Fs=2000;T=1/Fs;Tp=N*T;
t=0:T:(N-1)*T;k=0:N-1;f=k/Tp;
figure(4)
subplot(2,1,1);stem(f,abs(YK1)/max(abs(YK1)),'.');
title('重叠相加法计算y1(t)的频谱');
xlabel('f/Hz'); ylabel('幅度');axis([0,600,0,1]);grid;
subplot(2,1,2);stem(f,abs(YK2)/max(abs(YK2)),'.');
xlim([0,1000]);xlabel('f/Hz'); ylabel('幅度');
title('重叠保留法计算y2(t)的频谱');



