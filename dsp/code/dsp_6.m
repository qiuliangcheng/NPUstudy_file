clear all;close all;
%%
fs = 1000;    %采样频率
f1 = 100;   %细化起始频率
f2 = 200;  %细化结束频率
ND=50'
wp=2*f1/fs;ws=2*f2/fs;
wc=(wp+ws)/2;
h=fir1(ND-1,wc,blackman(ND));
% d = designfilt('lowpassfir','PassbandFrequency',0.25, ...
%          'StopbandFrequency',0.35,'PassbandRipple',0.5, ...
%          'StopbandAttenuation',65,'DesignMethod','kaiserwin');
% h = tf(d);   %系统传递函数
m = 1024;    %变换点数
y = fft(h,m);  % 直接FFT结果

w = exp(-1j*2*pi*(f2-f1)/(m*fs));   %螺旋轮廓点之间的比值
a = exp(1j*2*pi*f1/fs);             %螺旋轮廓起点
z = czt(h,m,w,a);                   %CZT变换
fn = (0:m-1)'/m;
fy = fs*fn;                         %频率
fz = (f2-f1)*fn + f1;
figure(1);
subplot(2,1,1)
plot(fy,abs(y))
xlim([0 500])
legend('FFT')
subplot(2,1,2)
plot(fz,abs(z),'b')
xlim([100 200])
legend('CZT')
xlabel('Frequency (Hz)')
grid on;
figure(2);
plot(fy,abs(y),'r-*',fz,abs(z),'b')
xlim([100 200])
legend('FFT','CZT')
xlabel('Frequency (Hz)')
grid on;
%%
%---------第二题--------------
A0=0.6;W_0=1.2;a_0=pi/3;w_0=2*pi/100;
A=A0*exp(1i*a_0);W=W_0*exp(-1i*w_0);
N=8;M=64;L=128;%给定参数指标
k=0:0.001:M-1;
zk=A.*W.^(-1j*k);
n=0:7;xn=ones(1,8);
A1=power(A,-n);W1=power(W,n.*n/2);
yn1=xn.*A1.*W1;yn=[yn1,zeros(1,L-N)];
yk=fft(yn);%产生y(n)序列并进行频谱分析
n=0:M-1;W2=power(W,-n.*n/2);
n=M:L-1;W3=W.^(-(L-n).*(L-n)/2);
hn=[W2,W3];hk=fft(hn);%产生h(n)序列并进行频谱分析
mk=ifft(yk.*hk);%产生v（k）序列
k=0:L-1;
Wk=W.^(k.^2/2);
k=0;
while(k<M)
    Xzk(k+1)=mk(k+1)*Wk(k+1);
    k=k+1;
end
figure(1)
subplot(1,1,1);stem(abs(mk),'.');axis auto normal;
title('序列m(k)波形图');xlabel('k');ylabel('m(k)');
figure(2)
subplot(2,1,1),stem(abs(Xzk),'.');axis auto;
title('Xzk频谱图');xlabel('k');ylabel('|Xzk|');
subplot(2,1,2),stem(20*log10(abs(Xzk)/max(abs(Xzk))),'.');
title('X(zk)分贝值');xlabel('k');axis auto;
ylabel('dB');
figure(3)
polar(angle(zk),abs(zk),'b');