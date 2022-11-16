function xt=xtg(N)
%信号x(t)产生,并显示信号的幅频特性曲线
%xt=xtg(N) 产生一个长度为N,有加性高频噪声的单频调幅信号xt,采样频率Fs=1000Hz
%载波频率fc=Fs/10=100Hz,调制正弦波频率f0=fc/10=10Hz.
Fs=1000;T=1/Fs;Tp=N*T;
t=0:T:(N-1)*T;
fc=Fs/10;f0=fc/10; %载波频率fc=Fs/10，单频调制信号频率为f0=Fc/10;
mt=cos(2*pi*f0*t);  %产生单频正弦波调制信号mt，频率为f0
ct=cos(2*pi*fc*t);  %产生载波正弦波信号ct，频率为fc
xt=mt.*ct;          %相乘产生单频调制信号xt
nt=2*rand(1,N)-1;   %产生随机噪声nt
%=======设计高通滤波器hn,用于滤除噪声nt中的低频成分,生成高通噪声=======
fp=150; fs=200;Rp=0.1;As=70;    % 滤波器指标
fb=[fp,fs];m=[0,1];           % 计算remezord函数所需参数f,m,dev
dev=[10^(-As/20),(10^(Rp/20)-1)/(10^(Rp/20)+1)];
[n,fo,mo,W]=remezord(fb,m,dev,Fs);  % 确定remez函数所需参数
hn=remez(n,fo,mo,W);      % 调用remez函数进行设计,用于滤除噪声nt中的低频成分
yt=filter(hn,1,10*nt);      %滤除随机噪声中低频成分，生成高通噪声yt
 
%================================================================
xt=xt+yt;           %噪声加信号
fst=fft(xt,N);k=0:N-1;f=k/Tp;
figure(1)
subplot(2,1,1);plot(t,xt);grid;
xlabel('t/s');ylabel('x(t)');
axis([0,Tp/5,min(xt),max(xt)]);
title('(1) 具有加性噪声的xt')
subplot(2,1,2);plot(f,abs(fst)/max(abs(fst)));grid;
xlabel('f/Hz');ylabel('幅度')
axis([0,Fs/2,0,1.2]);
title('(2) 信号xt的频谱')
end