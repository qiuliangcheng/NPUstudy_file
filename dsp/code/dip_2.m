clc;
clear;

%================时域采样验证==============
Tp=64/1000; %观察时间Tp=64微秒 
%产生M长采样序列x(n)
Fs=1000;T=1/Fs;
M=Tp*Fs;n=0:M-1;
A=444.128;alph=pi*50*2^0.5;omega=pi*50*2^0.5; xnt=A*exp(-alph*n*T).*sin(omega*n*T); Xk=T*fft(xnt,M);%M点FFT[xnt)] yn='xa(nT)';
figure(1);
subplot(3,2,1);
stem(xnt,'.'); %调用自编绘图函数tstem绘制序列图 
box on;
title('(a)Fs=1000Hz');
k=0:M-1;
fk=k/Tp;
subplot(3,2,2);plot(fk,abs(Xk));title('(a)T*FT[xa(nT)],Fs=1000Hz'); xlabel('f(Hz)');ylabel('幅度');
axis([0,Fs,0,1.2*max(abs(Xk))])

%Fs=300Hz
Fs=300;
T=1/Fs;
M=Tp*Fs;n=0:M-1;
A=444.128;alph=pi*50*2^0.5;omega=pi*50*2^0.5; xnt=A*exp(-alph*n*T).*sin(omega*n*T); Xk=T*fft(xnt,M);%M点FFT[xnt)] yn='xa(nT)';
figure(1);
subplot(3,2,3);
stem(xnt,'.'); %调用自编绘图函数tstem绘制序列图 
box on;
title('(b)Fs=300Hz');
k=0:M-1;
fk=k/Tp;
subplot(3,2,4);plot(fk,abs(Xk));title('(b)T*FT[xa(nT)],Fs=300Hz'); xlabel('f(Hz)');ylabel('幅度');
axis([0,Fs,0,1.2*max(abs(Xk))])

%和Fs=200Hz
Fs=200;
T=1/Fs;
M=Tp*Fs;n=0:M-1;
A=444.128;alph=pi*50*2^0.5;omega=pi*50*2^0.5; xnt=A*exp(-alph*n*T).*sin(omega*n*T); Xk=T*fft(xnt,M);%M点FFT[xnt)] yn='xa(nT)';
figure(1);
subplot(3,2,5);
stem(xnt,'.'); %调用自编绘图函数tstem绘制序列图 
box on;
title('(c)Fs=300Hz');
k=0:M-1;
fk=k/Tp;
subplot(3,2,6);plot(fk,abs(Xk));title('(c)T*FT[xa(nT)],Fs=200Hz'); xlabel('f(Hz)');ylabel('幅度');
axis([0,Fs,0,1.2*max(abs(Xk))])



%===================频率采样验证
figure(1);
M=27;n=0:M-1;
xa=1:ceil(M/2); xb=ceil(M/2)-1:-1:1; xn=[xa,xb];
Xk=fft(xn,1024);
X32k=fft(xn,32);
x32n=ifft(X32k);
X16k=X32k(1:2:32);
x16n=ifft(X16k,16);
subplot(1,2,2); stem(n,xn,'.');
title('(1)序列x(n)'); 
xlabel('n'); ylabel('x(n)'); axis([0,32,0,20]);
k=0:1023; wk=2*k/1024;
subplot(1,2,1); plot(wk,abs(Xk)); 
title('(1)FFT[x(n)]'); 
xlabel('\omega/\pi'); ylabel('|X(e^j^\omega)|');
axis([0,1,0,200]);
k=0:32-1;
figure(2)
subplot(1,2,1); stem(k,abs(X32k),'.');
 title('(2)32点频率采样'); 
xlabel('k'); ylabel('|X_3_2(k)|');
axis([0,16,0,200]);
n1=0:32-1;
subplot(1,2,2); stem(n1,x32n,'.'); 
title('(2)32点IDFT[X32(k)]'); 
xlabel('n'); ylabel('x_3_2(n)'); 
axis([0,32,0,20]);
k=0:32/2-1;
figure(3)
subplot(1,2,1); stem(k,abs(X16k),'.');
title('(3)16点频率采样'); 
xlabel('k'); ylabel('|X_1_6(k)|');
 axis([0,8,0,200]);
n2=0:32/2-1;
subplot(1,2,2);stem(n2,x16n,'.'); 
title('(3)16点IDFT[X16(k)]'); 
xlabel('n'); ylabel('x_1_6(n)'); axis([0,32,0,20]);