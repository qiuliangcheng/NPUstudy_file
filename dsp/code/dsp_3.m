% %��FFT���ź���Ƶ�׷��� 
%%
clear all;clc;
M=8;
N=16;
n=0:M;
x1n=ones(1,4);
xa=1:M/2;xb=M/2:-1:1;
x2n=[xa,xb];
xc=M/2:-1:1;xd=1:M/2;
x3n=[xc,xd];
X1k8=fft(x1n,M);
X2k8=fft(x2n,M);
X3k8=fft(x3n,M);
X1k16=fft(x1n,N);   
X2k16=fft(x2n,N);
X3k16=fft(x3n,N);
figure(1);
k=0:M-1;
subplot(1,2,1);stem(k,abs(X1k8),'.');
title('(1)x1n��8��DFT');
xlabel('k');ylabel('X1K-8');
k=0:N-1;
subplot(1,2,2);stem(k,abs(X1k16),'.');
title('(2)x1n��16��DFT');
xlabel('k');ylabel('X1K-16');
k=0:M-1;
figure(2)
subplot(1,2,1);stem(k,abs(X2k8),'.');
title('(3)x2n��8��DFT');
xlabel('k');ylabel('X2K-8');
k=0:N-1;
subplot(1,2,2);stem(k,abs(X2k16),'.');
title('(4)x2n��16��DFT');
xlabel('k');ylabel('X2K-16');
k=0:M-1;
figure(3)
subplot(1,2,1);stem(k,abs(X3k8),'.');
title('(5)x3n��8��DFT');
xlabel('k');ylabel('X3K-8');
k=0:N-1;
subplot(1,2,2);stem(k,abs(X3k16),'.');
title('(6)x3n��16��DFT');
xlabel('k');ylabel('X3K-16');

%%
%ʵ����--�ڶ���
clear all;clc;
N=16;n=0:N-1;
x4n=cos(pi*n/4);
x5n=cos(pi*n/4)+cos(pi*n/8);
N=8;n=0:N-1;
X4k8=fft(x4n,8);
X5k8=fft(x5n,8);
figure(1)
k=0:N-1;subplot(2,1,1);stem(k,abs(X4k8),'.');
title('(1)x4n��8��DFT');
xlabel('k');ylabel('X4k-8');
figure(2)
subplot(2,1,1);stem(k,abs(X5k8),'.');
title('(3)x5n��8��DFT');
xlabel('k');ylabel('X5k-8');
N=16;n=0:N-1;
X4k16=fft(x4n,16);
X5k16=fft(x5n,16);
figure(1)
k=0:N-1;subplot(2,1,2);stem(k,abs(X4k16),'.');
title('(2)x4n��16��DFT');
xlabel('k');ylabel('X4k-16');
figure(2)
subplot(2,1,2);stem(k,abs(X5k16),'.');
title('(4)x5n��16��DFT');
xlabel('k');ylabel('X5k-16');
%%

%ʵ����--������
clear all;clc;
Fs=64;T=1/Fs;
N=16;n=0:N-1;
Tp=N*T;F=1/Tp;
k=-N/2:N/2-1;fk=k*F;
x6aT=cos(8*pi*n*T)+cos(16*pi*n*T)+cos(20*pi*n*T);
X6k16=fft(x6aT);
X6k16=fftshift(X6k16);
subplot(3,1,1);stem(k,abs(X6k16),'.');
title('(1x6n��16��DFT');
xlabel('k');ylabel('X6k-16');
N=32;n=0:N-1;
Tp=N*T;F=1/Tp;
k=-N/2:N/2-1;fk=k*F;
x6aT=cos(8*pi*n*T)+cos(16*pi*n*T)+cos(20*pi*n*T);
X6k32=fft(x6aT);
X6k32=fftshift(X6k32);
subplot(3,1,2);stem(k,abs(X6k32),'.');
title('(2)x6n��32��DFT');
xlabel('k');ylabel('X6k-32');
 N=64;n=0:N-1;
Tp=N*T;F=1/Tp;
k=-N/2:N/2-1;fk=k*F;
x6aT=cos(8*pi*n*T)+cos(16*pi*n*T)+cos(20*pi*n*T);
X6k64=fft(x6aT);
X6k64=fftshift(X6k64);
subplot(3,1,3);stem(k,abs(X6k64),'.');
title('(3)x6n��64��DFT');
xlabel('k');ylabel('X6k-64');