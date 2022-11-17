function[y]=add(x,h)
M=length(h);N=M+1;
L=N+M-1;
Lx=length(x);
T=ceil(Lx/N);
t=zeros(1,M-1);
x=[x,zeros(1,(T+1)*N-Lx)];
y=zeros(1,(T+1)*N);
for i=0:1:T
    xi=i*N+1;
    xseg=x(xi:xi+N-1);
    yseg=circular_conv(xseg,h,L);
    yseg(1:M-1)=yseg(1:M-1)+t(1:M-1);
    t(1:M-1)=yseg(N+1:L);
    y(xi:xi+N-1)=yseg(1:N);
end
y=y(1:Lx+M-1);
figure(5);
subplot(2,1,1); plot(y);title('重叠相加法计算滤波输出yt的波形');grid;
xlabel('t/s'); ylabel('y(t)'); axis([0,1100,-1,1]);