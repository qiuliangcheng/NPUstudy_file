function[y]=over_save(x,h)
M=length(h);N=M+1;
L=N+M-1;
Lx=length(x);
T=ceil(Lx/N);
t=zeros(1,M-1);
x=[x,zeros(1,(T+1)*N-Lx)];
y=zeros(1,(T+1)*N);
for i=0:1:T
    xi=i*N+1;
    xseg=[t,x(xi:xi+N-1)];
    t=xseg(N+1:N+M-1);
    yseg=circular_conv(xseg,h,L);
    y(xi:xi+N-1)=yseg(M:M+N-1);
end
y=y(1:Lx+M-1);
figure(5);
subplot(2,1,2); plot(y);title('重叠保留法计算滤波输出yt的波形');grid;
xlabel('t/s'); ylabel('y(t)');axis([0,1100,-1,1]);