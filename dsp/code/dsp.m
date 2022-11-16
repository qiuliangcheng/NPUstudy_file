close all; clear all

%内容1:调用filter解差分方程,由系统对u(n)的响应判断稳定性
A=[1,-0.9];
B=[0.05,0.05]; %系统差分方程系数向量B和A 
xln=[1 1 1 1 1 1 1 1 zeros(1,50)]; %产生信号xln=R8n 
x2n=ones(1,128); %产生信号x2n=un 
hn=impz(B,A,58); %求系统单位脉冲响应h(n) 
subplot(3,1,1);y='h(n)';
stem(hn,'.');%调用函数tstem绘图 
xlabel('n')
ylabel('hn')
title('(1)系统单位脉冲响应h(n)');
yln=filter(B,A,xln); %求系统对x1n的响应yln 
subplot(3,1,2);
y='y1(n)';stem(yln,'.'); title('(2)系统对R8(n)的响应yl(n)')
xlabel('n');
ylabel('y1n');
y2n=filter(B,A,x2n); %求系统对x2n的响应y2n 
subplot(3,1,3);
y='y2(n)';
stem(y2n,'.'); 
xlabel('n');
ylabel('y2n');
title('(3)系统对u(n)的响应v2(n)')
%=-===-==%内容2:调用conv函数计算卷积
x1n=[1 1 1 1 1 1 1 1]; %产生信号xln=R8n 
hln=[ones(1,10) zeros(1,10)]; h2n=[1 2.5 2.5 1 zeros(1,10)]; y21n=conv(hln,x1n); y22n=conv(h2n,xln);
figure(2)
subplot(2,2,1);stem(hln,'.');%调用函数tstem绘医 
title('(1)系统单位脉冲响应h1(n)')
subplot(2,2,2);stem(y21n,'.'); title('(2)h1(n)与R8(n)的卷积v21(n)');
xlim([0,25]);ylim([0,8]);
subplot(2,2,3);stem(h2n,'.');%调用函数tstem绘图 
ylim([0,3]);title('(3)系统单位脉冲响应h2(n)')
subplot(2,2,4);stem(y22n,'.'); 
xlim([0,25]);ylim([0,8]);
title('(4)h2(n)与R8(n)的卷积y22(n)')
%一======%内容3:谐振器分析
un=ones(1,256);%产生信号ur 
n=0:255;
xsin=sin(0.014*n)+sin(0.4*n); %产生正弦信号
A=[1,-1.8237,0.9801];
B=[1/100.49,0,-1/100.49];
y31n=filter(B,A,un);
y32n=filter(B,A,xsin);
figure(3);
subplot(2,1,1);
y='y31(n)';
stem(y31n,'.');
title('输入信号位u(n)时的波形');
xlim([1,256]);
subplot(2,1,2);
y='y32(n)';
stem(y32n,'.');
title('输入信号位x(n)时的波形');xlim([1,256]);