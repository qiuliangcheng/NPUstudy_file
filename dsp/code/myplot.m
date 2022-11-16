function myplot(B,A)
%时域离散系统损耗函数绘图
%B为系统函数分子多项式系数向量
%A为系统函数分母多项式系数向量
[H,W]=freqz(B,A,1000);  %freqz函数是用来求离散系统频率响应特性
m=abs(H);        %取幅度值实部
plot(W/pi,20*log10(m/max(m)));grid on;
xlabel('\omega/\pi');ylabel('幅度（dB）')
axis([0,1,-80,5]);
end