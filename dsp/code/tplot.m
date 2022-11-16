function tplot(xn,T,yn)
%时域序列连续曲线绘图函数
%xn:信号数据序列，yn：绘图信号的纵坐标名称（字符串）
%T为采样间隔
n=0:length(xn)-1;t=n*T;
plot(t,xn);
xlabel('t/s');ylabel(yn);
axis([0,t(end),min(xn),1.2*max(xn)])
end 