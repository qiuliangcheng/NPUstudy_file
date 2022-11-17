num=[5*20 5];
den=conv(conv(conv([1,1],[0.25 1]),[166.7 1]),[1 0 ]);
sys=tf(num,den);
bode(sys);
grid


