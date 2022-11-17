t=0:0.01:180;
t1=t*pi/180;
n=length(t);
for i=1:n
 y(i)=(1+sin(t1(i)))/(1-sin(t1(i)));
end
plot(y)