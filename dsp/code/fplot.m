function fplot(B,A)
[H,W]=freqz(B,A,1000);
m=abs(H);
plot(W/pi*5000,20*log10(m/max(m)));
xlabel('f/Hz'); ylabel('·ù¶È/dB'); axis([0,4000,-80,5]);
end