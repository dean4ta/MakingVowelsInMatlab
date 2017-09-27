w = 0:pi/100:pi;
z = exp(j*w);
H1 = (z-exp(j*pi/4))./z;
absH1 = abs(H1);
H2 = (z-exp(j*7*pi/8))./z;
absH2 = abs(H2);
subplot(211)
plot(w,absH1)
subplot(212)
plot(w,absH2)