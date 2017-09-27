
keynum = 44;
dur = 1;
fs = 44100;
tt = 0:1/fs:dur;
A = 1;
f0 = 150;
xx = 0;
for i = 1:30
    xx = xx + A*cos(i*f0*2*pi*tt);
end
% soundsc(xx, fs)
ww = 0:pi/256:pi;
Hx = freqz(xx, 1, ww);
% subplot(211)
% plot(abs(Hx))
% pause(2)

OO = filter(B_pz, A_pz, xx);
% soundsc(OO, fs);

%A and B from pole zero place
[Hoo, W] = freqz(B_pz, A_pz, 256);
% subplot(212)
% plot(abs(Hoo))
