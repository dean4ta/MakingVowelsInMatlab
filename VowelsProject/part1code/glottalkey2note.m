function xx = glottalkey2note(keynum, dur)
fs = 44100;
tt = 0:1/fs:dur;
A = 1;
f0 = 220*(2^((keynum-49)/12));
xx = 0;
for i = 1:10
    xx = xx + A*cos(i*f0*2*pi*tt);
end
soundsc(xx, fs)
end
%frequency response of xx
