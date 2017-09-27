fs = 11250;
tt = 0:1/fs:1;
A = 1;
f0 = 150;
xx = 0;
for i = -30:30
    xx = xx + A*cos(i*f0*2*pi*tt);
end
figure
stem(xx)
