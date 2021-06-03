clc;clear;
rng default
Fs = 1000;
t = 0:1/Fs:1-1/Fs;
x = cos(2*pi*100*t) + randn(size(t));
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(N*2*pi)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
p = periodogram(x,rectwin(length(x)),length(x));
