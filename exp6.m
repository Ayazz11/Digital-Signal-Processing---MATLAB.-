clear;
f=3; %signal frequency = 3hz.i.e 3 cycles per sec
fs=60; % sampling freq = 40hz. i.e 40 samples per sec
Ts=1/fs;
N=60; % number of samples 
n=0:1:N-1;
t=n.*Ts;
x_n=1.*sin(2.*pi.*f.*t);
subplot(2,2,1);
stem(t,x_n);
xlabel('Time (s)');
ylabel('Amplitude');
title('Discrete Sine Wave');
Y=fft(x_n);
n1=0:N-1; % N/2 + 1 point DFT as N is even 
f1=(n1.*fs)/N;
subplot(2,2,2);
stem(f1,abs(Y));
