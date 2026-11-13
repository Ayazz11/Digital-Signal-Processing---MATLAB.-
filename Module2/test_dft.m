clc;
clear;
x = [1 2 0 1 5 6 1 8 5 9 12 10];
N=length(x);
n=(0:N-1);
[X,k]=dft(x,8);
X2=fft(x,8);
k2 = 0:length(X2)-1;
[x_r,n_r]=idft(X);

subplot(4,2,3);
stem(k,real(X));
title("Real part of 8 POINT(M<N) DFT implementation - for Truncated sequence");
xlabel("Bins");
ylabel("DFT samples");


subplot(4,2,4);
stem(k,imag(X));
title("img part of 8 POINT DFT(M<N) implementation - for Truncated sequence");
xlabel("Bins");
ylabel("samples");

subplot(4,2,5);
stem(k2,real(X2));
title("Real part of 8 POINT DFT using inbuilt function fft()");
xlabel("Bins");
ylabel("DFT samples");


subplot(4,2,6);
stem(k2,imag(X2));
title("img part of 8 POINT DFT using inbuilt function fft()");
xlabel("Bins");
ylabel("samples");

subplot(4,2,1);
stem(n,abs(x));
title("original sequence 12 point sequence");
xlabel("time index");
ylabel("samples");

subplot(4,2,7);
stem(n_r,abs(x_r));
title("recovered shorter length sequence from 8 point DFT samples using IDFT");
xlabel("time index");
ylabel("samples");
