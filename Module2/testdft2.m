clc;
clear;
x = [1 2 0 1 5 6 1 8 5 9 12 10];
N = length(x);
M = 8;                                  

[Xf,k]  = dft(x, M);                  
[Xt,k]  = dft(x(1:M), M);              

% DTFT of the ORIGINAL sequence at w_k = 2*pi*k/M
wk = 2*pi*k/M;
Xd = zeros(1, M);
for n = 0:N-1
    Xd = Xd + x(n+1)*exp(-1j*wk*n);
end

% DTFT of the original sequence on a dense grid (for the curve)
w  = linspace(0, 2*pi, 1000);
Xw = zeros(size(w));
for n = 0:N-1
    Xw = Xw + x(n+1)*exp(-1j*w*n);
end

subplot(2,2,1);
stem(k, real(Xd), 'k', 'filled');
hold on
stem(k, real(Xf), 'g');
stem(k, real(Xt), 'r');
title("Real part");
xlabel("Bins");
ylabel("samples");
legend("DTFT samples", "folded DFT", "truncated DFT");
subplot(2,2,2);
stem(k, imag(Xd), 'k', 'filled');
hold on
stem(k, imag(Xf), 'g');
stem(k, imag(Xt), 'r');
title("Imaginary part");
xlabel("Bins");
ylabel("samples");
subplot(2,2,3);
plot(w/pi, abs(Xw));
hold on
stem(wk/pi, abs(Xf), 'g');
stem(wk/pi, abs(Xt), 'r');
title("Magnitude against DTFT of original x");
xlabel("\omega / \pi");
ylabel("|X|");
legend("DTFT of x", "folded DFT", "truncated DFT");
