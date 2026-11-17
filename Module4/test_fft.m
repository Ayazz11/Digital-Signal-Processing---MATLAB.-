clc;
clear;

N = 16;  n = 0:N-1;  k = 0:N-1;
x  = cos(2*pi*3*n/N) + 0.5*sin(2*pi*5*n/N);
X  = myfft(x);
X2 = fft(x);
xr = myifft(X);

subplot(2,2,1);
stem(n, x, 'filled');
grid on
title("Input x(n), N = 16");
xlabel("time index");
ylabel("Amplitude");

subplot(2,2,2);
stem(k, abs(X), 'filled');
hold on
plot(k, abs(X2), 'ro', 'MarkerSize', 9);
hold off
grid on
legend("myfft", "fft");
title("|X(k)|");
xlabel("Bin k");
ylabel("Magnitude");

subplot(2,2,3);
stem(n, real(xr), 'filled');
hold on
plot(n, x, 'ro', 'MarkerSize', 9);
hold off
grid on
legend("myifft(myfft(x))", "x(n)");
title("Recovered sequence");
xlabel("time index");
ylabel("Amplitude");