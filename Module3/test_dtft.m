clc;
clear;
nx = -2:2;
%% signal 1
x  = [1 2 0 1 1];
w  = linspace(-2*pi, 2*pi, 2001);
[X,w] = dtft(x,nx,w);

%% IDTFT 
wp = linspace(-pi, pi, 2001);
Xp = dtft(x,nx,wp);
n_r = -4:4;
[x_r,n_r] = idtft(Xp,wp,n_r);

%% Signal 2 
alpha = 0.8;
n2 = 0:199;
x2 = alpha.^n2;
X2 = dtft(x2,n2,w);

%% Signal 3 
wc = pi/4;
Xl = double(abs(wp) < wc);
nl = -20:20;
[xl,nl] = idtft(Xl,wp,nl);

%% plot
subplot(3,2,1);
stem(nx, x, 'filled');
grid on;
title("Input sequence x(n)");
xlabel("time index");
ylabel("Amplitude");

subplot(3,2,2);
plot(w/pi, abs(X));
grid on
title("|X(e^{j\omega})| : continuous, period 2\pi");
xlabel("\omega / \pi");
ylabel("Magnitude");

subplot(3,2,3);
plot(w/pi, angle(X));
grid on
title("Phase of X(e^{j\omega})");
xlabel("\omega / \pi");
ylabel("Phase (rad)");

subplot(3,2,4);
stem(n_r, real(x_r), 'filled');
grid on
title("IDTFT from one period [-\pi, \pi]");
xlabel("time index");
ylabel("Amplitude");

subplot(3,2,5);
plot(w/pi, abs(X2));
grid on
title("\alpha^n u(n), \alpha = 0.8");
xlabel("\omega / \pi");
ylabel("Magnitude");

subplot(3,2,6);
stem(nl, real(xl), 'filled');
grid on
title("IDTFT of ideal low-pass, \omega_c = \pi/4");
xlabel("time index");
ylabel("Amplitude");
