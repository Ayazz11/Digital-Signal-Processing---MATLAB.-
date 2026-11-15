% dtft_dft_relation.m - the M-point DFT is M samples of the DTFT at w_k = 2*pi*k/M
% Requires dtft_direct.m and your dft.m ([X,k] = dft(x,M), folding when M < N).
clc; clear; close all;

x = [1 3 2 4 1 2 3 1];               % N = 8, starts at n = 0
N = length(x);
w = linspace(0, 2*pi, 1000);         % dense grid, stands in for the continuous DTFT
Xw = dtft(x, w);

%% 1) DFT samples lie on the DTFT; zero-padding only makes the grid denser
Ms = [8 16 32 64];
figure(1);
for p = 1:4
    M  = Ms(p);
    X  = dft(x, M);
    wk = 2*pi*(0:M-1)/M;             % bin k <-> frequency 2*pi*k/M
    err = max(abs(X - dtft(x, wk)));
    fprintf('M = %2d : max |DFT - DTFT(w_k)| = %.2e\n', M, err);
    subplot(2,2,p);
    plot(w/pi, abs(Xw)); hold on;
    stem(wk/pi, abs(X), 'filled'); hold off;
    xlabel("\omega / \pi"); ylabel("|X|"); grid on;
    title(sprintf("M = %d (N = %d)", M, N));
end

%% 2) M < N: folding keeps the samples ON the DTFT of x; truncation does not
M  = 5;
Xa = dft(x, M);                      % folds x modulo M
Xt = dft(x(1:M), M);                 % keeps only the first M samples
wk = 2*pi*(0:M-1)/M;
Xd = dtft(x, wk);             % DTFT of the FULL sequence at w_k
fprintf('M = %d < N: fold error = %.2e, truncate error = %.2e\n', ...
        M, max(abs(Xa - Xd)), max(abs(Xt - Xd)));
figure(2);
plot(w/pi, abs(Xw)); hold on;
stem(wk/pi, abs(Xa), 'filled');
stem(wk/pi, abs(Xt), 'r');
hold off; grid on;
legend("DTFT of x", "fold (alias)", "truncate");
xlabel("\omega / \pi"); ylabel("|X|");
title(sprintf("M = %d < N = %d", M, N));

%% 3) Magnitude and phase for M = N
M  = N;
X  = dft(x, M);
wk = 2*pi*(0:M-1)/M;
figure(3);
subplot(2,1,1);
plot(w/pi, abs(Xw)); hold on; stem(wk/pi, abs(X), 'filled'); hold off;
xlabel("\omega / \pi"); ylabel("|X|"); title("Magnitude"); grid on;
subplot(2,1,2);
plot(w/pi, angle(Xw)); hold on; stem(wk/pi, angle(X), 'filled'); hold off;
xlabel("\omega / \pi"); ylabel("phase (rad)"); title("Phase"); grid on;

%% 4) For M >= N the DFT samples determine the whole DTFT (interpolation)
%   X(w) = sum_k X[k] * P(w - 2*pi*k/M),   P(t) = (1/M) * sum_{n=0}^{M-1} exp(-j*t*n)
M  = N;
X  = dft(x, M);
Xr = zeros(size(w));
for k = 0:M-1
    t = w - 2*pi*k/M;
    P = zeros(size(w));
    for n = 0:M-1
        P = P + exp(-1j*t*n);
    end
    Xr = Xr + X(k+1) * P / M;
end
fprintf('Reconstruction from %d DFT samples: max |DTFT - rebuilt| = %.2e\n', ...
        M, max(abs(Xw - Xr)));
figure(4);
plot(w/pi, abs(Xw)); hold on;
plot(w/pi, abs(Xr), '--');
stem((2*pi*(0:M-1)/M)/pi, abs(X), 'filled'); hold off; grid on;
legend("DTFT", "rebuilt from DFT samples", "DFT samples");
xlabel("\omega / \pi"); ylabel("|X|");
title("Interpolating the DFT samples gives back the DTFT");