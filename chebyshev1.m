
clc;
clear;
close all;

fp = input("Enter pass band cutoff frequency (Hz): ");
fs = input("Enter stop band cutoff frequency (Hz): ");

wp = 2*pi*fp;
ws = 2*pi*fs;

Gp = input("Enter pass band attenuation in dB: ");
Gs = input("Enter stop band attenuation in dB: ");

[n,wc] = cheb1ord(wp,ws,Gp,Gs,'s');

[b,a] = cheby1(n,Gp,wc,'s');

omega = linspace(0,12000*pi,2000);

h = freqs(b,a,omega);

gain = 20*log10(abs(h));

% Plot
figure;
plot(omega/(2*pi),gain,'b','LineWidth',1.5);
grid on;
hold on;

% Passband and stopband frequency markings
xline(fp,'--r','fp');
xline(fs,'--m','fs');
xline(wc/(2*pi),'--y','fc');

% Attenuation markings
yline(-Gp,'--g','G_p');
yline(-Gs,'--k','G_s');

xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Chebyshev Type-I Low Pass Filter Response');

legend('Magnitude Response','Passband Edge','Stopband Edge',...
       'Cutoff Frequency','Passband Attenuation','Stopband Attenuation');

hold off;
