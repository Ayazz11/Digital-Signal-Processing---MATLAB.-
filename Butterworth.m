clc;
clear;
fp=input("Enter pass band cutt off");
fs=input("Enter stop band cutt off");
wp=2*pi*fp;
ws=2*pi*fs;
Gp=input("Enter pass band attenuation in db");
Gs=input("Enter stop band attenuation in db");
[n,wc]=buttord(wp,ws,Gp,Gs,'s'); %% does not design theb filter , only determines the order and cutt off freq of the filter 
[b,a]=butter(n,wc,'s'); %% Gives actual butterworth filter , s denotes that it is analog filter, n=Numerator polynomial coeff, a= deno...
omega=linspace(0,12000*pi,2000);
h=freqs(b,a,omega); %% freq response from the transfer function coefficients
gain=20*log10(abs(h));
figure;
plot(omega/(2*pi),gain,'b','LineWidth',1.5);
grid on;
hold on;

% Passband and stopband frequency markings
xline(fp,'--r','fp');
xline(fs,'--m','fs');
xline(wc/(2*pi),'--y','fc');

% Passband and stopband attenuation markings
yline(-Gp,'--g','G_p');
yline(-Gs,'--k','G_s');

xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Butterworth Low Pass Filter Response');

legend('Magnitude Response','Passband Edge','Stopband Edge',...
       'Passband Attenuation','Stopband Attenuation');

hold off;
