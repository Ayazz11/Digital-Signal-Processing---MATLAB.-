clear;
% Generate a unit impulse sequence
N=input("Enter the length");
n1=(-N:1:N);
y1=zeros(1,length(n1));
a=(length(n1)+1)/2;
y1(a)=1;
subplot(3,3,1);
stem(n1,y1);
xlabel("time");
ylabel("amplitude");
title("impulse sequence");
% Generate a unit step sequence
y2 = zeros(1, length(n1));
y2(a:end) = 1;
subplot(3,3,2);
stem(n1, y2);
xlabel("time");
ylabel("amplitude");
title("step sequence");
%Generate a unit ramp signal
y3 = zeros(1, length(n1)); % Corrected the typo from 'zeroes' to 'zeros'
y3(a:end) = n1(a:end); % or y3=y2.*n1;
subplot(3,3,3);
stem(n1, y3);
xlabel("time");
ylabel("amplitude");
title("ramp signal");
%Generate a signum sequence
y5= sign(n1);
subplot(3,3,4);
stem(n1,y5);
xlabel("time");
ylabel("amplitude");
title("ramp sequence");
%Generate a sinusoidal discrete signal
y4 = 4*sin(((2*pi)/10)*n1); % Corrected 'Sin' to 'sin'
subplot(3,3,5);
stem(n1, y4);
xlabel("time");
ylabel("amplitude");
title("sinusoidal signal");
% Generate a cosine discrete signal
y7 = 4*cos(((2*pi)/10)*n1);
subplot(3,3,6);
stem(n1, y7);
xlabel("time");
ylabel("amplitude");
title("cosine signal");
% Generate a complex exponential discrete signal
%y8=y7+1j*y4;
n2= -20:20;
y8=exp(1j*((2*pi)/10)*n2);
subplot(3,3,7);
stem(n2, angle(y8));
xlabel("time");
ylabel("phase ");
title("complex exponential signal");
% Generate a triangular signal
y6 = 3.*(1-abs(n1/N)); % Create a triangular waveform
subplot(3,3,8);
stem(n1, y6);
xlabel("time");
ylabel("amplitude");
title("triangular signal");
% Generate a rectangular discrete time signal
y9 = (n1>=-N) & (n1<=N);
subplot(3,3,9);
stem(n1,y9);
xlabel("time");
ylabel("amplitude");
title("square");

