clear;
N=input("Enter the length of input sequence");
n=0:N-1;
%define the LTI system Transfer function
b=[0.8 -0.44 0.16 0.02]; %coefficients in numerator
a=[1 0.7 -0.45 -0.6]; %coefficients in denominator
x=[1, zeros(1,N-1)]; %input is del(n), for impulse response
h=filter(b,a,x);
N1 = input("Enter the first index value: ");
N2 = input("Enter the second index value: ");
n1 = N1:N2;
y1=input("Enter the sequence");
if length(y1)<length(n1)
    disp("y(n) length is smaller then specified index, hence zero padded");
    y1=[y1, zeros(length(n1)-length(y1))];
elseif length(y1)>length(n1)
    y1=y1(1:length(n1));
else
    disp("Lengths match. No modification needed.");
end
n_o=min(n)+min(n1):max(n)+max(n1);
y_o=conv(y1,h);    
% Plot the impulse response of the LTI system
subplot(2,2,1);
stem(n, h);
xlabel('n');
ylabel('Impulse Response');
title('Impulse Response of the LTI System');
% Plot the input sequence to the LTI system
subplot(2,2,2);
stem(n1, y1);
xlabel('n1');
ylabel('input sequence');
title('Input to LTI system');
% Plot the output sequence of the LTI system
subplot(2,2,3);
stem(n_o, y_o);
xlabel('n_o');
ylabel('output sequence');
title('output of LTI system');
%deconvolution
y_deconv=deconv(y_o,h);
subplot(2,2,4);
stem(n1, y_deconv);
xlabel('n_o');
ylabel('Deconvolution output');
