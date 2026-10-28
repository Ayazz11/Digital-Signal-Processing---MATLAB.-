clc;
clear;
nx = -2:2;
nh = 0:3;
x = [1 2 0 1,1];
h = [2 2 1 1];
% for N point circular conv we need to match the size of both sequence
Lx=length(x);
Lh=length(h);
N=max(Lx,Lh); 
% padding with zeros in case we have length shortage
xc = [x, zeros(1, N - Lx)];         
hc = [h, zeros(1, N - Lh)];
%update the index array
nx_c = [nx, nx(end)+1 : nx(end)+(N-Lx)];
nh_c = [nh, nh(end)+1 : nh(end)+(N-Lh)];
% for Linear convolution using circular convolution
Ly=Lx+Lh-1; % This will be the N in case we want to implement linear conv, using circular conv.
%% Linear Convolution
[yl,n1] = linear_conv(x,nx,h,nh);
yl_ref = conv(x, h);
%% Circular Convolution
[yc,n2] = circular_convolution(xc,nx_c,hc,nh_c);
yc_ref = cconv(xc,hc,N);
%% Circular Convolution
[ylc,n] = circular_convolution([x, zeros(1,Ly-Lx)],nx,[h, zeros(1,Ly-Lh)],nh);
ylc_ref = cconv(x, h,Ly);


subplot(3,3,1);
stem(nx,x);
title("Input sequence");
xlabel("time index");
ylabel("Amplitude");

subplot(3,3,2);
stem(nh,h);
title("Filter");
xlabel("time index");
ylabel("Amplitude");

subplot(3,3,3);
stem(n1,yl);
title("Linear convolution");
xlabel("time index");
ylabel("Amplitude");

subplot(3,3,4);
stem(n2,yc);
title("Circular convolution");
xlabel("time index");
ylabel("Amplitude");

subplot(3,3,5);
stem(n2,yc_ref);
title("Ref Circular convolution");
xlabel("time index");
ylabel("Amplitude");

subplot(3,3,6);
stem(n1,yl_ref);
title("Ref Linear Convolution");
xlabel("time index");
ylabel("Amplitude");

subplot(3,3,9);
stem(n,ylc);
title("Linear Convolution using circular conv");
xlabel("time index");
ylabel("Amplitude");