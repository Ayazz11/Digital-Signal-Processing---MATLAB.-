function [X,nmul] = myfft(x)
x = x(:).';
N = length(x);
if N == 1
    X = x;
    nmul = 0;
    return
end
if mod(N,2) ~= 0
    error('N must be a power of 2');
end
[E,nE] = myfft(x(1:2:end));              % x(0), x(2), x(4), ...
[O,nO] = myfft(x(2:2:end));              % x(1), x(3), x(5), ...

W = exp(-1j*2*pi*(0:N/2-1)/N);   % defined as a row for direct computation of T
T = W .* O;                              
X = [E + T, E - T];                     
nmul = nE + nO + N/2;
end
