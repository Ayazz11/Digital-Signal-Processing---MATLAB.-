function x = myifft(X)
% Inverse FFT using the forward FFT:  x = conj( FFT( conj(X) ) ) / N
N = length(X);
x = conj(myfft(conj(X))) / N;
end