function [x,n] = idtft(X,w,n)
dw = w(2) - w(1);
x = zeros(1, length(n));
for i = 1:length(n)
    f = X .* exp(1j*w*n(i));                                % integrand on the grid
    x(i) = (sum(f) - 0.5*(f(1) + f(end))) * dw / (2*pi);    % trapezoidal rule
end
end
