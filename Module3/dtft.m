function [X,w] = dtft(x,n,w)     
X = zeros(1, length(w));
for i = 1:length(x)
    X = X + x(i)*exp(-1j*w*n(i)); 
end
end
