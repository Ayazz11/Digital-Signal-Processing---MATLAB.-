function [X,k] = dft(x,M)
N=length(x);
if N<M
    x=[x,zeros(1,M-N)];
end
if N > M
    for i = 1:N-M
        x(i) = x(M+i) + x(i);
    end
    x = x(1:M);
end

k=0:M-1;
X=zeros(1,M);
omega=(2*pi)/M;
W=exp(-1j*omega);

for i = 0:M-1
    for l = 0:M-1
        X(i+1)=X(i+1)+x(l+1)*(W.^(i*l));
    end
end
