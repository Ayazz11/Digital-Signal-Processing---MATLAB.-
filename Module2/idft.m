function [x,n]=idft(X)
M=length(X);
n=0:M-1;
x=zeros(1,M);
omega=(2*pi)/M;
W=exp(j*omega);
for i = 0:M-1
    for k = 0:M-1
        x(i+1)=x(i+1)+X(k+1)*W.^(k*i);
    end
    x(i+1) = x(i+1)/M;
end

