function [y,ny] = circular_convolution(x,nx,h,nh)
Lx=length(x);
Lh=length(h);
N=max(Lx,Lh);
x = [x, zeros(1, N - Lx)];          % pad at the end only
h = [h, zeros(1, N - Lh)];
display(x);
display(h);
n1=min(nx)+min(nh);
n2=n1+N-1;
ny=n1:n2;
y=zeros(1,N);
for i=1:1:N
    s=i-1;
    l=0;
    u=N-1;
    for j=l:1:u
        ix=j+1;
        ih=mod((s-j),N) +1;
        y(i)=y(i)+x(ix)*h(ih);
    end 
end
end
      

