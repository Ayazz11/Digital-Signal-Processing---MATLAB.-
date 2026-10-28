    function [y,ny] = linear_conv(x,nx,h,nh)
    %% define upper and lower index for linear conv
    Lx=length(x);
    Lh=length(h);
    Ly = Lx+Lh-1;
    n1 = (min(nx) + min(nh));
    n2 = (max(nx) + max(nh));
    ny = n1:n2;
    y=zeros(1,Ly);
    for i=1:1:Ly
       s=i-1;
       l=max(0, s-Lh+1);
       u=min(Lx-1,s);
       for j=l:1:u
       ix = j + 1; %ix is for matlab index 
       ih = (s-j) + 1; %iy is for matlab index 
       y(i) = y(i) + x(ix)*h(ih);
       end
    end
    end
