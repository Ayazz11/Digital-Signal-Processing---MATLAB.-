clc;
clear;
while(1)
type = input("Enter signal type","s");
switch lower(type)
    case 'impulse'
        n1=input("Enter min range :\n");
        n2=input("Enter max range :\n");
        no=input("Enter Delay");
        A=input("Enter Amplitude");
        [y1,n]=gen_seq(type,n1,n2,no,A);
        subplot(3,3,1);
        stem(n,y1);
        xlabel("Amplitude");
        ylabel("time index");
        title("impulse response");
    case 'step'
        n1=input("Enter min range :\n");
        n2=input("Enter max range :\n");
        no=input("Enter Delay");
        A=input("Enter Amplitude");
        [y2,n]=gen_seq(type,n1,n2,no,A);
        subplot(3,3,2);
        stem(n,y2);
        xlabel("Amplitude");
        ylabel("time index");
        title("Step response");

    case 'ramp'
        n1=input("Enter min range :\n");
        n2=input("Enter max range :\n");
        no=input("Enter Delay");
        A=input("Enter Amplitude");
        [y3,n]=gen_seq(type,n1,n2,no,A);
        subplot(3,3,3);
        stem(n,y3);
        xlabel("Amplitude");
        ylabel("time index");
        title("ramp response");
    case 'rect'
        n1=input("Enter min range :\n");
        n2=input("Enter max range :\n");
        no=input("Enter starting index for rectangle \n");
        A=input("Enter Amplitude");
        W=input("Enter the Width");
        if no+W>n2
            error('Length of signal out of the range');
        end
        [y4,n]=gen_seq(type,n1,n2,no,W,A);
        subplot(3,3,4);
        stem(n,y4);
        xlabel("Amplitude");
        ylabel("time index");
        title("rectangular");

    case 'sin'
        n1=input("Enter min range :\n");
        n2=input("Enter max range :\n");
        no=input("Enter Delay \n");
        A=input("Enter Amplitude");
        N=input("Enter period");
        [y5,n]=gen_seq(type,n1,n2,no,A,N);  
        subplot(3,3,5);
        stem(n,y5);
        xlabel("Amplitude");
        ylabel("time index");
        title("sinusoidal signal");

    case 'exp'
        n1=input("Enter min range :\n");
        n2=input("Enter max range :\n");
        no=input("Enter Delay \n");
        A=input("Enter Amplitude");
        N=input("Enter period");
        [y6,n]=gen_seq(type,n1,n2,no,A,N);
        subplot(3,3,6);
        stem(n,real(y6));
        hold on;
        stem(n,imag(y6));
        hold off;
        xlabel("Amplitude");
        ylabel("time index");
        title("complex exponential");
    case 'sinc'
        n1=input("Enter min range :\n");
        n2=input("Enter max range :\n");
        no=input("Enter Delay \n");
        A=input("Enter Amplitude");
        N=input("Enter period");
        [y7,n]=gen_seq(type,n1,n2,no,A,N);
        subplot(3,3,7);
        stem(n,y7);
        xlabel("Amplitude");
        ylabel("time index");
        title("sinc");
     case 'tri'
        n1=input("Enter min range :\n");
        n2=input("Enter max range :\n");
        no=input("Enter Delay \n");
        A=input("Enter Amplitude");
        W=input("Enter width");
        [y8,n]=gen_seq(type,n1,n2,no,A,W);
        subplot(3,3,8);
        stem(n,y8);
        xlabel("Amplitude");
        ylabel("time index");
        title("traingular");
    case 'sign'
        N=input("Enter length");
        n1=-N;
        n2=+N;
        [y9,n] = gen_seq(type,n1,n2,0);
        subplot(3,3,9);
        stem(n,y9);
        xlabel("Amplitude");
        ylabel("time index");
        title("signum");
    otherwise
        error('gen_seq: unknown sequence type "%s"', type);
end      
end