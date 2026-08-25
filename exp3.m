clear;
N1 = input("Enter the first index value: ");
N2 = input("Enter the second index value: ");
n = N1:N2;
y = input("Enter discrete sample values as a vector: ");
if length(y) ~= length(n)
    disp("The length of the input vector must match the specified index range.");
    if length(y)<length(n)
         y = [y, zeros(1, length(n) - length(y))]; % Pad with zeros if y is shorter
    end
    if length(y)>length(n)
        y = y(1:length(n)); % Truncate y if it is longer
    end
end

subplot(2,2,1)
stem(n, y);
%Delay and advance
k=input("enter delay");
delay_n=n+k;
subplot(2,2,2)
stem(delay_n,y);
advance_n=n-k;
subplot(2,2,3)
stem(advance_n,y);
