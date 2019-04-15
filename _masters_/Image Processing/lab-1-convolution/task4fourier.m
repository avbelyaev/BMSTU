clear all;
clc;

pi=3.14;
n = 5;
sum=0;
y=sinh(x); % function we want 
for n=1:n
    a0=(1/pi)*int(y,x,-pi,pi);
    an=(1/pi)*int(y*cos(n*x),x,-pi,pi); 
    bn=(1/pi)*int(y*sin(n*x),x,-pi,pi); 
    sum=sum+(an*cos(n*x)+bn*sin(n*x));
end

ezplot(x,y,[-pi,pi]);
grid on;
hold on; 
ezplot(x,(sum+a0/2),[-pi,pi]);
