clear all;
clc;

clear all;
clc;
step = 500;
x = 0.01:2/step:50;
n = length(x);
y_real = zeros(1, n);
y_img = zeros(1, n);
for i = 1:n
y_real(i) = abs((-sin(x(i)) - sin(2*x(i))) / (x(i) * sqrt(2*pi)));
y_img(i) = abs((-1 + 2*cos(x(i)) - cos(2*x(i))) / (x(i) * sqrt(2*pi)));
end
grid on
plot(x, y_real);
hold on
plot(x, y_img);
legend('abs(real)', 'abs(img)');
