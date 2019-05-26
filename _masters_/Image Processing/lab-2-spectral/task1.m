clear all;
clc;

step = 200;
xs = -1:2/step:1;
n = length(xs);

y_real = zeros(1, n);
y_img = zeros(1, n);
y = zeros(1, n);

for i = 1:n
    x = xs(i);
    y_real(i) = (-1 + cos(x) + 2*x*sin(x)) / (sqrt(2*pi)*x*x);
    y_img(i) = (x - sin(x)) / (sqrt(2*pi)*x*x);
    y(i) = func(x);
end

grid on
hold on
% plot(xs, y);
plot(xs, y_real);
plot(xs, y_img);
legend('abs(real)', 'abs(img)');


function y = func(x)
    if x >= -1 && x < 0
        y = 1;
    elseif x == 0
        y = 1/2;
    elseif x > 0 && x <= 1
        y = x;
    end
end
