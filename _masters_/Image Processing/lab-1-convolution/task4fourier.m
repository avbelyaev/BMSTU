clear all;

N = 100;    % steps

x = -pi;    % from
step = 0.05;
x_max = pi; % to

i = 1;
while x < x_max
    ys(i) = fourier(x);
    xs(i) = x;
    
    x = x + step;
    i = i + 1;
end

% draw convolution plot
f = figure();
plot(xs, ys);
filename = ['fr-' num2str(N) '.png'];
saveas(gcf, filename);


function sum = fourier(x)
    sum = 0;
    for n = 1:N
        a0 = (pi*x - x^2 / 2) / (2*pi);
        an = (2*sin(pi*n / 2)) / n;
        bn = (2*pi*n*cos(pi*n / 2) - 4*sin(pi*n / 2)) / (pi*n^2);
        sum = sum + a0 + (an*cos(n*x) + bn*sin(n*x));
    end
end
