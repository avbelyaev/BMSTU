clear all;

N = 100;    % steps

x = -pi;    % from
step = 0.01;
x_max = pi; % to

i = 1;
while x < x_max
    ys(i) = fourier(x, N);
    xs(i) = x;
    
    x = x + step;
    i = i + 1;
end

% draw convolution plot
f = figure();
plot(xs, ys);
filename = ['fr-' num2str(N) '.png'];
saveas(gcf, filename);


function sum = fourier(x, N)
    sum = (pi*x - x^2 / 2) / (2*pi);  % a0
    for n = 1:N
        an = (2*sin(pi*n / 2)) / n;
        bn = (2*pi*n*cos(pi*n / 2) - 4*sin(pi*n / 2)) / (pi*n^2);
        sum = sum + (an*cos(n*x) + bn*sin(n*x));
    end
end
