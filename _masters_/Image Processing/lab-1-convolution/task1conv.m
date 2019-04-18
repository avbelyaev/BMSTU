clear all; % clear workspace

x = -30:0.1:30; 
N = length(x);

for i = 1:N
    f(i) = func(x(i));
    g(i) = 1 / (pi*(1 + x(i)^2));
end

w = conv(f, g, 'same');
fprintf('max: %f', max(w));

% draw F, G
figure(1)
plot(x, g);
xlim([-10 10]);
hold on     % draw multiple plots in same window
plot(x, f);
legend('g(x)', 'f(x)');
hold off
saveas(gcf, 'fg.png');

% draw convolution
figure(2)
plot(x, w);
xlim([-10 10]);
saveas(gcf, 'fg-conv.png');


function y = func(x)
    if x >= -1 && x <= 1
        y = exp(-1 * x^2);
    else
        y = 0;
    end
end
