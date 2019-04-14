clear all; % clear workspace

x = -30:0.1:30; 
N = length(x);

for i = 1:N
    f(i) = func(x(i));
    g(i) = sin(x(i));
end

w = conv(f, g);
fprintf('max: %s', max(w));

% plot F, G
figure(1)
plot(x, g);
hold on % draw multiple plots in same window
plot(x, f);
legend('g(x)', 'f(x)');
hold off

% plot convolution
figure(2)
plot(x, f);


function y = func(x)
    if x >= -1 && x <= 1
        y = exp(-1 * x^2);
    else
        y = 0;
    end
end