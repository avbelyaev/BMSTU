clear all;
clc;

% NUM_OF_PERIODS = 100;
NUM_OF_PERIODS = 200;

PERIOD = 2*pi;
MIN = 0;
MAX = NUM_OF_PERIODS * PERIOD;
% step = 2 * PERIOD;
% step = 1 * PERIOD;
step = 0.5 * PERIOD;

x = MIN:step:MAX;
% y = sin(x);
y = 1 + cos(pi + x) + cos(x - pi);
y1 = fft(y);

grid on;
hold on;
% plot(x, y);     
plot(x, abs(imag(y1)));     
legend(sprintf('Step: %.2f\nNumOfPeriods: %d',step, NUM_OF_PERIODS));
