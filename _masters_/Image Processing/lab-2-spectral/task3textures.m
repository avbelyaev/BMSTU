clear all;
clc;

REGULAR = 'regular_19.jpg';
NOT_REG = 'stoh_19.jpg';
img = imread(REGULAR);

WINDOW_SIZE = 50;

grayImage = rgb2gray(img);
F = fft2(double(grayImage));
S = fftshift(fftshift(F), WINDOW_SIZE);
A = abs(log2(S));

% plot(imagesc(img));

figure();
mesh(A);
% mesh(real(A));
% mesh(imag(A));
