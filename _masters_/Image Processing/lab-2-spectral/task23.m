clear all;
clc;

REGULAR = 'regular.jpg';
NOT_REG = 'stoh.jpg';
img = imread(NOT_REG);

WINDOW_SIZE = 50;

grayImage = rgb2gray(img);
F = fft2(double(grayImage));
S = fftshift(fftshift(F), WINDOW_SIZE);
A = abs(log2(S));

subplot(2,2, 1), imagesc(img);
subplot(2,2, 2), imagesc(A);
subplot(2,2, 3), imagesc(real(A));
subplot(2,2, 4), imagesc(imag(A));
