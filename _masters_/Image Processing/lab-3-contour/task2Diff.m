clear all;

img = imread('regular_002.jpg');
img = rgb2gray(img);

D = diff(img);
% imshow(D);

binarized = imbinarize(D);
imshow(binarized);
