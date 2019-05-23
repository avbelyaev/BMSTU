clear all;

img = imread('regular_39.jpg'); 
img = rgb2gray(img);
img = double (img);

sobelm_x = [ 1 2 1; 0 0 0; -1 -2 -1];
sobelm_y = [ -1 0 1; -2 0 2; -1 0 1];

% convolute with sobel kernel
w_x = conv2(img, sobelm_x); 
w_y = conv2(img, sobelm_y);
w = sqrt(double(w_x.^2 + w_y.^2));
imshow(w);

imwrite(w, gray(256), 'conv.png');
img = imread('conv.png');

binarized = imbinarize(img);
% imshow(binarized);

skeletonized = bwskel(binarized);
imshow(skeletonized);
