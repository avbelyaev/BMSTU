clear all;

img = imread('BioID_0003.pgm'); 
[img_h, img_w, dim] = size(img);
% imshow(img);

% 3x3 kernels
sobelm_x = [ 1 2 1; 0 0 0; -1 -2 -1];
sobelm_y = [ -1 0 1; -2 0 2; -1 0 1];

prewit_x = [ -1 1 1; -1 -2 1; -1 1 1];
prewit_y = [ 1 1 1; -1 -2 1; -1 -1 1];

kirsch_x = [ -3 -3 5; -3 0 5; -3 -3 5];
kirsch_y = [ -3 5 5; -3 0 5; -3 -3 -3];

robin3_x = [ -1 0 1; -1 0 1; -1 0 1];
robin3_y = [ 0 1 1; -1 0 1; -1 -1 1];

robin5_x = [ -1 0 1; -2 0 2; -1 0 1];
robin5_y = [ 0 1 2; -1 0 1; -2 1 0];
    

w_x = conv2(img, robin5_x); 
w_y = conv2(img, robin5_y);
W = sqrt(double(w_x.^2 + w_y.^2)); % wiki formula

% draw image convolution
imwrite(W, gray(256), 'robin5.png');
% imshow(W, []);

% draw convolution plot
f = figure();
plot(W);
saveas(gcf, 'robin5-p.png');