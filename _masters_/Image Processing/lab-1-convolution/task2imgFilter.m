import Conv.*

clear all;

img = imread('BioID_0003.pgm'); 
[img_h, img_w, dim] = size(img);
% imshow(img);


% 3x3 kernels
sobelm_x = [ 1 2 1; 0 0 0; -1 -2 -1];
sobelm_y = [ -1 0 1; -2 0 2; -1 0 1];
c1 = Conv('sobel', sobelm_x, sobelm_y);
convolute(c1, img);

prewit_x = [ -1 1 1; -1 -2 1; -1 1 1];
prewit_y = [ 1 1 1; -1 -2 1; -1 -1 1];
c2 = Conv('prew', prewit_x, prewit_y);
convolute(c2, img);

kirsch_x = [ -3 -3 5; -3 0 5; -3 -3 5];
kirsch_y = [ -3 5 5; -3 0 5; -3 -3 -3];
c3 = Conv('kir', kirsch_x, kirsch_y);
convolute(c3, img);

robin3_x = [ -1 0 1; -1 0 1; -1 0 1];
robin3_y = [ 0 1 1; -1 0 1; -1 -1 1];
c4 = Conv('rob3', robin3_x, robin3_y);
convolute(c4, img);

robin5_x = [ -1 0 1; -2 0 2; -1 0 1];
robin5_y = [ 0 1 2; -1 0 1; -2 1 0];
c5 = Conv('rob5', robin5_x, robin5_y);
convolute(c5, img);


