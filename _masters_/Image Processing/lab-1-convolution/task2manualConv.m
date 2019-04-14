clc;

img = imread('BioID_0003.pgm'); 
[img_h, img_w, dim] = size(img);

% 3x3 sobel kernel
sob_x = [ 1 2 1; 0 0 0; -1 -2 -1];
sob_y = [ -1 0 1; -2 0 2; -1 0 1];
    
% image submatrices
top_left = img(1:3, 1:3);
top_right = img(1:3, img_w-2:img_w);
bot_left = img(img_h-2:img_h, 1:3);
bot_right = img(img_h-2:img_h, img_w-2:img_w);
center = img(img_h/2-1:img_h/2+1, img_w/2-1:img_w/2+1);

fprintf('top left: %f\n', convolute(top_left,  sob_x, sob_y));
fprintf('top rght: %f\n', convolute(top_right, sob_x, sob_y));
fprintf('bot left: %f\n', convolute(bot_left,  sob_x, sob_y));
fprintf('bot rght: %f\n', convolute(bot_right, sob_x, sob_y));
fprintf('center  : %f\n', convolute(center,    sob_x, sob_y));


function res = conv_mask(m, mask)
    res = double(m(1,1))*mask(1,1) + double(m(1,2))*mask(1,2) + double(m(1,3))*mask(1,3) + double(m(2,1))*mask(2,1) + double(m(2,2))*mask(2,2) + double(m(2,3))*mask(2,3) + double(m(3,1))*mask(3,1) + double(m(3,2))*mask(3,2) + double(m(3,3))*mask(3,3);
end

function res = convolute(m, mask_x, mask_y)
    z_x = conv_mask(m, mask_x);
    z_y = conv_mask(m, mask_y);
    res = sqrt(z_x^2 + z_y^2);
end
