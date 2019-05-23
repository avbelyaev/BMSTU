clear all;

BACKGROUND = 0;
CONTOUR = 255;
% own value for each image
THRESHOLD = 190;

img = imread('regular_33.jpg');
img = rgb2gray(img);
[img_h, img_w, dim] = size(img);

for i = 1:img_h
    for j = 1:img_w-1
        if img(i, j) > THRESHOLD
            img(i, j) = CONTOUR;
        else
            img(i, j) = BACKGROUND; 
        end
    end
end 
imshow(img);
