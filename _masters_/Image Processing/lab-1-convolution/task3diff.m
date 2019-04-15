clear all; % clear workspace

img = imread('BioID_0003.pgm');
[img_h, img_w, dim] = size(img);

for i = 1:img_h
    for j = 1:img_w-1
        img(i,j) = img(i,j) - img(i, j+1);
    end
end 

for i = 1:img_h-1
    for j = 1:img_w
        img(i,j) = img(i,j) - img(i+1, j);
    end
end

imwrite(img, 'finite-difference.png');
