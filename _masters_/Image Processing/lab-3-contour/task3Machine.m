clear all;

% mark boundaries with red color
RED_COLOR = 'r';
COMPONENTS = 8;

img = imread('regular_002.jpg');
img = rgb2gray(img);
bin = imbinarize(img);

[B,L] = bwboundaries(bin, COMPONENTS, 'noholes');
imshow(label2rgb(L, @jet, [.5 .5 .5]))
hold on
for k = 1:length(B)
   boundary = B{k};
   x = boundary(:,2);
   y = boundary(:,1);
   plot(x, y, RED_COLOR, 'LineWidth', 2);
end
