function [label, n] = count_dots(z_matrix, dim, sigma)

% laplacian-of-gaussian filter
dat = (z_matrix);
w = fspecial('log', dim, sigma);
filtered = imfilter(dat, w, 'replicate');
filtered = filtered * (-1);
thresh = graythresh(filtered); % otsu thresholding 
bin = thresh < filtered;
bin_2 = imfill(bin, 'holes');
bin_3 = bwareafilt(bin_2, [5 200]); % filer dots based on area
[label, old_n] = bwlabel(bin_3);
n = old_n - 1; 