allowCache = true;
kCacheFile = 'cache/k';
betaCacheFile = 'cache/beta';

updateKidx = 1;
updateGMMidx = 2;
updateAlphaIdx = 3;

% init variables
K = 5;
numAlphaValues = 2;
numColors = 3;
fg_idx = 1;
bg_idx = 2;
fg_val = 1;
bg_val = 0;
r_idx = 1;
g_idx = 2;
b_idx = 3;

lambda = 50; % author stated 50 was best

random_pixel_image_max_size = 10000;
% type 1 - neighboring 4
% type 2 - neighboring 8
Vtype = 1;
Voffsets = [];
Vdim = 0;
if Vtype == 1
    Vdim = 4;
    Voffsets = [1, 0; 0, 1; -1, 0; 0, -1];
end

epsilon = 10^-100;