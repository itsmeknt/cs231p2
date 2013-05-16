data_dir = 'data';
im_dir = [data_dir '/' 'data_GT'];
seg_dir = [data_dir '/' 'seg_GT'];
bbox_dir = [data_dir '/' 'ICCV09_new_bounding_boxes'];

allowCache = true;
kCacheFile = 'cache/k';
betaCacheFile = 'cache/beta';

updateKidx = 1;
updateGMMidx = 2;
updateAlphaIdx = 3;

% init variables
numK_fg = 5;
numK_bg = 5;
numAlphaValues = 2;
fg_idx = 1;
bg_idx = 2;
fg_val = 1;
bg_val = 0;
r_idx = 1;
g_idx = 2;
b_idx = 3;

lambda = 50; % author stated 50 was best

MAX_ITER = 30;
CONVERGENCE_DISTANCE_THRESHOLD = 0.001;

random_pixel_image_max_size = 30000;
% type 1 - neighboring 4
% type 2 - neighboring 8
Vtype = '8pixel';
if strcmp(Vtype, 'none') == 0
    Voffsets = [];
elseif strcmp(Vtype, '8pixel')
    Voffsets = [1, 0; 1, 1; 0, 1; -1, 1; -1, 0; -1, -1; 0, -1; 1, -1];
elseif strcmp(Vtype, '16pixel')
    Voffsets = [1, 0; 1, 1; 0, 1; -1, 1; -1, 0; -1, -1; 0, -1; 1, -1; 2,0;2,2;0,2;-2,2;-2,0;-2,-2;0,-2;2,-2];
end

epsilon = 10^-100;
delta = 10^-6;

%{
ng kmeans
Improper assignment with rectangular empty matrix.

Error in kMeans (line 21)
        centers(i,:)=minData+interval/2;

Error in updateK>initK (line 62)
[kBg, ~] = kMeans(im_data_vectorized(:,bbox_vectorized==0)', K);
% k is [1 x numpixel]

Error in updateK (line 4)
    k = initK(im_data_vectorized, alpha, bbox_vectorized);

Error in grabcut (line 88)
        k = updateK(im_data_vectorized, alpha, bbox_vectorized, pi, mu, sigma,
        startingStep==updateKidx && iter == 1);
        %}