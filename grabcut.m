function grabcut(im_name)

initGlobalVariables;

% convert the pixel values to [0,1] for each R G B channel.
im_data = double(imread(im_name)) / 255;
% display the image
imagesc(im_data);

% a bounding box initialization
disp('Draw a bounding box to specify the rough location of the foreground');
set(gca,'Units','pixels');
ginput(1);
p1=get(gca,'CurrentPoint');fr=rbbox;p2=get(gca,'CurrentPoint');
p=round([p1;p2]);
xmin=min(p(:,1));xmax=max(p(:,1));
ymin=min(p(:,2));ymax=max(p(:,2));
[im_height, im_width, channel_num] = size(im_data);
xmin = max(xmin, 1);
xmax = min(im_width, xmax);
ymin = max(ymin, 1);
ymax = min(im_height, ymax);

bbox = [xmin ymin xmax ymax];
line(bbox([1 3 3 1 1]),bbox([2 2 4 4 2]),'Color',[1 0 0],'LineWidth',1);
if channel_num ~= 3
    disp('This image does not have all the RGB channels, you do not need to work on it.');
    return;
end

% grabcut algorithm
disp('grabcut algorithm');

% reshape image
im_data = reshape(im_data, im_height*im_width, numColors)';     %im_data is [color x vectorized pixel]

% init alpha
alpha = bg_val*ones(im_height, im_width);                       % alpha is [1 x numpixel]
alpha(ymin:ymax, xmin:xmax) = fg_val;
alpha = reshape(alpha, im_height*im_width, 1)';

% init k
kmap = ones(1, length(alpha));                                  % kmap is [1 x numpixel]

% RANDOM INITIALIZATION
kmap = randi([1, k], 1, length(alpha));

while true
    
    % step 2
    [pi mu sigma] = updateGMM(alpha, kmap, im_data);
    
    % step 3
%     MAX-FLOW/MIN-CUT ENERGY MINIMIZATION
     
     if hasConverged()
         break;
     end
     
     % step 1
     
end
end

function converge = hasConverged()
converge = false;
end

function [pi mu sigma] = updateGMM(alpha, kmap, im_data)
initGlobalVariables;
pi = zeros(numAlphaValues, k);
mu = zeros(numColors, numAlphaValues, k);
sigma = zeros(numColors, numColors, numAlphaValues, k);
for i = 1:k
    pixelWithK = kmap==i;
    fgPixelWithK = pixelWithK & alpha==fg_val;
    bgPixelWithK = pixelWithK & alpha==bg_val;
    
    pi(fg_idx, i) = sum(fgPixelWithK)/sum(pixelWithK);
    pi(bg_idx, i) = sum(bgPixelWithK)/sum(pixelWithK);
    
    mu(:, fg_idx, i) = mean(im_data(:, fgPixelWithK),2);
    mu(:, bg_idx, i) = mean(im_data(:, bgPixelWithK),2);
    
    mu_fg_i = mu(:, fg_idx, i);
    sigma(:, :, fg_idx, i) = (mu_fg_i - mean(mu_fg_i)) * (mu_fg_i - mean(mu_fg_i))';
    
    mu_bg_i = mu(:, bg_idx, i);
    sigma(:, :, bg_idx, i) = (mu_bg_i - mean(mu_bg_i)) * (mu_bg_i - mean(mu_bg_i))';
end
end


