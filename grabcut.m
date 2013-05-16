function alpha = grabcut(im_name)

initGlobalVariables;
% convert the pixel values to [0,1] for each R G B channel.
im_data = double(imread([im_dir '/' im_name])) / 255;
periodIdx = find(im_name == '.');
seg_name = im_name;
if (~isempty(periodIdx))
    seg_name = [im_name(1:periodIdx(end)) 'bmp'];
end

true_seg_data = double(imread([seg_dir '/' seg_name]))/255;
true_seg_data_vectorized = true_seg_data(:)';
assert(sum(true_seg_data_vectorized == 1) > 0);
% display the image
imagesc(im_data);

%{
% a bounding box initialization
disp('Draw a bounding box to specify the rough location of the foreground');

%
set(gca,'Units','pixels');
ginput(1);
p1=get(gca,'CurrentPoint');fr=rbbox;p2=get(gca,'CurrentPoint');
p=round([p1;p2]);
xmin=min(p(:,1));xmax=max(p(:,1));
ymin=min(p(:,2));ymax=max(p(:,2));
%
 
%{
 xmin = 10;
 xmax = 412;
 ymin = 10;
 ymax = 412;
%}

xmin = max(xmin, 1);
xmax = min(im_width, xmax);
ymin = max(ymin, 1);
ymax = min(im_height, ymax);

bbox = [xmin ymin xmax ymax];

if channel_num ~= 3
    disp('This image does not have all the RGB channels, you do not need to work on it.');
    return;
end
%}

[im_height, im_width, channel_num] = size(im_data);
bbox = getBbox(im_name);
xmin = bbox(1);
ymin = bbox(2);
xmax = bbox(3);
ymax = bbox(4);
line(bbox([1 3 3 1 1]),bbox([2 2 4 4 2]),'Color',[1 0 0],'LineWidth',1);


% grabcut algorithm
% disp('grabcut algorithm');

fileSlashIdx = find(im_name=='/');
filePeriodIdx = find(im_name=='.');
im_file_name = im_name;
if ~isempty(filePeriodIdx)
    im_file_name = im_file_name(1:filePeriodIdx(end)-1);
end
if ~isempty(fileSlashIdx)
    im_file_name = im_file_name(fileSlashIdx(end)+1:end);
end

% reshape image
im_data_vectorized = reshape(im_data, im_height*im_width, numColors)';     %im_data is [color x vectorized pixel]
bbox_vectorized = zeros(im_height, im_width);
bbox_vectorized(bbox(2):bbox(4), bbox(1):bbox(3)) = 1;
bbox_vectorized = reshape(bbox_vectorized, im_height*im_width, 1)';     %im_data is [color x vectorized pixel]

% init alpha
alpha = bg_val*ones(im_height, im_width);                       % alpha is [1 x numpixel]
alpha(ymin:ymax, xmin:xmax) = fg_val;
alpha = reshape(alpha, im_height*im_width, 1)';



% init k, GMM params to 0
k = zeros(1, length(alpha));
pi = zeros(numAlphaValues, K);
mu = zeros(numColors, numAlphaValues, K);
sigma = zeros(numColors, numColors, numAlphaValues, K);
oldPi = pi;
oldMu = mu;
oldSigma = sigma;

% compute beta
randomPixels = im_data_vectorized;
if (size(im_data_vectorized,2)  > random_pixel_image_max_size)
    randomPixels = randomPixelSample(im_data_vectorized, random_pixel_image_max_size);
end
beta = computeBeta(randomPixels, im_file_name);
%beta=0.01;

% starting step - decides which step we start the initialization process
startingStep = 1;

converge = false;
iter = 1;
skipStep = true;
scores = [];
energies = [];
while (iter<MAX_ITER)
%    iter
    figure;
    imshow(reshape(alpha,im_height,im_width));
    
    % step 1
    if (~skipStep || startingStep==updateKidx)
        k = updateK(im_data_vectorized, alpha, bbox_vectorized, pi, mu, sigma, startingStep==updateKidx && iter == 1);
        if (startingStep == updateKidx)
            skipStep = false;
        end
%         fg=alpha==fg_val;
%         bg=alpha==bg_val;
%         figure;
%         kfg=fg.*k;
%         kbg=bg.*k;
%         imshow(reshape(kfg/5,im_height,im_width));
%         figure;
%         imshow(reshape(kbg/5,im_height,im_width));
    end
    
    % step 2
    if (~skipStep || startingStep==updateGMMidx)
        oldPi = pi;
        oldMu = mu;
        oldSigma = sigma;
        [pi mu sigma] = updateGMM(im_data_vectorized, alpha, k, startingStep==updateGMMidx && iter == 1);
        
        %{
        if (iter > 1 && hasConverged(pi, mu, sigma, oldPi, oldMu, oldSigma))
            converge = true;
        else
            converge = false;
        end
        %} 
        
        if (startingStep == updateGMMidx)
            skipStep = false;
        end
        
        % plot GMM params
        %{
        figure; hold on;
        redFgMu = mu(1,1,:);
        greenFgMu = mu(2,1,:);
        redBgMu = mu(1,2,:);
        greenBgMu = mu(2,2,:);
        scatter(redFgMu, greenFgMu, 'r');
        scatter(redBgMu, greenBgMu, 'b');
        %}
    end
    
    % step 3
    if (~skipStep || startingStep==updateAlphaIdx)
        [alpha energy] = updateAlpha(im_data, im_data_vectorized, k, alpha, pi, mu, sigma, lambda, beta, startingStep==updateAlphaIdx && iter == 1, iter, bbox_vectorized);
        energies = [energies, energy];
        score = computeScore(alpha, true_seg_data_vectorized);
        scores = [scores, score];
        
        % energy
        % score
        
        if (startingStep == updateAlphaIdx)
            skipStep = false;
        end
    end
    
%     if (converge)
%         break;
%     end
    
    skipStep = false;
    iter = iter+1;
end

figure;
imshow(reshape(alpha,im_height,im_width));

figure;
plot(1:iter-1, energies);
title('energy vs iteration');

figure;
plot(1:iter-1, scores);
title('scores vs iteration');
final_score = scores(end)

end

function score = computeScore(alpha, true_seg_data_vectorized)
initGlobalVariables;
alpha_bg = alpha==bg_val;
alpha_fg = alpha==fg_val;
true_bg = true_seg_data_vectorized==0;
true_fg = true_seg_data_vectorized.*(true_seg_data_vectorized>0);
totalScore = sum(alpha_bg.*true_bg + alpha_fg.*true_fg);
score = totalScore/length(alpha);
end

function randomPixels = randomPixelSample(im_data_vectorized, n)
idx = randsample(size(im_data_vectorized,2), n);
randomPixels = im_data_vectorized(:, idx);
end

function beta = computeBeta(im_data_vectorized, im_file_name)
initGlobalVariables;
filename = [betaCacheFile '-' im_file_name];
% disp('computing beta');
if (allowCache)
    try
        load([filename]);
        return;
    catch
    end
end

numPixels = size(im_data_vectorized, 2);
% rep_im_data = repmat(im_data_vectorized, 1, 2);
% term1 = im_data_vectorized;

sumSq = 0;
numPixelDiff = 0;
for i=1:numPixels-1
    term1 = im_data_vectorized(:,1+i:end);
    term2 = im_data_vectorized(:,1:end-i);
    
    diff = term1-term2;
    numPixelDiff = numPixelDiff + size(diff, 2);
    sumSq = sumSq + sum(sum(diff.*diff));
end
%{
for i = 1:numPixels
    term2 = rep_im_data(:, i:i+numPixels-1);
    
    truncateLength = i-1;
    term2 = [zeros(3, truncateLength), term2(:, truncateLength+1:end)];
    
    assert(size(term1,1) == size(term2,1));
    assert(size(term1,2) == size(term2,2));
    
    diff = term1(:, truncateLength+1:end)-term2(:, truncateLength+1:end);
    numPixelDiff = numPixelDiff + size(diff,2);
    sumSq = sumSq + sum(sum(diff.*diff));
end
%}
sumSq = sumSq/numPixelDiff;

beta = 1/(2*sumSq);
% disp('computing beta done!');
save(filename, 'beta');
end

function converge = hasConverged(pi, mu, sigma, oldPi, oldMu, oldSigma)
initGlobalVariables;
distance = euclidean_distance(pi, mu, sigma, oldPi, oldMu, oldSigma)
converge = distance < CONVERGENCE_DISTANCE_THRESHOLD;
end


function distance = euclidean_distance(pi, mu, sigma, oldPi, oldMu, oldSigma)
sumSq= 0 ;

piDiff = pi-oldPi;
piDiffSq = piDiff.*piDiff;
sumSq = sumSq + sum(piDiffSq(:));

muDiff = mu - oldMu;
muDiffSq = muDiff.*muDiff;
sumSq = sumSq + sum(muDiffSq(:));

sigmaDiff = sigma - oldSigma;
sigmaDiffSq = sigmaDiff.*sigmaDiff;
sumSq = sumSq + sum(sigmaDiffSq(:));

distance = sqrt(sumSq);
end

function [bbox] = getBbox(im_file_name)
initGlobalVariables;
periodIdx = find(im_file_name=='.');
if (~isempty(periodIdx))
    im_file_name = im_file_name(1:periodIdx(end)-1);
end
im_file_name = [im_file_name '.txt'];
fid = fopen([bbox_dir '/' im_file_name]);
bbox = fscanf(fid, '%g', [1 inf]);
fclose(fid);
end