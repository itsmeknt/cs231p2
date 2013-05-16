function [alphaMatrix kMatrix_fg kMatrix_bg] = grabcut(im_name, rng_seed)

initGlobalVariables;
% convert the pixel values to [0,1] for each R G B channel.
im_data = double(imread([im_dir '/' im_name])) / 255;
periodIdx = find(im_name == '.');

%im_data = rgb2hsv(im_data);
% im_data = im_data(:,:,1);

% colorTransform = makecform('srgb2lab');
% im_data = applycform(im_data, colorTransform);

seg_name = im_name;
if (~isempty(periodIdx))
    seg_name = [im_name(1:periodIdx(end)) 'bmp'];
end

[im_height, im_width, channel_num] = size(im_data);
true_seg_data = double(imread([seg_dir '/' seg_name]))/255;
true_seg_data_vectorized = true_seg_data(:)';
assert(sum(true_seg_data_vectorized == 1) > 0);
figure;
% display the image
imagesc(im_data);

bbox_type = 'true';

if (strcmp(bbox_type, 'select'))
    disp('Draw a bounding box to specify the rough location of the foreground');
    
    set(gca,'Units','pixels');
    ginput(1);
    p1=get(gca,'CurrentPoint');fr=rbbox;p2=get(gca,'CurrentPoint');
    p=round([p1;p2]);
    xmin=min(p(:,1));xmax=max(p(:,1));
    ymin=min(p(:,2));ymax=max(p(:,2));
    
    xmin = max(xmin, 1);
    xmax = min(im_width, xmax);
    ymin = max(ymin, 1);
    ymax = min(im_height, ymax);
    
    bbox = [xmin ymin xmax ymax];
elseif (strcmp(bbox_type, 'true'))
    bbox = getBbox(im_name);
    xmin = bbox(1);
    ymin = bbox(2);
    xmax = bbox(3);
    ymax = bbox(4);
elseif (strcmp(bbox_type, 'all'))
    xmin = 1;
    xmax = im_width;
    ymin = 1;
    ymax = im_height;
    bbox = [xmin, ymin, xmax, ymax];
end

%{
if channel_num ~= 3
    disp('This image does not have all the RGB channels, you do not need to work on it.');
    return;
end
%}

%
% end

% use given bbox
% start
%
%
% end

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
numColors = size(im_data, 3);
im_data_vectorized = reshape(im_data, im_height*im_width, numColors)';     %im_data is [color x vectorized pixel]
bbox_vectorized = zeros(im_height, im_width);
bbox_vectorized(bbox(2):bbox(4), bbox(1):bbox(3)) = 1;
bbox_vectorized = reshape(bbox_vectorized, im_height*im_width, 1)';     %im_data is [color x vectorized pixel]

% init alpha
alpha = bg_val*ones(im_height, im_width);                       % alpha is [1 x numpixel]
alpha(ymin:ymax, xmin:xmax) = fg_val;
alpha = reshape(alpha, im_height*im_width, 1)';



% init k, GMM params to 0
k_fg = zeros(1, length(alpha));
k_bg = zeros(1, length(alpha));

pi_fg = zeros(1, numK_fg);
mu_fg = zeros(numColors,  numK_fg);
sigma_fg = zeros(numColors, numColors, numK_fg);

pi_bg = zeros(1, numK_bg);
mu_bg = zeros(numColors, numK_bg);
sigma_bg = zeros(numColors, numColors, numK_bg);

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
oldAlpha = -inf*ones(size(alpha));
while (sum((alpha-oldAlpha).^2)/length(alpha) > CONVERGENCE_DISTANCE_THRESHOLD && iter < MAX_ITER)
%    iter
    %figure;
    %imshow(reshape(alpha,im_height,im_width));
    
    % step 1
    if (~skipStep || startingStep==updateKidx)
        [k_fg] = updateK(im_data_vectorized, alpha, fg_val, bbox_vectorized, numK_fg, pi_fg, mu_fg, sigma_fg, startingStep==updateKidx && iter == 1, rng_seed);
        [k_bg] = updateK(im_data_vectorized, alpha, bg_val, bbox_vectorized, numK_bg, pi_bg, mu_bg, sigma_bg, startingStep==updateKidx && iter == 1, rng_seed);
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
        [pi_fg mu_fg sigma_fg] = updateGMM(im_data_vectorized, alpha, fg_val, k_fg, numK_fg, numColors, startingStep==updateGMMidx && iter == 1);
        [pi_bg mu_bg sigma_bg] = updateGMM(im_data_vectorized, alpha, bg_val, k_bg, numK_bg, numColors, startingStep==updateGMMidx && iter == 1);
        
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
        D_fg = computeD(im_data_vectorized, numK_fg, pi_fg, mu_fg, sigma_fg);
        D_bg = computeD(im_data_vectorized, numK_bg, pi_bg, mu_bg, sigma_bg);
        oldAlpha = alpha;
        [alpha energy] = updateAlpha(im_data, bbox_vectorized, D_fg, D_bg, alpha, beta, lambda);
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
alphaMatrix = reshape(alpha,im_height,im_width);
imshow(alphaMatrix);
kMatrix_fg = reshape(k_fg,im_height,im_width);
kMatrix_bg = reshape(k_bg,im_height,im_width);
figure;
imagesc((alphaMatrix==fg_val).*kMatrix_fg);


%{
[c_fg_energies c_fg_components] = sortComponentsByEnergy(im_data_vectorized, numK_fg, pi_fg, mu_fg, sigma_fg);
c_top_fg = chooseTopComponentsByEnergy(c_fg_energies, c_fg_components, 1/numK_fg)
c_fg_vectorized = zeros(size(k_fg));
for i = 1:length(c_top_fg)
    c_fg_vectorized = c_fg_vectorized | (k_fg==c_top_fg(i));
end

[c_bg_energies c_bg_components] = sortComponentsByEnergy(im_data_vectorized, numK_bg, pi_bg, mu_bg, sigma_bg);
c_top_bg = chooseTopComponentsByEnergy(c_bg_energies, c_bg_components, 1/numK_bg)
c_bg_vectorized = zeros(size(k_bg));


c_fg = reshape(c_fg_vectorized, im_height, im_width);
figure;
imshow((alphaMatrix==fg_val).*c_fg);
%}

%{
for i = 1:numK_fg
    c = reshape((k_fg==c_fg(i)), im_height, im_width);
    figure;
    imshow((alphaMatrix==fg_val).*c);
end
c_bg = sortComponentsByEnergy(im_data_vectorized, numK_bg, pi_bg, mu_bg, sigma_bg);
for i = 1:numK_bg
    c = reshape((k_bg==c_bg(i)), im_height, im_width);
    figure;
    imshow((alphaMatrix==bg_val).*c);
end
%}

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
%{
alpha_bg = alpha==bg_val;
alpha_fg = alpha==fg_val;
true_bg = true_seg_data_vectorized==bg_val;
true_fg = true_seg_data_vectorized==fg_val;
true_ambig = true_seg_data_vectorized > 0.3 & true_seg_data_vectorized < 0.7;
totalScore = sum(alpha_bg.*true_bg + alpha_fg.*true_fg + alpha_fg.*true_ambig);
score = totalScore/length(alpha);
%}
score = sum(alpha == true_seg_data_vectorized)/length(true_seg_data_vectorized);
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

function [energies sortedComponents] = sortComponentsByEnergy(im_data_vectorized, numK_fg, pi_fg, mu_fg, sigma_fg)
D = computeD(im_data_vectorized, numK_fg, pi_fg, mu_fg, sigma_fg);
Davg = mean(D, 2);
[energies, sortedComponents] = sort(Davg)
if (min(energies) < 0)
    energies = energies + abs(min(energies));
end
if (min(energies) < 1)
    energies = energies + (1-min(energies));
end
energies = log(energies)
end

function [components] = chooseTopComponentsByEnergy(energies, sortedComponents, threshold)
cumEnergy = cumsum(energies)/sum(energies);
endIdx = 1;
while (true)
    if (endIdx+1 > length(energies))
        break;
    end
    if (cumEnergy(endIdx+1) > threshold)
        break;
    end
    endIdx = endIdx+1;
end
components = sortedComponents(1:endIdx);
end