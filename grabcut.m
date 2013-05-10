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
bbox_vectorized = zeros(im_height, im_width);
bbox_vectorized(bbox(2):bbox(4), bbox(1):bbox(3)) = 1;
bbox_vectorized = reshape(bbox_vectorized, im_height*im_width, 1)';     %im_data is [color x vectorized pixel]

% init alpha
alpha = bg_val*ones(im_height, im_width);                       % alpha is [1 x numpixel]
alpha(ymin:ymax, xmin:xmax) = fg_val;
alpha = reshape(alpha, im_height*im_width, 1)';

% init k by k-means
[kmapFg, ~] = kMeans(im_data(:,bbox_vectorized==1)', k);                                   % kmap is [1 x numpixel]
[kmapBg, ~] = kMeans(im_data(:,bbox_vectorized==0)', k);                                   % kmap is [1 x numpixel]
kmap = zeros(1, length(alpha));
fgIdx = 1;
bgIdx = 1;
for i=1:length(alpha)
    if (alpha(i) == fg_val)
        kmap(i) = kmapFg(fgIdx);
        fgIdx = fgIdx+1;
    elseif (alpha(i) == bg_val)
        kmap(i) = kmapBg(bgIdx);
        bgIdx = bgIdx+1;
    end
end

% init GMM params to 0
pi = zeros(numAlphaValues, k);
mu = zeros(numColors, numAlphaValues, k);
sigma = zeros(numColors, numColors, numAlphaValues, k);

iter = 1;
while true
    iter
    oldPi = pi;
    oldMu = mu;
    oldSigma = sigma;
    % step 2
    [pi mu sigma] = updateGMM(alpha, kmap, im_data);
    
    % step 3
%     MAX-FLOW/MIN-CUT ENERGY MINIMIZATION
     
     if hasConverged(pi, mu, sigma, oldPi, oldMu, oldSigma)
         break;
     end
     
     % step 1
     D = computeD(alpha, im_data, pi, mu, sigma);
     [~, kmapUpdated] = min(D); 
     kmap = kmapUpdated;
     
     iter = iter+1;
end
end

function converge = hasConverged(pi, mu, sigma, oldPi, oldMu, oldSigma)
distance = euclidean_distance(pi, mu, sigma, oldPi, oldMu, oldSigma)
converge = distance < 1;
end

% step 2
function [pi mu sigma] = updateGMM(alpha, kmap, im_data)
initGlobalVariables;
pi = zeros(numAlphaValues, k);                              % pi is [2 x k]
mu = zeros(numColors, numAlphaValues, k);                   % mu is [3 x 2 x k]
sigma = zeros(numColors, numColors, numAlphaValues, k);     % sigma is [3 x 3 x 2 x k]
for i = 1:k
    kmapAtI = kmap == i;
    pixelInFg = alpha==fg_val;
    pixelInFgWithK = pixelInFg & kmapAtI;
    pixelInBg = alpha==bg_val;
    pixelInBgWithK = pixelInBg & kmapAtI;
    
    pi(fg_idx, i) = sum(pixelInFgWithK)/sum(pixelInFg);
    pi(bg_idx, i) = sum(pixelInBgWithK)/sum(pixelInBg);
    
    mu(:, fg_idx, i) = mean(im_data(:, pixelInFgWithK),2);
    mu(:, bg_idx, i) = mean(im_data(:, pixelInBgWithK),2);
    
    mu_fg_i = mu(:, fg_idx, i);
    sigma(:, :, fg_idx, i) = (mu_fg_i - mean(mu_fg_i)) * (mu_fg_i - mean(mu_fg_i))';
    
    mu_bg_i = mu(:, bg_idx, i);
    sigma(:, :, bg_idx, i) = (mu_bg_i - mean(mu_bg_i)) * (mu_bg_i - mean(mu_bg_i))';
    
    
    sigma(:, :, fg_idx, i) = eye(3);
    sigma(:, :, bg_idx, i) = eye(3);
end

redFgMu = mu(1,1,:);
greenFgMu = mu(2,1,:);
redBgMu = mu(1,2,:);
greenBgMu = mu(2,2,:);
figure(1);
scatter(redFgMu, greenFgMu, 'r');
hold on;
scatter(redBgMu, greenBgMu, 'b');
end


% step 1
function D = computeD(alpha, im_data, pi, mu, sigma)       % D is [k x numpixels]
initGlobalVariables;
numPixels = length(alpha);
D = zeros(k, numPixels);

% cache computations
term1cache = -log(pi);                                         % 2 x k
term2cache = zeros(numAlphaValues, k);                         % 2 x k
sigmaInvCache = zeros(numColors, numColors, numAlphaValues, k);      % 3 x 3 x 2 x k
for i=1:k
    term2cache(fg_idx, i) = 0.5*log(det(sigma(:,:,fg_idx, i)));
    term2cache(bg_idx, i) = 0.5*log(det(sigma(:,:,bg_idx, i)));
    
    sigmaInvCache(:, :, fg_idx, k) = pinv(sigma(:,:,fg_idx,i));
    sigmaInvCache(:, :, bg_idx, k) = pinv(sigma(:,:,bg_idx,i));
end


aa = tic;
a = tic;
for c = 1:k
    term1 = term1cache(fg_idx,c)*(alpha==fg_val) + term1cache(bg_idx,c)*(alpha==bg_val);    % 1 x pixel
    term2 = term2cache(fg_idx,c)*(alpha==fg_val) + term2cache(bg_idx,c)*(alpha==bg_val);    % 1 x pixel
    
    muForAllPixel = bsxfun(@times, mu(:,fg_idx,c), alpha==fg_val) + bsxfun(@times, mu(:,bg_idx,c), alpha==bg_val);    % 3 x pixel
    firstMoment = im_data - muForAllPixel;                                          % 3 x pixel
    
    term3fgHalf = 0.5*firstMoment'*sigmaInvCache(:,:,fg_idx,c);                          % pixel x 3
    term3fg = zeros(1, numPixels);
    for i=1:numPixels
        term3fg(i) = term3fgHalf(i,:) * firstMoment(:,i);
    end
    term3bgHalf = 0.5*firstMoment'*sigmaInvCache(:,:,bg_idx,c);                          % pixel x 3
    term3bg = zeros(1, numPixels);
    for i=1:numPixels
        term3bg(i) = term3bgHalf(i,:) * firstMoment(:,i);
    end
    term3 =  term3fg.*(alpha==fg_val) + term3bg.*(alpha==bg_val);                               % 1 x pixel
    D(c, :) = term1+term2+term3;
end
bb = toc(aa)

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