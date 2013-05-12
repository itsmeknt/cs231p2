function [k] = updateK(im_data, alpha, bbox_vectorized, pi, mu, sigma, init )
disp('step 1: finding k');
if init
    k = initK(im_data, alpha, bbox_vectorized);
else
    D = computeD(im_data, alpha, pi, mu, sigma);
    [~, kUpdated] = min(D);
    k = kUpdated;
end
disp('step 1 done!');
end



function D = computeD(im_data, alpha, pi, mu, sigma)       % D is [k x numpixels]
initGlobalVariables;
numPixels = length(alpha);
D = zeros(K, numPixels);

% cache computations
term1cache = -log(pi);                                         % 2 x K
term2cache = zeros(numAlphaValues, K);                         % 2 x K
sigmaInvCache = zeros(numColors, numColors, numAlphaValues, K);      % 3 x 3 x 2 x K
for i=1:K
    term2cache(fg_idx, i) = 0.5*log(det(sigma(:,:,fg_idx, i)));
    term2cache(bg_idx, i) = 0.5*log(det(sigma(:,:,bg_idx, i)));
    
    sigmaInvCache(:, :, fg_idx, i) = pinv(sigma(:,:,fg_idx,i));
    sigmaInvCache(:, :, bg_idx, i) = pinv(sigma(:,:,bg_idx,i));
end


aa = tic;
a = tic;
for c = 1:K
    term1 = term1cache(fg_idx,c)*(alpha==fg_val) + term1cache(bg_idx,c)*(alpha==bg_val);    % 1 x pixel
    term2 = term2cache(fg_idx,c)*(alpha==fg_val) + term2cache(bg_idx,c)*(alpha==bg_val);    % 1 x pixel
    
    muForAllPixel = bsxfun(@times, mu(:,fg_idx,c), alpha==fg_val) + bsxfun(@times, mu(:,bg_idx,c), alpha==bg_val);    % 3 x pixel
    firstMoment = im_data - muForAllPixel;                                          % 3 x pixel
    
    term3fgHalf = firstMoment'*sigmaInvCache(:,:,fg_idx,c);                          % pixel x 3
    term3fg = 0.5*sum(term3fgHalf' .* firstMoment);
    term3bgHalf = firstMoment'*sigmaInvCache(:,:,bg_idx,c);                          % pixel x 3
    term3bg = 0.5*sum(term3bgHalf' .* firstMoment);
    term3 =  term3fg.*(alpha==fg_val) + term3bg.*(alpha==bg_val);                               % 1 x pixel
    D(c, :) = term1+term2+term3;
end
bb = toc(aa)
end



function k = initK(im_data, alpha, bbox_vectorized)       % D is [k x numpixels]
initGlobalVariables;
k = zeros(1, size(im_data,2));
try
    load([kCacheFile]);
catch
    disp('initializing k using kmeans');
    [kFg, ~] = kMeans(im_data(:,bbox_vectorized==1)', K);                                   % k is [1 x numpixel]
    [kBg, ~] = kMeans(im_data(:,bbox_vectorized==0)', K);                                   % k is [1 x numpixel]
    k = zeros(1, length(alpha));
    fgIdx = 1;
    bgIdx = 1;
    for i=1:length(alpha)
        if (alpha(i) == fg_val)
            k(i) = kFg(fgIdx);
            fgIdx = fgIdx+1;
        elseif (alpha(i) == bg_val)
            k(i) = kBg(bgIdx);
            bgIdx = bgIdx+1;
        end
    end
    save(kCacheFile, 'k');
end
end