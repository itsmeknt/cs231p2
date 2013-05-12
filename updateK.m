function [k] = updateK(im_data_vectorized, alpha, bbox_vectorized, pi, mu, sigma, init )
disp('step 1: finding k');
if init
    k = initK(im_data_vectorized, alpha, bbox_vectorized);
else
    D = computeD(im_data_vectorized, alpha, pi, mu, sigma);
    [~, kUpdated] = min(D);
    k = kUpdated;
end
disp('step 1 done!');
end



function D = computeD(im_data_vectorized, alpha, pi, mu, sigma)       % D is [k x numpixels]
initGlobalVariables;
numPixels = length(alpha);
D = zeros(K, numPixels);

% cache computations
term1cache = -log(pi);                                         % 2 x K

start = tic;
alpha_fg_idx = alpha==fg_val;
alpha_bg_idx = alpha==bg_val;
for c = 1:K
    term1 = term1cache(fg_idx,c)*alpha_fg_idx + term1cache(bg_idx,c)*alpha_bg_idx;    % 1 x pixel
    
    term2_fg = 0.5*log(det(sigma(:,:,fg_idx, c)));
    term2_bg = 0.5*log(det(sigma(:,:,bg_idx, c)));
    term2 = term2_fg*alpha_fg_idx + term2_bg*alpha_bg_idx;    % 1 x pixel
    
    muForAllPixel = bsxfun(@times, mu(:,fg_idx,c), alpha_fg_idx) + bsxfun(@times, mu(:,bg_idx,c), alpha_bg_idx);    % 3 x pixel
    firstMoment = im_data_vectorized - muForAllPixel;                                          % 3 x pixel
    
    sigmaInv_fg = pinv(sigma(:,:,fg_idx,c));
    sigmaInv_bg = pinv(sigma(:,:,bg_idx,c));
    term3fgHalf = (firstMoment'*sigmaInv_fg)';                          % 3 x pixel
    term3fg = 0.5*sum(term3fgHalf .* firstMoment);
    term3bgHalf = (firstMoment'*sigmaInv_bg)';                          % 3 x pixel
    term3bg = 0.5*sum(term3bgHalf .* firstMoment);
    term3 =  term3fg.*alpha_fg_idx + term3bg.*alpha_bg_idx;                               % 1 x pixel
    D(c, :) = term1+term2+term3;
end
computeDinComputeK_time = toc(start)
end



function k = initK(im_data_vectorized, alpha, bbox_vectorized)       % D is [k x numpixels]
initGlobalVariables;
k = zeros(1, size(im_data_vectorized,2));
if (allowCache)
    try
        load([kCacheFile]);
        return
    catch
    end
end
disp('initializing k using kmeans');
[kFg, ~] = kMeans(im_data_vectorized(:,bbox_vectorized==1)', K);                                   % k is [1 x numpixel]
[kBg, ~] = kMeans(im_data_vectorized(:,bbox_vectorized==0)', K);                                   % k is [1 x numpixel]
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