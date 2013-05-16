function [k] = updateK(im_data_vectorized, alpha, alphaVal, bbox_vectorized, numK, pi, mu, sigma, init, rng_seed )
% disp('step 1: finding k');
if init
    k = initK(im_data_vectorized, alpha, alphaVal, bbox_vectorized, numK, rng_seed);
else
    D = computeD(im_data_vectorized, numK, pi, mu, sigma);
    [~, kUpdated] = min(D);
    k = kUpdated;
end
% disp('step 1 done!');
end

function k = initK(im_data_vectorized, alpha, alphaVal, bbox_vectorized, numK, rng_seed)       % D is [k x numpixels]
initGlobalVariables;
rng(rng_seed);

if (strcmp(initType, 'random'))
    k = randi(numK, [1 size(im_data_vectorized, 2)]);
elseif (strcmp(initType, 'kmeans'))
    %if (allowCache)
    %    try
    %        load([kCacheFile]);
    %        return
    %    catch
    %    end
    %end
    
    % disp('initializing k using kmeans');
    numPixels = size(im_data_vectorized, 2);
    [kFg, ~] = kMeans(im_data_vectorized(:,bbox_vectorized==1)', numK);                                   % k is [1 x numpixel]
    [kBg, ~] = kMeans(im_data_vectorized(:,bbox_vectorized==0)', numK);                                   % k is [1 x numpixel]
    k = zeros(1, numPixels);
    fgIdx = 1;
    bgIdx = 1;
    for i=1:numPixels
        if (alpha(i) == fg_val)
            k(i) = kFg(fgIdx);
            fgIdx = fgIdx+1;
        elseif (alpha(i) == bg_val)
            k(i) = kBg(bgIdx);
            bgIdx = bgIdx+1;
        end
    end
    % save(kCacheFile, 'k');
end
end