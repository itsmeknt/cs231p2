function [alpha] = updateAlpha(k, pi, mu, alph, gamma, beta, init)
disp('step 3: finding alpha');
alpha = zeros(1, length(k));
disp('step 3 done!');
end


function D = computeD(im_data, k, pi, mu, sigma)       % D is [numAlpha x numpixels]
initGlobalVariables;
numPixels = length(k);
D = zeros(numAlphaValues, numPixels);

% cache computations
term1cache = -log(pi);                                         % a x K
term2cache = zeros(numAlphaValues, K);                         % a x K
sigmaInvCache = zeros(numColors, numColors, numAlphaValues, K);      % 3 x 3 x a x K
for i=1:K
    term2cache(fg_idx, i) = 0.5*log(det(sigma(:,:,fg_idx, i)));
    term2cache(bg_idx, i) = 0.5*log(det(sigma(:,:,bg_idx, i)));
    
    sigmaInvCache(:, :, fg_idx, i) = pinv(sigma(:,:,fg_idx,i));
    sigmaInvCache(:, :, bg_idx, i) = pinv(sigma(:,:,bg_idx,i));
end


aa = tic;
for a = 1:numAlphaValues
    term1 = zeros(1, numPixels);
    term2 = zeros(1, numPixels);
    term3 = zeros(1, numPixels);
    
    muForAllPixels = zeros(1, numPixels);
    for c = 1:K
        k_idx = k==c;
        term1 = term1 + term1cache(a,c)*k_idx;    % 1 x pixel
        term2 = term2 + term2cache(a,c)*k_idx;    % 1 x pixel
        
        muForAllPixels = muForAllPixels + bsxfun(@times, mu(:,a,c), k_idx);    % 3 x pixel
    end
    
    
        firstMoment = bsxfun(@minus, im_data, mu(:,a,c), alpha==fg_val)im_data - muForAllPixels;                                          % 3 x pixel
        
        term3fgHalf = firstMoment'*sigmaInvCache(:,:,fg_idx,c);                          % pixel x 3
        term3fg = 0.5*sum(term3fgHalf' .* firstMoment);
        term3bgHalf = firstMoment'*sigmaInvCache(:,:,bg_idx,c);                          % pixel x 3
        term3bg = 0.5*sum(term3bgHalf' .* firstMoment);
        term3 =  term3fg.*(alpha==fg_val) + term3bg.*(alpha==bg_val);                               % 1 x pixel
    
    D(c, :) = term1+term2+term3;
end
bb = toc(aa)
end