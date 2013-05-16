function [pi mu sigma] = updateGMM(im_data_vectorized, alpha, alphaVal, k, numK, numColors, init)
% disp('step 2: finding theta');
initGlobalVariables;
pi = zeros(1, numK);                           % pi is [1 x k]
mu = zeros(numColors, numK);                   % mu is [3 x k]
sigma = zeros(numColors, numColors, numK);     % sigma is [3 x 3 x k]
for i = 1:numK
    idx_k = k == i;
    idx_alpha = alpha==alphaVal;
    idx_alpha_k = idx_alpha & idx_k;
    
    numPixel_alpha_k = sum(idx_alpha_k);
    numPixel_alpha = sum(idx_alpha);
    
    % compute pi
    pi(i) = numPixel_alpha_k/numPixel_alpha;
    %{
    if (numPixel_alpha_k == 0 || numPixel_alpha == 0 || abs(pi(i)) < epsilon)
        pi(i) = epsilon;
    end
    %}
    
    % comptue mu
    im_data_alpha_k = im_data_vectorized(:, idx_alpha_k);
    if numPixel_alpha_k == 0
        mu(:, i) = zeros(numColors, 1);
    else
        mu(:, i) = mean(im_data_alpha_k,2);
    end
    
    % compute sigma
    %if numPixel_alpha_k <= 1
    %    sigma(:, :, i) = eye(numColors)*delta;
    %else
        first_central_moment_full = bsxfun(@minus, im_data_vectorized, mu(:, i));
        first_central_moment = first_central_moment_full(:, idx_alpha_k);
        sigma(:, :, i) =  (first_central_moment * first_central_moment')/(numPixel_alpha_k-1) + eye(numColors)*delta;
    %end
end
% disp('step 2 done!');


end
