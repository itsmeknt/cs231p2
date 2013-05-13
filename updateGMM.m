function [pi mu sigma] = updateGMM(im_data_vectorized, alpha, k, init)
disp('step 2: finding theta');
initGlobalVariables;
pi = zeros(numAlphaValues, K);                              % pi is [2 x k]
mu = zeros(numColors, numAlphaValues, K);                   % mu is [3 x 2 x k]
sigma = zeros(numColors, numColors, numAlphaValues, K);     % sigma is [3 x 3 x 2 x k]
for i = 1:K
    idx_k_i = k == i;
    idx_fg = alpha==fg_val;
    idx_fg_i = idx_fg & idx_k_i;
    idx_bg = alpha==bg_val;
    idx_bg_i= idx_bg & idx_k_i;
    
    numPixelFg_i = sum(idx_fg_i);
    numPixelFg = sum(idx_fg);
    pi(fg_idx, i) = numPixelFg_i/numPixelFg;
    if (numPixelFg_i == 0 || numPixelFg == 0 || abs(pi(fg_idx, i)) < epsilon)
        pi(fg_idx, i) = epsilon;
    end
    
    numPixelBg_i = sum(idx_bg_i);
    numPixelBg = sum(idx_bg);
    pi(bg_idx, i) = numPixelBg_i/numPixelBg;
    if (numPixelBg_i == 0 || numPixelBg == 0 || abs(pi(bg_idx, i)) < epsilon)
        pi(bg_idx, i) = epsilon;
    end
    
    im_data_fg_i = im_data_vectorized(:, idx_fg_i);
    im_data_bg_i = im_data_vectorized(:, idx_bg_i);
    
    if numPixelFg_i == 0
        mu(:, fg_idx, i) = zeros(numColors, 1, 1);
        sigma(:, :, fg_idx, i) = eye(numColors)*epsilon;
    else
        mu(:, fg_idx, i) = mean(im_data_fg_i,2);
        
        temp_fg_sigma_vector = bsxfun(@minus, im_data_fg_i, mu(:, fg_idx, i));
        sigma(:, :, fg_idx, i) =  temp_fg_sigma_vector * temp_fg_sigma_vector';
    end
    
    if numPixelBg_i == 0
        mu(:, bg_idx, i) = zeros(numColors, 1, 1);
        sigma(:, :, bg_idx, i) = eye(numColors)*epsilon;
    else
        mu(:, bg_idx, i) = mean(im_data_bg_i,2);
        
        temp_bg_sigma_vector = bsxfun(@minus, im_data_bg_i, mu(:, bg_idx, i));
        sigma(:, :, bg_idx, i) = temp_bg_sigma_vector * temp_bg_sigma_vector';
    end
    
    
    
    
end
disp('step 2 done!');


end
