function [pi mu sigma] = updateGMM(im_data, alpha, k, init)
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
    
    pi(fg_idx, i) = sum(idx_fg_i)/sum(idx_fg);
    pi(bg_idx, i) = sum(idx_bg_i)/sum(idx_bg);
    
    im_data_fg_i = im_data(:, idx_fg_i);
    im_data_bg_i = im_data(:, idx_bg_i);
    mu(:, fg_idx, i) = mean(im_data_fg_i,2);
    mu(:, bg_idx, i) = mean(im_data_bg_i,2);
    
    temp_fg_sigma_vector = bsxfun(@minus, im_data_fg_i, mu(:, fg_idx, i));
    sigma(:, :, fg_idx, i) =  temp_fg_sigma_vector * temp_fg_sigma_vector';
    
    temp_bg_sigma_vector = bsxfun(@minus, im_data_bg_i, mu(:, bg_idx, i));
    sigma(:, :, bg_idx, i) = temp_bg_sigma_vector * temp_bg_sigma_vector';
end
disp('step 2 done!');

% plot
%{
redFgMu = mu(1,1,:);
greenFgMu = mu(2,1,:);
redBgMu = mu(1,2,:);
greenBgMu = mu(2,2,:);
figure(1);
scatter(redFgMu, greenFgMu, 'r');
hold on;
scatter(redBgMu, greenBgMu, 'b');
%}
end
