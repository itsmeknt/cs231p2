function [ alphas ] = grabcut( im_name, numSamples)
initGlobalVariables;

seg_name = im_name;
periodIdx = find(im_name == '.');
if (~isempty(periodIdx))
    seg_name = [im_name(1:periodIdx(end)) 'bmp'];
end

true_seg_data = double(imread([seg_dir '/' seg_name]))/255;
true_seg_data_vectorized = true_seg_data(:)';


alphas = [];
for i = 1:numSamples;
    newAlpha = grabcutRun(im_name, i);
    if isempty(alphas)
        alphas = newAlpha;
    else
        alphas = alphas & newAlpha;
    end
end


score = computeScore(alphas(:)', true_seg_data_vectorized)
figure;
imshow(alphas);
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

