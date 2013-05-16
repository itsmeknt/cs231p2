function [alpha energy] = updateAlpha(im_data, bbox_vectorized, D_fg, D_bg, alpha, beta, gamma)
initGlobalVariables;
% disp('step 3: finding alpha');
% Dsliced = D(fg_idx,:).*(alpha==fg_val) + D(bg_idx,:).*(alpha==bg_val);
%tU = find(bbox_vectorized==1);
%D = D(:,tU);
% Dsliced = Dsliced(:,tU);
% alpha_fg_idx=alpha==fg_val;
% alpha_bg_idx=alpha==bg_val;
% D(1,:)=D(1,:).*alpha_fg_idx;
% D(2,:)=D(2,:).*alpha_bg_idx;

D_fg_vector = D_fg;
if (size(D_fg, 1) > 1)
    D_fg_vector = min(D_fg);
end

D_bg_vector = D_bg;
if (size(D_bg, 1) > 1)
    D_bg_vector = min(D_bg);
end

D = [D_fg_vector; D_bg_vector];
V = computeV(im_data, alpha, beta, gamma);

%w = size(im_data,2);
%h = size(im_data,1);
%figure;
%imagesc(reshape(V(1,:,1), h, w));
%figure;
%imagesc(reshape(V(2,:,1), h, w));
%pause;
%V = V(:,tU,:);
% maxV = max(sum(V(:,:,1)));
%alphatU = minCut(D, V, tU, numPixels, maxV+1);
%alpha(tU) = alphatU;

newAlpha = segment_full(D, V);

%{
Dviz = D(fg_idx,:).*(alpha==fg_val) + D(bg_idx,:).*(alpha==bg_val);
figure;
imagesc(reshape(Dviz, h, w));
DvizInv = D(fg_idx,:).*(alpha==bg_val) + D(bg_idx,:).*(alpha==fg_val);
figure;
imagesc(reshape(DvizInv, h, w));
DvizDiff = Dviz - DvizInv;
figure;
imagesc(reshape(DvizDiff.*(DvizDiff>0), h, w));
figure;
imagesc(reshape(DvizDiff.*(DvizDiff<0), h, w));
%}

bbox_vectorized_idx = find(bbox_vectorized==1);
alpha(bbox_vectorized_idx) = newAlpha(bbox_vectorized_idx);
energy = computeEnergy(D, V, alpha);
% disp('step 3 done!');
end

% function D=computeD(im_data_vectorized, alpha, k, pi, mu, sigma)
% 
% initGlobalVariables;
% 

% function D = computeD(im_data_vectorized, k, pi, mu, sigma,alpha)       % D is [numAlpha x numpixels]
% initGlobalVariables;
% numPixels = length(k);
% D = zeros(numAlphaValues, numPixels);
% alpha_fg_idx=alpha==fg_val;
% alpha_bg_idx=alpha==bg_val;
% 
% % cache computations
% term1cache = -log(pi);                                         % a x K
% 
% start = tic;
% for a = 1:numAlphaValues
%     term1 = zeros(1, numPixels);
%     term2 = zeros(1, numPixels);
%     term3 = zeros(1, numPixels);
%     
%     for c = 1:K
%         k_idx = k==c;
%         term1 = term1 + term1cache(a,c)*k_idx;    % 1 x pixel
%         
%         term2logDet = 0.5*log(det(sigma(:,:,a, c)));
%         term2 = term2 + term2logDet*k_idx;    % 1 x pixel
%         
%         muForAllPixels = bsxfun(@times, mu(:,a,c), k_idx);    % 3 x pixel
%         firstMoment = bsxfun(@times, im_data_vectorized - muForAllPixels, k_idx) ;                                          % 3 x pixel
%         sigmaInv = pinv(sigma(:,:,a,c));
%         term3half = (firstMoment'*sigmaInv)';                          % 3 x pixel
%         term3 = term3 + 0.5*sum(term3half .* firstMoment);
%     end
%     D(a, :) = term1+term2+term3;
% end
% computeDinComputeAlpha_time = toc(start)
% end

% function V=getV(im_data,alpha,beta,gamma)
% 
% 
% end

function V = computeV(im_data, alpha, beta, gamma)       % V is [numV x numpixels x 2]
initGlobalVariables;
V = zeros(size(Voffsets,1), size(im_data,1)*size(im_data,2), 2);
for i = 1:size(Voffsets,1)
    offset = Voffsets(i, :);
    V(i,:,:) = computeVrow(im_data, alpha, beta, gamma, offset);
end

end

function Vrow = computeVrow(im_data, alpha, beta, gamma, offset)
initGlobalVariables;
w = size(im_data, 2);
h = size(im_data, 1);
im_data_offset = offsetMatrix(im_data, offset, inf);                % pixels out of boundary will be infinity, so V at that entry will be 0
alphaMatrix = reshape(alpha, h, w);
alphaMatrixOffset = offsetMatrix(alphaMatrix, offset, NaN);         % pixels out of boundary will be NaN so equality test always fail
im_idx = reshape(1:h*w, h, w);
im_idx_offset = offsetMatrix(im_idx, offset, -1);                   % pixels out of boundary will be -1

diff = im_data - im_data_offset;
dist = sum(diff.*diff, 3);

weight = 1; % sqrt(sum(offset.*offset));
weightedExpDist = weight*gamma*exp(-1*beta*dist);

alphaDiffIdx = alphaMatrix ~= alphaMatrixOffset;
Vmatrix = alphaDiffIdx.*weightedExpDist;
Vrow1 = reshape(Vmatrix, 1, h*w);

Vrow = zeros(1, h*w, 2);
Vrow(:,:,1) = Vrow1;
Vrow(:,:,2) = reshape(im_idx_offset, 1, h*w);
end

function matrixOffset = offsetMatrix(matrix, offset, outOfBoundaryValue)        % offset should be [yoffset, xoffset] integers
offset_y_dist = abs(offset(1));
offset_x_dist = abs(offset(2));

matrixOffset = outOfBoundaryValue*ones(size(matrix,1),size(matrix,2), size(matrix,3));
if offset(1) >=0 && offset(2) >= 0
    matrixOffset(offset_y_dist+1:end, offset_x_dist+1:end, :) = matrix(1:end-offset_y_dist, 1:end-offset_x_dist, :);
elseif offset(1) < 0 && offset(2) >= 0
    matrixOffset(1:end-offset_y_dist, offset_x_dist+1:end, :) = matrix(offset_y_dist+1:end, 1:end-offset_x_dist, :);
elseif offset(1) >= 0 && offset(2) < 0
    matrixOffset(offset_y_dist+1:end, 1:end-offset_x_dist, :) = matrix(1:end-offset_y_dist, offset_x_dist+1:end, :);
elseif offset(1) < 0 && offset(2) < 0
    matrixOffset(1:end-offset_y_dist, 1:end-offset_x_dist, :) = matrix(offset_y_dist+1:end, offset_x_dist+1:end, :);
end
end

function matrix = stripOffset(matrixOffset, offset)
offset_y_dist = abs(offset(1));
offset_x_dist = abs(offset(2));

if offset(1) >=0 && offset(2) >= 0
    matrix = matrixOffset(1+offset_y_dist:end, 1+offset_x_dist:end);
elseif offset(1) < 0 && offset(2) >= 0
    matrix = matrixOffset(1:end-offset_y_dist, 1+offset_x_dist:end);
elseif offset(1) >= 0 && offset(2) < 0
    matrix = matrixOffset(1+offset_y_dist:end, 1:end-offset_x_dist);
elseif offset(1) < 0 && offset(2) < 0
    matrix = matrixOffset(1:end-offset_y_dist, 1:end-offset_x_dist);
end
end

function energy = computeEnergy(D, V, alpha) 
initGlobalVariables;
DatAlpha = D(fg_idx, :).*(alpha==fg_val) + D(bg_idx,:).*(alpha==bg_val);
energy = sum(DatAlpha) + sum(sum(V(:, :, 1)));
end

