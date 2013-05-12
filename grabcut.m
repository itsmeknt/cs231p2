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



% init k, GMM params to 0
k = zeros(1, length(alpha));
pi = zeros(numAlphaValues, K);
mu = zeros(numColors, numAlphaValues, K);
sigma = zeros(numColors, numColors, numAlphaValues, K);
oldPi = pi;
oldMu = mu;
oldSigma = sigma;

% starting step - decides which step we start the initialization process
startingStep = 1;

converge = false;
iter = 1;
skipStep = true;
beta = 0;
gamma = 0;
while true
    iter
    
    % step 1
    if (~skipStep || startingStep==updateKidx)
        k = updateK(im_data, alpha, bbox_vectorized, pi, mu, sigma, startingStep==updateKidx && iter == 1);
        if (startingStep == updateKidx)
            skipStep = false;
        end
    end
    
    % step 2
    if (~skipStep || startingStep==updateGMMidx)
        oldPi = pi;
        oldMu = mu;
        oldSigma = sigma;
        [pi mu sigma] = updateGMM(im_data, alpha, k, startingStep==updateGMMidx && iter == 1);
        
        if (iter > 1 && hasConverged(pi, mu, sigma, oldPi, oldMu, oldSigma))
            converge = true;
        else
            converge = false;
        end
        
        if (startingStep == updateGMMidx)
            skipStep = false;
        end
    end
    
    % step 3
    if (~skipStep || startingStep==updateAlphaIdx)
        alpha = updateAlpha(k, pi, mu, sigma, gamma, beta, startingStep==updateAlphaIdx && iter == 1);
        
        if (startingStep == updateAlphaIdx)
            skipStep = false;
        end
    end
    
    if (converge)
        break;
    end
    
    skipStep = false;
    iter = iter+1;
end
end

function converge = hasConverged(pi, mu, sigma, oldPi, oldMu, oldSigma)
disp('testing convergence')
distance = euclidean_distance(pi, mu, sigma, oldPi, oldMu, oldSigma)
converge = distance < 1;
disp('testing convergence done!')
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