function [D] = computeD(im_data_vectorized, numK, pi, mu, sigma)               % D_part is [k x numPixels]
% computes either fg or bg D, depending on if numK, pi, mu, sigma refers to
% fg or bg
D = zeros(numK, size(im_data_vectorized,2));
for c = 1:numK
    term1 = computeTerm1(pi, c);
    term2 = computeTerm2(sigma, c);
    term3 = computeTerm3(im_data_vectorized, mu, sigma, c);
    
    D(c, :) = term1+term2+term3;
end
end

function term1 = computeTerm1(pi, c)            % 1 x 1
term1 = -log(pi(c));
end

function term2 = computeTerm2(sigma, c)         % 1 x 1
term2 = 0.5*log(det(sigma(:,:,c)));
end

function term3 = computeTerm3(im_data_vectorized, mu, sigma, c)         % 1 x pixel
firstMoment = bsxfun(@minus, im_data_vectorized, mu(:,c));
sigmaInv = pinv(sigma(:,:,c));
term3half = (firstMoment'*sigmaInv)';
term3 = 0.5*sum(term3half .* firstMoment);
end
