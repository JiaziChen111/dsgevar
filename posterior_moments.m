function [post_mean,HDImin,HDImax] = posterior_moments(xx,mh_conf_sig)

sortedVec = sort(xx);

post_mean = mean(sortedVec);

ciIdx = ceil(mh_conf_sig*length(sortedVec));
nCIs = length(sortedVec) - ciIdx;  % number of vector elements that make HDI

% Determine middle of HDI to get upper and lower bound
ciWidth = zeros(nCIs,1); 
for ind = 1:nCIs 
    ciWidth(ind) = sortedVec(ind + ciIdx) - sortedVec(ind);
end
    [~,idxMin] = min(ciWidth);
    HDImin = sortedVec(idxMin);
    HDImax = sortedVec(idxMin + ciIdx);
end
 
