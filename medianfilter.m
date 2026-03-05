function out = medianfilter(d,n)

% d: Values to be filtered
% n: Maximum excursion for filter

% Apply a 1-dimensional median filter with maximum excursion n, conditional
% on outlier determination

M      = numel(d);
out    = d;

if n == 0; return; end
   
% Find an estimate of the local median using the maximum excursion
out    = zeros(1,M);
for k = 1:M
   kfirst = max(1,k-n);
   klast  = min(M,k+n);
   out(k) = median(d(kfirst:klast));
end
 