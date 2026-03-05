function c = c1smooth(c,niters)

n = numel(c);

% Smooth complex time series niters times with a three point smoother

b = c;

k = 2:n-1;

if niters > 0

   for iter = 1:niters
        
      b(k) = 0.25*(c(k+1)+c(k-1)+2.0*c(k));
        
      c(k) = 0.25*(b(k+1)+b(k-1)+2.0*b(k));
                
   end
    
end
