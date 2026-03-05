function [corr] = correlation(c1,c2)

c1bar   = mean(c1(:));
c2bar   = mean(c2(:));
c1prime = c1-c1bar;
c2prime = c2-c2bar;

cov     = c1prime.*conj(c2prime);

varc1p  = c1prime.*conj(c1prime);
c1sigma = sqrt(mean(varc1p(:)));
varc2p  = c2prime.*conj(c2prime);
c2sigma = sqrt(mean(varc2p(:)));

corr    = mean(cov(:))/(c1sigma*c2sigma);
