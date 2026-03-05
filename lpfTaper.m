function [clpf] = lpfTaper(c,dt,fcut,shift)
%
%--------------------------------------------------------------------------
%
% Low-pass filter with:
%
%     1. A linear trend removal
%
%     2. An optional shift in time
%
%     3. A tapered frequency response
%
%--------------------------------------------------------------------------
%
% Fundamental parameters
%
n        = numel(c);
df       = 1/(n*dt);
f        = df*(-n/2:-1+n/2);

%--------------------------------------------------------------------------
%
% Remove linear trend, apply time shift and FFT
%
clin     = c;

for k = 1:n
   clin(k) = c(1)+(c(end)-c(1))*(k-1)/(n-1);
end

clpf     = fft(c-clin).*fftshift(exp(2i*pi*shift*(-n/2:-1+n/2)/n));

%--------------------------------------------------------------------------
%
% Create mask function to achieve a tapered frequency response
%
cmask    = ones(size(c));

for k = 1:n
   if abs(f(k)) < 0.75*fcut
      cmask(k) = 1;
   elseif abs(f(k)) > 1.25*fcut
      cmask(k) = 0;
   else
      cmask(k) = 1-(abs(f(k))-0.75*fcut)/(0.5*fcut);
   end
end
%
%--------------------------------------------------------------------------
%
% Add back the linear term, multiply by mask and inverse FFT to get the
% final answer

clpf     = clin+ifft(clpf.*fftshift(cmask));
