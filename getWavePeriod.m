function [wavePeriod,aOut] = getWavePeriod(a,time,lowF,highF,fignum,figTitle)
%
%--------------------------------------------------------------------------
%
% Low-pass filter with:
%
%     1. A linear trend removal
%
%     2. Find best frequency
%
%--------------------------------------------------------------------------
%
% Fundamental parameters
%
n          = numel(a);             % Number of time steps
dt         = mean(diff(time));     % Mean time step
df         = 1/(n*dt);             % Frequency step for normal FFT

%--------------------------------------------------------------------------
%
% Remove linear trend
%
alin       = a;

for k = 1:n
   alin(k) = a(1)+(a(end)-a(1))*(k-1)/(n-1);
end

aDetrend   = a-alin;

aScoreF    = zeros(1,n);
F          = (df/8)*(0:n-1);

kmin       = round(lowF/(df/8));
kmax       = min(n,round(highF/(df/8)));

for k = kmin:kmax
   f          = F(k);
   aCos       = aDetrend.*cos(2*pi*f*time);
   aSin       = aDetrend.*sin(2*pi*f*time);
   aScoreF(k) = mean(aCos).^2+mean(aSin).^2;
end

if fignum > 0
   figure(2*fignum-1);clf
   plot(F,10*log10(aScoreF+1e-6*max(aScoreF)),'b')
end

% First pass - use the peak frequency
[~,kFreq]  = max(aScoreF);
waveFreq   = F(kFreq);

dfreq      = 0.1*waveFreq;
kmin       = max(1,round((waveFreq-dfreq)/(df/8)));
kmax       = min(n,round((waveFreq+dfreq)/(df/8)));

aScoreF    = zeros(1,n);
for k = kmin:kmax
   f          = F(k);
   aCos       = aDetrend.*cos(2*pi*f*time);
   aSin       = aDetrend.*sin(2*pi*f*time);
   aScoreF(k) = mean(aCos)^2+mean(aSin)^2;
end

% Second pass - use a weighted version near the peak
waveFreq   = sum(aScoreF.*F)/sum(aScoreF);

% Diagnostic - estimate the output values according to the peak model
aOut       = alin+2*mean(aDetrend.*cos(2*pi*waveFreq*time),2)*cos(2*pi*waveFreq*time)+ ...
                  2*mean(aDetrend.*sin(2*pi*waveFreq*time),2)*sin(2*pi*waveFreq*time);
               
% Convert to wave period
wavePeriod = min(1/(0.000001+waveFreq),n*dt);

if fignum > 0
   
   fprintf(1,['\n Wave Period Estimate (sec) for ' figTitle ' %8.4f\n'],wavePeriod);

   hold on
   plot(F,10*log10(aScoreF+1e-6*max(aScoreF)),'r')
   xlabel('Frequency (Hz)')
   ylabel('Power (dB)')
   title (figTitle);
   legend('First Pass','Second Pass')
   xlim([0 highF])
   grid on; hold off

   aFit = aOut-alin;
   
   figure(2*fignum);clf
   plot(time,a,'b')
   hold on
   plot(time,aFit,'r')
   plot(time,a-aFit,'g')
   xlabel('Time (s)')
   ylabel('Original and Fit')
   title (figTitle);
   legend('Original','Data fit','Diff')
   grid on; hold off; axis tight

end