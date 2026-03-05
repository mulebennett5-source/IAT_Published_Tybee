function [meanHeadingAcc,meanLengthAcc] = TrackAccuracyTest(icase,R,dr,PulseSNRdb,ShipAspectDeg,Sweep)

% Test code for track accuracy function

set(0,'DefaultAxesTitleFontWeight','bold')
set(0,'DefaultAxesFontWeight','bold')
set(0,'defaultLineLineWidth',2)

%--------------------------------------------------------------------------
%
% This program evaluates the track accuracy for three radar modes:
%
%    1. Circular Scan Radar (CSR), used by the system to identify targets
%       for ISAR imaging
%
%    2. ISAR Wiggle Mode (IWM), used during ISAR data collection to keep
%       the fat part of the beam on the target and to estimate the aspect
%       angle for 3-D ISAR processing
%
%    3. A 2-beam mode (TwoBeam) that alternates between two beam angles
%       for even and odd pulses. This mode would likely require an AESA.
%
%--------------------------------------------------------------------------

% Arguments to override the defaults for Range, Range Resolution,
% PulseSNRdb and Aspect
if nargin < 2 || isempty(R);             R             = 60000; end
if nargin < 3 || isempty(dr);            dr            = 5;     end
if nargin < 4 || isempty(PulseSNRdb);    PulseSNRdb    = 14;    end
if nargin < 5 || isempty(ShipAspectDeg); ShipAspectDeg = 45;    end
if nargin < 6 || isempty(Sweep);         Sweep         = 1;     end

% Default parameters - for all modes
D                    = 0.9;          % Antenna width (m)
lambda               = 0.03;         % Wavelength (m)
RadarPRF             = 512;          % Pulse repetition frequency (Hz)
Nscatterers          = 10;           % Number of scatterers
Npresum              = 4;            % Pre-sum factor
ShipSpeed            = 7.716;        % Ship Speed (m/s) [15 knots]

% The conceptual model for SNR is that there are 10 major scatterers on the
% target ship and that each scatterer has a single-pulse SNR of 14 dB.
% There is also some type of pre-summing done to sample the beamwidth
% properly (perhaps four effective samples per beam). The factor of 4 is
% based on the CSR scan rate of 180 deg/sec or about 90 beamwidths/second.
% At a PRF of 512 Hz this yields 5.7 samples per beamwidth. Allowing for
% some small oversampling, we can call this 4.
%
%--------------------------------------------------------------------------

%     Results April 25, 2020 for CSR and IWM


% CSR mode is limited by azimuth sample error

% Scan Rate (degree/sec):           180.00000
% Azimuth Sample (degree):            0.35156
% Total Time:                       200.00000
% Sample Time:                        2.00000
% Total Samples:                    100.00000
% Azimuth Error (degree):             0.10149
% 
% Mean heading accuracy (deg):        0.88083
% Mean length accuracy (percent):     1.53587


% ISAR Wiggle mode is limited by SNR Azimuth Error

% Scan Rate (degree/sec):            64.00000
% Azimuth Sample (degree):            0.12500
% Total Time:                        28.00000
% Sample Time:                        0.03125
% Total Samples:                    896.00000
% Azimuth Error (degree):             0.04218
% 
% Mean heading accuracy (deg):        0.87341
% Mean length accuracy (percent):     1.52634

% Factor that favors CSR Mode
% (28/200)^2      (Track Time)        0.01960

% Factors that favor ISAR Wiggle Mode
% 0.35156/0.125   (Az sampling)       2.81250
% 896/100         (Total samples)     8.96000
% 0.10149/0.04228 (Azimuth Error)     2.40004
% All factors favoring Wiggle Mode   60.49570

% All factors plus and minus          1.18570
% 2.40004/2.81250 (switch Az limits)  1.01180

% The bottom line is that the Wiggle Mode has two advantages - Improved
% Azimuth Accuracy and More Samples

% These two factors are roughly compensated by the larger Track Time of
% the CSR mode

%--------------------------------------------------------------------------

%     Results April 27, 2020 for Two-Beam Mode

% Two-Beam Mode Results - parameters chosen to get comparable results
% to CSR and IWM

% Scan Rate (degree/sec):            64.00000
% Azimuth Sample (degree):            0.12500
% Total Time:                        28.00000
% Sample Time:                        0.01563
% Total Samples:                         1792
% Azimuth Error (degree):             0.06310
% 
% Mean heading accuracy (deg):        0.92309
% Mean length accuracy (percent):     1.61674

%--------------------------------------------------------------------------
%
if     icase == 1

   % Default parameters - for system circular scan mode
   TrackTime         = 200;          % Tracking Time (s)
   
   RevisitTime       = 2*Sweep;      % Revisit Time (s)
   scanWidthDeg      = 360*Sweep;    % Width for sector scan

   % Disable unused parameters for this mode
   WiggleFreq        = [];           % Wiggle frequency (Hz)
   WiggleWidthDeg    = 360;          % Wiggle width (Deg)
   TwoBeamDiffDeg    = [];           % Angle difference between beams (Deg)
   AmpPerDeg         = [];           % Derivative of RCS wrt angle offset
   
   modename          = 'SystemCSR';  % Mode name

elseif icase == 2

   % Default parameters - for ISAR Wiggle Mode
   TrackTime         = 28;           % Tracking Time (s)
   
   WiggleFreq        = 16;           % Wiggle frequency (Hz)
   WiggleWidthDeg    = 2;            % Wiggle width (Deg)
   
   % Disable unused parameters for this mode
   RevisitTime       = [];           % Revisit Time (s)
   scanWidthDeg      = [];           % Width for sector scan
   TwoBeamDiffDeg    = [];           % Angle difference between beams (Deg)
   AmpPerDeg         = [];           % Derivative of RCS wrt angle offset

   modename          = 'Wiggle';     % Mode name

elseif icase == 3

   % Default parameters - for 2-Beam Mode
   TrackTime         = 28;           % Tracking Time (s)
   
   TwoBeamDiffDeg    = 1;            % Angle difference between beams (Deg)
   AmpPerDeg         = 0.5;          % Derivative of RCS wrt angle offset
   
   % Disable unused parameters for this mode
   RevisitTime       = [];
   scanWidthDeg      = [];           % Width for sector scan
   WiggleFreq        = [];
   WiggleWidthDeg    = [];

   modename          = 'TwoBeam';    % Mode name

end

rng('default');                      % Reset random number generator

%--------------------------------------------------------------------------
%
% Course Accuracy versus Range
%
nR                = 20;              % Output dimension
Rtest             = zeros(1,nR);
HeadingR          = zeros(1,nR);
for k = 1:nR
   Rtest(k)      = 12500*k;
   [HeadingR(k)] = TrackAccuracy(D,lambda,dr,Rtest(k),TrackTime,RevisitTime,    ...
                                 WiggleFreq,WiggleWidthDeg,RadarPRF,PulseSNRdb, ...
                                 Nscatterers,Npresum,TwoBeamDiffDeg,AmpPerDeg,  ...
                                 scanWidthDeg,ShipSpeed,ShipAspectDeg,icase);
end

figure(1001);clf
plot(Rtest/1000,HeadingR,'b')
title (['Course Accuracy vs Range at Aspect ' sprintf('%i',round(ShipAspectDeg)) ' degrees'])
xlabel('Range (km)')
ylabel('Course Accuracy (deg)')
axis tight; grid on; hold off
saveas(gcf,[modename '_CourseAccuracy_vs_Range'],'png')

%--------------------------------------------------------------------------
%
% Course Accuracy versus Ship Course
%
nH                = 17;
Htest             = zeros(1,nH);
HeadingH          = zeros(1,nH);
for k = 1:nH
   Htest(k)      = 5*k;
   [HeadingH(k)] = TrackAccuracy(D,lambda,dr,R,TrackTime,RevisitTime,           ...
                                 WiggleFreq,WiggleWidthDeg,RadarPRF,PulseSNRdb, ...
                                 Nscatterers,Npresum,TwoBeamDiffDeg,AmpPerDeg,  ...
                                 scanWidthDeg,ShipSpeed,Htest(k),icase);
end

figure(1002);clf
plot(Htest,HeadingH,'b')
title (['Aspect Accuracy vs Aspect at ' sprintf('%i',round(R/1000)) ' km Range'])
xlabel('Ship Aspect (deg)')
ylabel('Aspect Accuracy (deg)')
axis tight; grid on; hold off
saveas(gcf,[modename '_AspectAccuracy_vs_Aspect'],'png')

LengthAcc  = (100*pi/180)*HeadingH.*tand(Htest);
LengthAcc0 = (100*pi/180)*3*tand(Htest);

figure(1003);clf
plot(Htest,LengthAcc,'b')
hold on
plot(Htest,LengthAcc0,'r')
title (['Length Accuracy vs Aspect at ' sprintf('%i',round(R/1000)) ' km Range'])
xlabel('Ship Aspect (deg)')
ylabel('Length Accuracy (Percent)')
legend('Full Calculation','3-Degree Model')
axis tight; grid on; hold off
saveas(gcf,[modename '_LengthAccuracy_vs_Aspect'],'png')

meanHeadingAcc = mean(HeadingH);
meanLengthAcc  = mean(LengthAcc);

fprintf(1,'\nMean heading accuracy (deg):     %10.5f',  meanHeadingAcc);
fprintf(1,'\nMean length accuracy (percent):  %10.5f\n',meanLengthAcc);
