function [HeadingAccuracy] = TrackAccuracy(D,lambda,dr,R,TrackTime,RevisitTime,   ...
                                           WiggleFreq,WiggleWidthDeg,RadarPRF,    ...
                                           PulseSNRdb,Nscatterers,Npresum,        ...
                                           TwoBeamDiffDeg,AmpPerDeg,scanWidthDeg, ...
                                           ShipSpeed,ShipHeadingDeg,icase)

%--------------------------------------------------------------------------
%
% Sample Parameters for Circular Scan Radar (CSR) or ISAR Wiggle Mode (IWM)
% or Two-Beam Mode
% 
% Antenna width:      0.9 m
% 
% Wavelength:         0.03 m
% 
% Beamwidth:          0.03/0.9 radians (~2 degrees)
% 
% Range resolution:   5 m
% 
% Mean range          60000 m
%
% Tracking time:      200 sec for CSR, 28 sec for IWM
% 
% Sample time:        2 sec (0.5 Hz) for CSR, 1/32 sec (32 Hz) for IWM,
%                     2 pulses for Two-Beam Mode
% 
% RadarPRF:           512 Hz (PRI ~0.002 Second)
% 
% PulseSNRdb:         14 dB
%
% Nscatterers         10 (10 dB)
%
% Npresum              4 (6 dB)
%
% Effective SNR for entire target ~30 dB
% 
% Ship speed/heading: 15 knots, 45T degree steady course
% 
% Assume the look direction of the radar is 0 degrees so that the North
% component is always in-range and the East component is always cross-range
%
% This ship and platform geometry is exactly correct for the case where the
% aircraft's East component of velocity matches that of the ship. It is
% approximately correct for a range much longer than the ship track. One
% could generalize this approximation by taking into account a general
% model of the aircraft track and correcting for the aircraft pointing. In
% this genearl model the heading is the sum of the aspect angle to the
% radar and the pointing angle. The heading is presumed to be constant
% while the aspect and pointing angles have the same time dependence but
% opposite signs.
%
% Source:
%
% Radar 2009 A_15 Parameter Estimation and Tracking Part 1.pdf
%
% Dr. Robert M. O'Donnell
% IEEE New Hampshire Section Guest Lecturer
% January 1, 2010
%
% This analysis follows Merrill Skolnik's 'Radar Handbook'
%
%--------------------------------------------------------------------------

if icase     == 1
   ScanWidth  = 2*pi*scanWidthDeg/360;
   SampleTime = RevisitTime;
elseif icase == 2
   ScanWidth  = WiggleWidthDeg*pi/180;
   SampleTime = 1/(2*WiggleFreq);
elseif icase == 3
   ScanWidth  = TwoBeamDiffDeg*pi/180;
   SampleTime = 2*Npresum/RadarPRF;
end

%--------------------------------------------------------------------------

BeamWidth       = lambda/D;                            % 3-dB Beamwidth (radians)

PulseSNR        = 10^(0.1*PulseSNRdb);
SNR             = Nscatterers*Npresum*PulseSNR;        % Signal to Noise Ratio

% Range error is the SNR-enhanced basic resolution of the radar

errorR          = dr/sqrt(SNR);                        % Range accuracy (m)

% For cross-range position error use the maximum of the SNR-enhanced value
% of the 3-dB beamwidth and the error due to azimuth sampling

% Theoretical azimuth accuracy versus SNR (Skolnik)
if icase == 1 || icase == 2
   accuracyAz = 0.7*BeamWidth/sqrt(SNR);               % Azimuth accuracy (radians)
else
   accuracyAz = (pi/180)/(sqrt(SNR)*AmpPerDeg);
end

ScanRate        = ScanWidth/SampleTime;                % Scan rate (radians/second)

sampleAz        = ScanRate/RadarPRF;                   % Azimuth sampling (radians)

sampleAzErr     = sampleAz/sqrt(12);

errorAz         = max(accuracyAz,sampleAzErr);         % Azimuth error (radians)

errorCR         = R*errorAz;                           % Cross-range accuracy (m)

%--------------------------------------------------------------------------

Nsamples        = floor(TrackTime/SampleTime);         % Samples per track

% fprintf(1,'\nScan Rate (degree/sec):          %10.5f'  ,ScanRate*180/pi);
% fprintf(1,'\nAzimuth Sample (degree):         %10.5f'  ,sampleAz*180/pi);
% fprintf(1,'\nTotal Time:                      %10.5f'  ,TrackTime);
% fprintf(1,'\nSample Time:                     %10.5f'  ,SampleTime);
% fprintf(1,'\nTotal Samples:                   %10i'    ,Nsamples);
% fprintf(1,'\nAzimuth Error (degree):          %10.5f\n',errorAz *180/pi);

time            = (0:Nsamples-1)*SampleTime;           % Pulse times

% True values for East and North displacements
ShipEast        = ShipSpeed*sind(ShipHeadingDeg)*time; % East component of ship track
ShipNorth       = ShipSpeed*cosd(ShipHeadingDeg)*time; % North component of ship track

% Realizations for the Monte Carlo simulation
nReal           = 10000;

% Velocity components for realizations
inRangeV        = zeros(1,nReal);
crossRangeV     = zeros(1,nReal);

for k = 1:nReal
   
   % Calculate each realization
   inRange        = ShipNorth+errorR *randn(1,Nsamples);
   crossRange     = ShipEast +errorCR*randn(1,Nsamples);
   
   % Estimate the velocity components for each realization
   inRangeV(k)    = regression(time,inRange);
   crossRangeV(k) = regression(time,crossRange);

   % Plot the displacements for each realization
   % figure(1000);clf
   % plot(time,inRange,'b')
   % hold on
   % plot(time,crossRange,'r')
   % xlabel('Time (s)')
   % ylabel('Displacement (m)')
   % legend('InRange','CrossRange')
   % axis tight; grid on; hold off

end

% Estimated ship heading for the realizations
ShipHeadingEst  = atan2d(crossRangeV,inRangeV);

% Mean and standard deviation over the realizations
HeadingMean     = mean(ShipHeadingEst);
HeadingAccuracy = std(ShipHeadingEst);

ShipSpeedEst    = sqrt(crossRangeV.^2+inRangeV.^2);

% Mean and standard deviation over the realizations
SpeedEstMean    = mean(ShipSpeedEst);
SpeedEstStd     = std(ShipSpeedEst);

% fprintf(1,'\nHeading - True, Estimated Mean and Std. Dev.: %10.5f %10.5f %10.5f',   ...
%              ShipHeadingDeg,HeadingMean,HeadingAccuracy);% Estimated ship heading for the realizations
% 
% fprintf(1,'\nSpeed - True, Estimated Mean and Std. Dev.:   %10.5f %10.5f %10.5f\n', ...
%              ShipSpeed,SpeedEstMean,SpeedEstStd);
