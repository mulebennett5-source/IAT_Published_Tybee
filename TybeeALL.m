
set(0,'DefaultAxesTitleFontWeight','bold')
set(0,'DefaultAxesFontWeight','bold')
set(0,'defaultLineLineWidth',2)

%-----------------------------------------------------------------------------------------------
%
% Compute the encounter frequencies for 8 cases of the USCG Tybee
% observed with the SCATR radar.
%
% We have estimates of these values from the adaptive motion compensation
% algorithm of the ISAR code. This calculation just tries to fit this
% information to an oceanographic model with a constant period deep water
% gravity wave at an angle to the ship's course for 8 cases where the
% Tybee was on an octagonal course.
%
% The period of 5.0 is a guess based on running the code for various
% values and comparing the encounter periods with those from the ISAR code.
%
%-----------------------------------------------------------------------------------------------

% Oceanographic deep water wave gravtheory: Wave length is 1.56*T^2,
% where 1.56 = g/(2*pi), g = gravity = 9.8 m/s^2

T            = 5.0;                                % Guess of wave period (seconds)
    
L            = 1.56*T^2;                           % Implied wavelength (m)

V_wave       = L/T;                                % Implied phase speed (m/s)

%-----------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------
%
% Mean aspect angle and initial radial velocities used in the 3D ISAR program and the ISAR processor 
% The aspect angle is for an average of 4 minutes after skipping the first and last 2 minutes of each
% octagon leg
%
% The initial velocity is the start mocomp velocity after skipping the first and last 2 minutes of each
% octagon leg
%
%                       Produced from the GPS data on October 13, 2025
%
AspGPS       = [ -72.1, -14.2,   34.0,   67.9,  99.5,  131.2,  168.8, -143.7];  % 244 second version
VrGPS        = [  0.75,  5.05,   4.24,   2.03, -1.13,  -3.96,  -6.34,  -4.49];
%
%-----------------------------------------------------------------------------------------------

% Heading and speed for the Tybee from the GPS data (not sure about dwell times)
HdgGPS       = [ -10.3, -52.1, -103.9, -145.7, 169.9,  127.1,   76.9,   33.3];
BearBar      = [ -89.6, -68.9,  -67.0,  -75.8, -87.9, -101.0, -111.2, -112.4];
SpdGPS       = [  5.50,  5.27,   5.46,   5.60,  5.83,   6.03,   5.77,   5.65];

% Observed encounter period (seconds-per-crest) and wave amplitude (m/s)
% from the Chapeau filter used in the ISAR AutoTrack code (60 second dwells)
TChapeau     = [  5.21,  3.15,  2.64,  4.42,  6.40, 13.41, 18.20,  8.78];
RMSChapeau   = [ 0.261, 0.171, 0.085, 0.215, 0.147, 0.197, 0.154, 0.138];

% LOA from the IAT algorithm - From Focus3D with 244 seconds data and aspect angle from AspGPS above
LOA_IAT      = [ 34.5, 32.2, 30.2, 36.8, 36.0, 31.1, 31.0, 31.1];

% Mean rotation rates from the 3D and GPS algorithms
RotRate3D    = [ -.1009, -.0511, -.0259,  .0458,  .0087, -.2189,  .2931, -.1152];
RotRateGPS   = [  .0520,  .1045,  .0135, -.0362, -.0370, -.0145,  .1064, -.0241];
RotRateHdg   = [  .0189, -.0881, -.0404, -.0009,  .0015, -.0158, -.1125,  .0603];
RotRateBear  = [  .0708,  .0164, -.0268, -.0353, -.0354, -.0303, -.0060,  .0362];

% 3D spectrum peak for tilt
TiltFreqHz   = [ 0.101, 0.220, 0.271, 0.271, 0.203, 0.085, 0.067, 0.101];
TiltPeriod   = 1./TiltFreqHz;

%-----------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------

numfiles  = 8;            % SCATR Data: 8 each cases for Point Stuart and Tybee

mocfiles  = dir('*.moc'); % Motion compensation data from the ISAR processor

pulseskip = 32;           % Skip the acquisition mode pulses used by the ISAR processor

for k = 1:numfiles

   fprintf(1,['\n\n ' mocfiles(k).name '\n']);
   
   filename      = mocfiles(k).name;

   % Get the parameters from the GPS routine that apply to all methods

   [RdotBar(k),Rbar(k),SpeedBar,VgpsR,VgpsC,AspBar(k),HdgBar,tgpsV,moc32,tmoc,dtmoc] = GetSCATR_GPS(filename,k,pulseskip);

end

% Set up arrays for search over the ocean wave direction of propagation
Angle_Wave = zeros(1,361);
vrelw      = zeros(8,361);
encfreq    = zeros(8,361);
encperiod  = zeros(8,361);
evar       = zeros(8,361);
sumvar     = zeros(1,361);

for k = 1:361               % Wave angle, 0-360 degrees by 1 degree increments

   Angle_Wave(k) = (k-1)*1;
   
   %
   %   The encounter frequency consists of two terms:
   %
   %   1. 1/T is the rate of encounter of crests for a stationary object
   %
   %   2. vrel/L is the rate that the ship crosses the wave crests due to its own motion
   %
   
   for m = 1:8
      
      vrelw(m,k)     = SpdGPS(m)*cosd(Angle_Wave(k)-HdgGPS(m));
      encfreq(m,k)   = (1/T)-vrelw(m,k)/L;
      encperiod(m,k) = 1/abs(encfreq(m,k));
   
      evar(m,k)      = (encperiod(m,k)-TChapeau(m))^2;

   end

end

fprintf(1,'\n RMS values of GPS, Hdg, Bearing, 3D Rotation Rates:   %10.4f %10.4f %10.4f %10.4f\n', ...
        rms(RotRateGPS),rms(RotRateHdg),rms(RotRateBear),rms(RotRate3D));

figure(401);clf
hold on
xlabel('WaveAngle (deg)')
ylabel('Variance')
title ('Variance for 8 Tybee Cases versus Wave Angle')
grid on
xlim([0 360])
legend('location','north')
for m = 1:8
   plot(Angle_Wave,evar(m,:))
   for k = 1:361
      sumvar(k) = sumvar(k)+evar(m,k);
   end
end

figure(402);clf
hold on
plot(Angle_Wave,1./sumvar)
xlabel('WaveAngle (deg)')
ylabel('Inverse of Summed Variance')
title ('Inverse of Summed Variance over 8 Tybee Cases ')
grid on; axis tight

figure(403);clf
scatter(TChapeau,encperiod(:,91),'b','filled')
hold on
redline = 0:1:20;
plot(redline,redline,'r')
xlabel('ISAR Wave Period (s)')
ylabel('39-meter Wave Model (s)')
title ('Fit of ISAR Wave Period Data to a 39-meter Ocean Wave')
grid on; axis tight

figure(404);clf
scatter(TChapeau,TiltPeriod,'b','filled')
xlabel('Peak Period from ISAR (s)')
ylabel('Tilt Period from 3D (s)')
hold on
redline = 0:1:18;
plot(redline,redline,'r')
grid on

figure(405);clf
scatter(7:14,RMSChapeau,'b','filled')
title ('Mocomp Amplitude vs Tybee Case')
xlabel('Tybee Case Number')
ylabel('RMSChapeau (m/s)')
grid on

figure(406);clf
scatter(HdgGPS,RMSChapeau,'b','filled')
title ('Wave Amplitude vs Heading')
xlabel('Heading (deg)')
ylabel('RMSChapeau (m/s)')
grid on

figure(407);clf
scatter(RotRateGPS,RotRate3D,'b','filled')
hold on
title ('Rotation Rates for GPS and 3D')
xlabel('GPS Rotation Rate (deg/s)')
ylabel('3D Rotation Rate (deg/s)')
xylim  = max(max(abs(RotRateGPS)),max(abs(RotRate3D)));
xyline = -xylim:0.01*xylim:xylim;
plot(xyline,xyline,'r')
xlim([-xylim xylim])
ylim([-xylim xylim])
grid on

figure(408);clf
scatter(1:8,RotRateGPS,'b','filled')
hold on
scatter(1:8,RotRateHdg,'r','filled')
scatter(1:8,RotRateBear,'g','filled')
scatter(1:8,RotRate3D,'k','filled')
title ('Rotation Rates for GPS, Heading, Bearing and 3D')
xlabel('Case Number')
ylabel('Rotation Rate (deg/s)')
legend('GPS','Heading','Bearing','3D')
grid on

