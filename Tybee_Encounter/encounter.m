
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

% Oceanographic deep water wave theory: Wave length is 1.56*T^2, where 1.56 = g/(2*pi),
% g = gravity = 9.8 m/s^2

T            = 5.0;                                % Guess of wave period (seconds)
    
L            = 1.56*T^2;                           % Implied wavelength (m)

V_wave       = L/T;                                % Implied phase speed (m/s)

%-----------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------

% Aspect angle and radial velocities are from analysis of the GPS and ISAR data
AspGPS       = [-79.3, -16.8,   36.9,   69.9, 102.2, 131.8, 171.9, -145.8];
VrGPS        = [ 1.02,  5.03,   4.37,   1.92, -1.23, -4.00, -5.70,  -4.67];

% Heading, Bearing and speed for the Tybee from the GPS data
HdgGPS       = [ -10.3, -52.1, -103.9, -145.7, 169.9,  127.1,   76.9,   33.3];
BearBar      = [ -89.6, -68.9,  -67.0,  -75.8, -87.9, -101.0, -111.2, -112.4];
SpdGPS       = [  5.50,  5.27,   5.46,   5.60,  5.83,   6.03,   5.77,   5.65];

% Observed encounter period (seconds-per-crest) and wave amplitude (m/s)
% from the Chapeau filter used in the ISAR AutoTrack code
TChapeau     = [  5.21,  3.15,  2.64,  4.42,  6.40, 13.41, 18.20,  8.78];
RMSChapeau   = [ 0.261, 0.171, 0.085, 0.215, 0.147, 0.197, 0.154, 0.138];

% LOA from the IAT algorithm
LOA_IAT      = [ 35.5, 32.7, 32.9, 33.6, 27.7, 31.7, 27.4, 31.1];

% Mean rotation rates from the 3D and GPS algorithms
RotRate3D    = [ -0.1009, -0.0511, -0.0259,  0.0458,  0.0087, -0.2189,  0.2931, -0.1152];
RotRateGPS   = [  0.0708,  0.0164, -0.0266, -0.0350, -0.0354, -0.0294, -0.0065,  0.0354];

% 3D spectrum peak for tilt
TiltFreqHz   = [ 0.101, 0.220, 0.271, 0.271, 0.203, 0.085, 0.067, 0.101];
TiltPeriod   = 1./TiltFreqHz;

%-----------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------

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

figure(1);clf
hold on
xlabel('Wave Angle (deg)')
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

figure(2);clf
hold on
plot(Angle_Wave,1./sumvar)
xlabel('Wave Angle (deg)')
ylabel('Inverse of Summed Variance')
title ('Inverse of Summed Variance over 8 Tybee Cases ')
grid on; axis tight

figure(3);clf
scatter(TChapeau,encperiod(:,91),'b','filled')
hold on
redline = 0:1:20;
plot(redline,redline,'r')
xlabel('ISAR Wave Period (s)')
ylabel('39-meter Wave Model (s)')
title ('Fit of ISAR Wave Period Data to a 39-meter Ocean Wave')
grid on; axis tight

figure(4);clf
scatter(TChapeau,TiltPeriod,'b','filled')
xlabel('Peak Period from ISAR (s)')
ylabel('Tilt Period from 3D (s)')
hold on
redline = 0:1:18;
plot(redline,redline,'r')
grid on

figure(5);clf
scatter(7:14,RMSChapeau,'b','filled')
title ('Mocomp Amplitude vs Tybee Case')
xlabel('Tybee Case Number')
ylabel('RMSChapeau (m/s)')
grid on

figure(6);clf
scatter(HdgGPS,RMSChapeau,'b','filled')
title ('Wave Amplitude vs Heading')
xlabel('Heading (deg)')
ylabel('RMSChapeau (m/s)')
grid on
