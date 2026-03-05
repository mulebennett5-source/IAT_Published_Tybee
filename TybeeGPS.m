
set(0,'DefaultAxesTitleFontWeight','bold')
set(0,'DefaultAxesFontWeight','bold')
set(0,'defaultLineLineWidth',2)

% The radar is at the Navy base at Point Loma at a location with a north/south
% orientation of the coast. Th Tybee starts in the morning from the US
% Coast Guard station near San Diego airport, whihc is inside San Diego harbor.
% Thus, the Tybee should start out heading west and then turn north to the
% radar site.

BaseLat      =   32.70326358;                     % Latitude of the radar
BaseLon      = -117.25496430;                     % Longitude of the radar

%------------------------------------------------------------------------------

% The 'Mission Bay Offshore' buoy Location in March 1998 according to the
% CDIP website, the archive for the Scripps Institution of Oceanography
BuoyLat      =   32+44.9/60;
BuoyLon      = -117-22.1/60;

dLatBuoy     = BuoyLat-BaseLat;
dLonBuoy     = BuoyLon-BaseLon;

Re           = 6378000;                           % Earth radius (m)
northBuoy    = dLatBuoy*Re*pi/180;                % North component (m)
eastBuoy     = dLonBuoy*Re*pi*cosd(BaseLat)/180;  % East component (m)

fprintf(1,'\n Mission Bay Offshore Buoy Location from Radar:     %6.2f km West %6.2f km North\n',-eastBuoy/1000,northBuoy/1000);

%------------------------------------------------------------------------------

% % Do a sanity check for the first 30 minutes of the file to make sure the
% % signs are correct
% %
% % The original data had the Tybee heading east into San Diego and then turning
% % south towards Mexico. But Jay Trischman figured out that this must have been
% % due to a swapping of the ship Lat/Lon and the reference for the differential GPS.
% % Thus, we correct the signs of the differences.
% 
% tybee = load('SC15GPS_First1800.txt','-ascii');
% dLat  = -(tybee(:,1)-BaseLat);
% dLon  = -(tybee(:,2)-BaseLon);
% 
% % Plot the two time series for the lat and lon offsets and then the 2-D
% % plot, skipping the first 100 points
% 
% figure(101);clf
% plot(dLat(101:end))
% title('Diff. Latitude')
% 
% figure(102);clf
% plot(dLon(101:end))
% title('Diff. Longitude')
% 
% figure(103);clf
% plot(dLon(101:end),dLat(101:end))
% title('Lat/Lon')

%------------------------------------------------------------------------------

% Two hours (7200 1-second samples) of GPS data from samples 1880 to 9079
% (7200 samples) of the file SC15GPS.CSV from NRL

% Read the 7200 samples
tybee        = load('SC15GPS_1880_7200.txt','-ascii')';

ngps         = 3000;            % Number of GPS samples to be used here

tybee        = tybee(:,1:ngps); % Shorten to the desired samples

% First sample for the 8 Tybee cases, according to the NRL spreadsheet
FirstSamp    = [1, 368, 705, 1068, 1423, 1787, 2150, 2513];

% Add 60 seconds to time align R0 and A0. These initial values were chosen
% by the analyst at NRL to avoid the time when the Tybee was changing its
% heading to the new required value for this section of the octagon path.
FirstSamp    = FirstSamp+60;

%------------------------------------------------------------------------------

% ISAR guidance information from the NRL short spreadsheet used here only
% for comparison with the GPS results

% Initial range, Tybee
R0           = [ 4253, 5120, 6822, 8380, 9065, 8649, 7136, 5126];

% Initial radial velocity, Tybee
Vr0          = [ 0.7, 4.8, 3.5,  2.9, -0.2, -4.4, -7.6, -5.0];

% Aspect angles, Tybee
A0           = [ -78.27, -16.25, 36.10, 70.01, 101.90, 130.96, 173.63, -145.62];

%------------------------------------------------------------------------------

dLat         = -(tybee(1,:)-BaseLat);
dLon         = -(tybee(2,:)-BaseLon);

rFile        = tybee(6,:)*1852; % Convert nautical miles to meters
HdgFile      = tybee(4,:);      % GPS Heading from file
AzGPS        = tybee(7,:);      % Antenna azimuth

% Plot the Lat/Lon path relative to the radar position. The path of the
% Tybee is counterclockwise.
figure(104);clf
plot(dLon,dLat,'b')
grid on
title ('Tybee GPS Path 09:26:17-10:16:17 March 6, 1998')
xlabel('Longitude Diff. (deg)')
ylabel('Latitude Diff. (deg)')

% Time steps and times for the GPS records relative to the initial time
dtgps        = 1;                                 % Sampling time for GPS data (1 second)
tgps         = dtgps*(0:ngps-1);                  % Time offsets for GPS data

% Convert the lat/lon coordinates to meters
Re           = 6378000;                           % Earth radius (m)
north0       = dLat*Re*pi/180;                    % North component (m)
east0        = dLon*Re*pi*cosd(BaseLat)/180;      % East component (m)

% Range offset from the Tybee to the buoy
RBuoykm      = sqrt((northBuoy-north0).^2+(eastBuoy-east0).^2)/1000;
fprintf(1,'\n Mean Mission Bay Offshore Buoy Range from Tybee:   %6.2f km \n',mean(RBuoykm));

% Smooth the raw position data slightly
north0       = c1smooth(north0,1);
east0        = c1smooth(east0,1);

% Ship velocity in north/east components
ndot         = diff(north0)/dtgps;
edot         = diff(east0)/dtgps;

% Smooth the raw velocity data slightly
ndot         = c1smooth(ndot,1);
edot         = c1smooth(edot,1);

% Center north and east to have the same times as ndot and edot
north        = north0(1:end-1);
east         = east0(1:end-1);
for k = 1:ngps-1
   north(k) = (north0(k)+north0(k+1))/2;
   east(k)  = (east0(k)+east0(k+1))/2;
end

% Calculate range versus time
RLL          = sqrt(north.^2+east.^2);

% Define the times for velocity estimates and the ship speed
tV           = zeros(1,ngps-1);
Spd          = zeros(1,ngps-1);
for k = 1:ngps-1
   tV(k)    = (tgps(k)+tgps(k+1))/2;
   Spd(k)   = sqrt(ndot(k)^2+edot(k)^2);
end

% Calculate the in-range and cross-range velocities and the in-track and
% cross-track velocities plus the aspect angle

[AspLL,BearingLL,HeadingLL,VIR,VCR] = GetAspect(edot,ndot,east,north);

save('Tybee_1880_3000.mat','FirstSamp','tV','edot','ndot','RLL','VIR','VCR','AspLL','Spd','HeadingLL','BearingLL');

% At this point all values for the test program for the ISAR AutoTrack algorithm are done.
% That program just reads the above mat-file.

% The following code contains the details of the evaluation the accuracy of data in the NRL 
% spreadsheets.

%---------------------------------------------------------------------------------------------
%
% Check the range versus the NRL file
%
% Range offset from the GPS data versus the NRL file
rFile        = rFile(1:end-1);

figure(105);clf
plot(tV,rFile,'b')
hold on
plot(tV,RLL,'r')
plot(tV,rFile-RLL,'g')
scatter(FirstSamp,R0,'k','filled')
xlabel('GPS Time (sec)')
ylabel('Range (m)')
title ('Range from the Lat/Lon in the GPS data vs the file value')
legend('GPS-calc','GPS-file','Diff','R0')
grid on; axis tight

figure(106);clf
plot(tV,rFile-RLL,'g')
xlabel('GPS Time (sec)')
ylabel('Range Diff. (m)')
title ('Diff. GPS-calc, GPS-file in meters')
grid on; axis tight

%---------------------------------------------------------------------------------------------
%
% Check the aspect angle calculated from the GPS Lat/Lon data and the NRL file
%

% Choose to make heading 0-360 to match the NRL GPS spreadsheet
HeadingFile2 = HeadingLL;
for k = 1:numel(HeadingFile2)
   if HeadingFile2(k) < 0; HeadingFile2(k) = HeadingFile2(k)+360; end
end

% Difference between this calculation from the GPS data and the values
% listed in the GPS spreadsheet
HeadingDif   = HeadingFile2-HdgFile(1:end-1);
for k = 1:numel(HeadingDif)
   if HeadingDif(k) < -180; HeadingDif(k) = HeadingDif(k)+360; end
   if HeadingDif(k) >  180; HeadingDif(k) = HeadingDif(k)-360; end
end

figure(107);clf
plot(tV,HeadingFile2,'b')
hold on
plot(tV,HdgFile(1:end-1),'r')
plot(tV,HeadingDif,'g')
xlabel('GPS Time (sec)')
ylabel('Heading (deg)')
title ('Heading from the Lat/Lon in the GPS data vs the file value')
legend ('Hdg-calc','Hdg-file','Diff')
grid on; axis tight

%---------------------------------------------------------------------------------------------
%
% Evaluate the system aspect, taking into account that the file value is defined so
% that zero aspect is a bow view
AspFile      = tybee(5,:)-180;
for k = 1:numel(AspFile)
   if AspFile(k) >  180; AspFile(k) = AspFile(k)-360; end
   if AspFile(k) < -180; AspFile(k) = AspFile(k)+360; end
end

AspectDif    = AspLL-AspFile(1:end-1);
for k = 1:numel(AspectDif)
   if AspectDif(k) >  180; AspectDif(k) = AspectDif(k)-360; end
   if AspectDif(k) < -180; AspectDif(k) = AspectDif(k)+360; end
end

figure(108);clf
plot(tV,AspLL,'b')
hold on
plot(tV,AspFile(1:end-1),'r')
plot(tV,AspectDif,'g')
scatter(FirstSamp,A0,'k','filled')
xlabel('GPS Time (sec)')
ylabel('Aspect (deg)')
title ('Aspect from the Lat/Lon in the GPS data vs the file value')
legend('Asp-calc','Asp-file','Diff','A0')
grid on; axis tight

%---------------------------------------------------------------------------------------------

figure(109);clf
plot(tV,AspLL,'b')
xlabel('GPS Time (sec)')
ylabel('Aspect (deg)')
legend('AspLL')
grid on; axis tight

figure(110);clf
plot(tV,VIR,'m')
hold on
scatter(FirstSamp,Vr0,'k','filled')
xlabel('GPS Time (sec)')
ylabel('In-range Velocity (m/s)')
title ('In-range velocity from the GPS data vs the NRL guidance')
legend('GPS','Vr0')
grid on; axis tight

figure(111);clf
plot(tV,VCR,'b')
hold on
xlabel('GPS Time (sec)')
ylabel('Cross-range Velocity (m/s)')
title ('Cross-range velocity from the GPS data')
legend('XrV from GPS')
grid on; axis tight; hold off

BearFirst = zeros(1,8);
for k = 1:8
   BearFirst(k) = BearingLL(FirstSamp(k));
end
figure(112);clf
plot(tV,BearingLL,'b')
hold on
plot(tgps,AzGPS-360,'r')
scatter(FirstSamp,BearFirst,'k','filled')
% plot(tgps,BearingLL-AzGPS+360,'g')
xlabel('GPS Time (sec)')
ylabel('Bearing/Azimuth (deg)')
% legend('BearingLL','AzGPS','Diff','FirstSamp')
legend('BearingLL','AzGPS','FirstSamp')
grid on; axis tight; hold off

fprintf(1,'\n RMS Bearing-Azimuth (deg):                       %8.4f\n',rms(BearingLL-AzGPS(1:end-1)+360));

figure(113);clf
plot(tV,east,'b')
hold on
plot(tV,north,'r')
xlabel('Time (s)')
ylabel('East/North (m)')
legend('East','North')
grid on; axis tight; hold off

figure(114);clf
plot(tV,edot,'b')
hold on
plot(tV,ndot,'r')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
legend('Edot','Ndot')
grid on; axis tight; hold off

fprintf(1,'\n Mean and standard deviations for East and North velocities: %9.3f %9.3f %9.3f %9.3f\n',mean(edot),mean(ndot),std(edot),std(ndot));

figure(115);clf
plot(tV,AspLL,'b')
hold on
plot(tV,atan2d(VCR,VIR),'r')
xlabel('Time (s)')
ylabel('Aspect (deg)')
grid on; axis tight; hold off

%---------------------------------------------------------------------------------------------

Adiff1 = 0*A0;
Adiff2 = 0*A0;

for m = 1:8

   Adiff1(m) = A0(m)-AspLL(FirstSamp(m));
   Adiff2(m) = A0(m)-AspFile(FirstSamp(m));
   
   fprintf(1,['\n m, FirstSamp, R0, RLL, R, A0, AspLL, AspFile: '            ...
              ' %8i %8i %8.0f %8.0f %8.0f %8.1f %8.1f %8.1f %8.1f %8.1f\n'], ...
              m,FirstSamp(m),                                                ...
              R0(m),RLL(FirstSamp(m)),rFile(FirstSamp(m)),                   ...
              A0(m),AspLL(FirstSamp(m)),AspFile(FirstSamp(m)),Adiff1(m),Adiff2(m));

end
