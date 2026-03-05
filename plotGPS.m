function [Rgps,VgpsR,VgpsC,tgpsV,AspGPS,Speed,tgps,HeadingLL,BearingLL] = plotGPS(filename)

BaseLat      = 32.70326358;                       % Latitude of the radar
BaseLon      = -117.2549643;                      % Longitude of the radar; but not used here

LL           = load([filename '.txt'],'-ascii');  % Selected GPS data, including the North and
                                                  % East offsets from the radar (degrees)

% Time steps and times for the GPS records relative to the initial time
dtgps        = 1;                                 % Sampling time for GPS data (1 second)
ngps         = size(LL,1);                        % Number of GPS pulses for this file
tgps         = dtgps*(0:ngps-1);                  % Time offsets for GPS data

% Convert the ship coordinates to meters taking into account the fact that
% the GPS data appear to be flipped in sign due to an apparent software glitch
Re           = 6378000;                           % Earth radius (m)
north        = -LL(:,7)'*Re*pi/180;               % North component (m)
east         = -LL(:,8)'*Re*pi*cosd(BaseLat)/180; % East compomnent (m)

% Smooth the raw data slightly
north        = c1smooth(north,2);
east         = c1smooth(east,2);

% Range offset and range velocity for raw data
Rgps         = sqrt(north.^2+east.^2);
VgpsR0       = diff(Rgps)/dtgps;

% Ship velocity in north/east components
ndot         = diff(north)/dtgps;
edot         = diff(east)/dtgps;

% Define the times for velocity estimates and the ship speed
tgpsV        = zeros(1,ngps-1);
Speed        = zeros(1,ngps-1);
for k = 1:ngps-1
   tgpsV(k) = (tgps(k)+tgps(k+1))/2;
   Speed(k) = sqrt(ndot(k)^2+edot(k)^2);
end

% Mean and difference for the north and east positions
MeanNorth    = mean(north);
MeanEast     = mean(east);

DifNorth     = mean(diff(north));
DifEast      = mean(diff(east));

% Check the aspect angle calculated from the GPS Lat/Lon data and the NRL file

% Azimuth and heading calculated from the GPS Lat/Lon info
BearingLL    = atan2d(MeanEast,MeanNorth);
HeadingLL    = atan2d(DifEast,DifNorth);

% In the file the first element is the estimated heading
xbar         = mean(sind(LL(:,1)));
ybar         = mean(cosd(LL(:,1)));
HdgFile      = atan2d(xbar,ybar);

HeadingDif   = HeadingLL-HdgFile;
if HeadingDif >  180; HeadingDif = HeadingDif-360; end
if HeadingDif < -180; HeadingDif = HeadingDif+360; end

% Check the aspect angle calculated from the GPS Lat/Lon data and the NRL file
AspGPS       = BearingLL-HeadingLL;
if AspGPS >  180; AspGPS = AspGPS-360; end
if AspGPS < -180; AspGPS = AspGPS+360; end

% Evaluate the system aspect, taking into account that the file value is defined so
% that zero aspect is a bow view
xbar         = mean(sind(LL(:,2)));
ybar         = mean(cosd(LL(:,2)));
AspFile      = atan2d(-xbar,-ybar);  % 180 degree shift by flipping the two signs
if AspFile >  180; AspFile = AspFile-360; end
if AspFile < -180; AspFile = AspFile+360; end

AspectDif    = AspGPS-AspFile;
if AspectDif >  180; AspectDif = AspectDif-360; end
if AspectDif < -180; AspectDif = AspectDif+360; end

fprintf(1,'\n Calculated Heading, File Heading, Diff:               %8.3f %8.3f %8.3f', ...
              HeadingLL,HdgFile,HeadingDif)

fprintf(1,'\n Calculated Aspect, File Aspect, Diff:                 %8.3f %8.3f %8.3f', ...
              AspGPS,AspFile,AspectDif)

fprintf(1,'\n Bearing from GPS data:                                %8.3f',BearingLL);

%---------------------------------------------------------------------------------------------

X            = zeros(1,ngps-1);
Y            = zeros(1,ngps-1);
for k = 1:ngps-1
   X(k) = (east(k)+east(k+1))/(Rgps(k)+Rgps(k+1));
   Y(k) = (north(k)+north(k+1))/(Rgps(k)+Rgps(k+1));
end

VgpsR        =  X.*edot+Y.*ndot;
VgpsC        = -Y.*edot+X.*ndot;

AspGPS2      = atan2d(mean(VgpsC),mean(VgpsR));

AspectDif    = AspGPS-AspGPS2;
if AspectDif >  180; AspectDif = AspectDif-360; end
if AspectDif < -180; AspectDif = AspectDif+360; end

fprintf(1,'\n Calculated Aspect, Mean Components Aspect, Diff:      %8.3f %8.3f %8.3f', ...
              AspGPS,AspGPS2,AspectDif)

compSpeed   = sqrt(mean(VgpsR)^2+mean(VgpsC)^2);

fprintf(1,'\n Mean(Speed), Speed of Mean Components, Diff.:         %8.3f %8.3f %8.3f',mean(Speed),compSpeed,mean(Speed)-compSpeed);

fprintf(1,'\n RMS for Rdot measures: %12.10f\n',sqrt(mean((VgpsR-VgpsR0).^2)));
