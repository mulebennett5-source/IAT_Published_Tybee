
% Sample execution of the tracker accuracy tool Sept. 21, 2020 with
% program output for the following defaults

% icase             = 1 - Circular Scan Mode
% Range             = 60, 90, 120, 150, 180, 210 km
% Range resolution  = 5 m
% PulseSNRdb        = 14 dB for each scatterer (but with 10 scatterers)
% Ship aspect angle = 30 degrees
% Sweep             = 1/3 - Width of scan, fraction of a rotation

% The mean values below are averages over aspect angles from 5 to 85
% degrees - the input aspect is only used for Plot 1

TrackAccuracyTest(1,60000,5,14,30,1/3);
% Mean heading accuracy (deg):        0.50912
% Mean length accuracy (percent):     0.88874


TrackAccuracyTest(1,90000,5,14,30,1/3);
% Mean heading accuracy (deg):        0.76379
% Mean length accuracy (percent):     1.33378


TrackAccuracyTest(1,120000,5,14,30,1/3);
% Mean heading accuracy (deg):        1.01863
% Mean length accuracy (percent):     1.77966


TrackAccuracyTest(1,150000,5,14,30,1/3);
% Mean heading accuracy (deg):        1.27367
% Mean length accuracy (percent):     2.22667


TrackAccuracyTest(1,180000,5,14,30,1/3);
% Mean heading accuracy (deg):        1.52897
% Mean length accuracy (percent):     2.67511


TrackAccuracyTest(1,210000,5,14,30,1/3);
% Mean heading accuracy (deg):        1.78460
% Mean length accuracy (percent):     3.12528
