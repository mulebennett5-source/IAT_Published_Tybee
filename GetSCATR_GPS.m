function [RdotBar,RLLBar,SpeedBar,VgpsR,VgpsC,AspBarLL,HdgBarLL,tgpsV,moc32,tmoc,dtmoc] = ...
          GetSCATR_GPS(filename,k,kFr,isarprf,FirstPulse,NPulses,toffsetGPS,FrTime,plotMoco)

%------------------------------------------------------------------------------

% Load the adaptive motion compensation data from the ISAR processor
moc           = load(filename,'-ascii')';
tmoc          = moc(1,FirstPulse:FirstPulse-1+NPulses);
dtmoc         = mean(diff(tmoc));

tmoc          = tmoc-tmoc(1);

moc32         = moc(3,FirstPulse:FirstPulse-1+NPulses);   % Radial velocity during ISAR
RMoco         = moc(2,FirstPulse:FirstPulse-1+NPulses);   % Range offset for motion compensation

RdotBar       = mean(moc32);

%---------------------------------------------------------------------------

fprintf(1,'\n Parameters from the GPS routine that apply to all methods:\n')

% Save statement from TybeeGPS.m:
% save('Tybee_1880_3000.mat','FirstSamp','tV','edot','ndot','RLL','VIR','VCR','AspLL','Spd','HeadingLL','BearingLL');

if filename(1:4) == 'sc01'; load('Point_Stuart.mat'); end
if filename(1:4) == 'sc15'; load('Tybee_1880_3000.mat'); end

SampOffset    = FirstSamp(1)-1;

% The first sample depends on the frame number
FSamp         = FirstSamp-SampOffset+toffsetGPS+(kFr-1)*FrTime;

% All 4 velocities plus Rbar and SpeedBar have the same length, FrTime

FrameSamples  = FSamp(k):FSamp(k)+FrTime-1;

RLLBar        = mean(RLL(FrameSamples));

SpeedBar      = mean(Spd(FrameSamples));

VgpsR         = VIR(FrameSamples);        % In-range velocity

VgpsC         = VCR(FrameSamples);        % Cross-range velocity

AspLL         = AspLL(FrameSamples);

BearLL        = BearingLL(FrameSamples);

HdgLL         = HeadingLL(FrameSamples);

edot          = edot(FrameSamples);       % East velocity

ndot          = ndot(FrameSamples);       % North velocity

fprintf(1,'\n Mean and standard deviations for East and North velocities: %9.3f %9.3f %9.3f %9.3f\n',mean(edot),mean(ndot),std(edot),std(ndot));

dtgps         = mean(diff(tV));

tgpsV         = 0.5*dtgps+tV(FrameSamples)-tV(FSamp(k));

% Time offset centered at zero
tgpsVp        = tgpsV-mean(tgpsV);

%---------------------------------------------------------------------------

% NRL guidance for the time mean of the aspect angle - do not use this
% AspBar        = A0(k);

% Calculate the time mean of the aspect angle from the GPS data
Xa            = mean(sind(AspLL));
Ya            = mean(cosd(AspLL));
AspBarLL      = atan2d(Xa,Ya);

% Calculate the rotation rate, ignoring possible 360 degree wraps, which do
% not exist for this data
[La,Ma,Aa]    = LinReg(AspLL,tgpsVp);

fprintf(1,'\n Mean Aspect (deg) and rate from GPS data (deg/s):      %9.4f %9.4f\n',Ma,Aa);

% Calculate the time mean of the GPS aspect angles
Xb            = mean(sind(BearLL));
Yb            = mean(cosd(BearLL));
BearBarLL     = atan2d(Xb,Yb);

% Calculate the rotation rate, ignoring possible 360 degree wraps, which do
% not exist for this data
[Lb,Mb,Ab]    = LinReg(BearLL,tgpsVp);

fprintf(1,'\n Mean bearing (deg) rate from GPS data (deg/s):         %9.4f %9.4f\n',BearBarLL,Ab);

% Calculate the time mean of the GPS heading
Xh            = mean(sind(HdgLL));
Yh            = mean(cosd(HdgLL));
HdgBarLL      = atan2d(Xh,Yh);

% Calculate the in-track and cross-track velocities, based on the mean
% heading

%--------------------------------------------------------------------------
%
% Evaluate variations in the ship heading
VIT           =  cosd(HdgBarLL)*ndot+sind(HdgBarLL)*edot;
VCT           = -cosd(HdgBarLL)*edot+sind(HdgBarLL)*ndot;

XIT           = dtgps*cumsum(VIT-mean(VIT));
XCT           = dtgps*cumsum(VCT-mean(VCT));

XIT           = XIT-XIT(floor(numel(XIT)/2)+1);
XCT           = XCT-XCT(floor(numel(XCT)/2)+1);

fprintf(1,'\n Std of In-track and Cross-track Displacements:         %9.4f %9.4f\n',std(XIT),std(XCT));

% Calculate the rotation rate, ignoring possible 360 degree wraps, which do
% not exist for this data
[Lh,Mh,Ah]    = LinReg(HdgLL,tgpsVp);

fprintf(1,'\n Mean heading (deg), Speed (m/s):                       %9.4f %9.4f\n',HdgBarLL,SpeedBar);
fprintf(1,'\n Mean heading (deg) rate from GPS data (deg/s):         %9.4f %9.4f\n',Mh,Ah);

eastBar      = Xb*RLLBar;
northBar     = Yb*RLLBar;

% Aspect rate based on the mean motion of the ship

% Mean east and north velocities
edotBar      = Xh*SpeedBar;
ndotBar      = Yh*SpeedBar;

% Time-dependent east and north positions
eastBar_t    = eastBar +edotBar*tgpsVp;
northBar_t   = northBar+ndotBar*tgpsVp;

% Time-dependent bearing and aspect
BearBar_t    = atan2d(eastBar_t,northBar_t);

%==========================================================================================================

% The calculations below use the mean heading and thus miss the normal time
% change of it. The true aspect angle and its time change are computed in
% the first section above.
 
AspBar_t     = BearBar_t-HdgBarLL;

for l = 1:numel(AspBar_t)
   if AspBar_t(l) >  180; AspBar_t(l) = AspBar_t(l)-360; end
   if AspBar_t(l) < -180; AspBar_t(l) = AspBar_t(l)+360; end
end

% Calculate the rotation rate, ignoring possible 360 degree wraps, which do
% not exist for this data
[Li,Mi,Ai]   = LinReg(AspBar_t,tgpsVp);

AaIdeal      = -(180/pi)*mean(VgpsC)/RLLBar;

fprintf(1,'\n Ideal Aspect (deg) and rates from GPS data (deg/s):    %9.4f %9.4f %10.4f\n',Mi,Ai,AaIdeal);

%==========================================================================================================

if plotMoco && (FrTime == 240) % Do only if the plotMoco flag is set and if we are processing the
                               % full 4-minute data set

   % These lines are only for advanced diagnostics - they have no impact on the results of the IAT
   
   figure(200+k);clf
   plot(tgpsV,HdgLL,'b')
   hold on
   plot(tgpsV,HdgBarLL-atan2d(VCT,VIT),'r')
   title ([filename(1:end-4) ': Heading versus Track Velocities'])
   xlabel('time (s)')
   ylabel('Heading (deg)')
   axis tight; grid on; hold off

   figure(208+k);clf
   plot(tgpsV,VIT,'b')
   hold on
   plot(tgpsV,VCT,'r')
   plot(tgpsV,VgpsR,'g')
   plot(tgpsV,VgpsC,'m')
   title ([filename(1:end-4) ': In/Cross Range vs In/Cross Track'])
   legend('In-track','Cross-track','In-range','Cross-range')
   xlabel('time (s)')
   ylabel('Velocity (m/s)')
   axis tight; grid on; hold off

   figure(216+k);clf
   plot(tgpsV,XIT,'b')
   hold on
   plot(tgpsV,XCT,'r')
   title ([filename(1:end-4) ': In-track vs Cross-track Displacement'])
   legend('In-track','Cross-track')
   xlabel('time (s)')
   ylabel('Displacement (m)')
   axis tight; grid on; hold off

   % Asp and atan2d(VgpsC,VgpsR) are equal so the green curve overlays the blue curve
   figure(224+k);clf
   plot(tgpsV,AspLL,'b')
   hold on
   plot(tgpsV,AspBar_t,'r')
   plot(tgpsV,atan2d(VgpsC,VgpsR),'g')
   legend('Asp','AspBar_t','AspAtan2')
   grid on; hold off
   
   figure(232+k);clf
   plot(tgpsV,HdgLL-HdgBarLL,'b')
   hold on
   plot(tgpsV,BearLL-BearBarLL,'r')
   plot(tgpsV,AspLL-AspBarLL,'g')
   title ([filename(1:end-4) ': GPS Heading, Bearing and Aspect offsets'])
   legend('Heading','Bearing','Aspect')
   xlabel('Time (s)')
   ylabel('Aspect (deg)')
   grid on; hold off; axis tight
   saveas(gcf,['Hdg_vs_Bear_' filename(1:end-4) ],'fig')

   figure(240+k);clf
   plot(tgpsV,edot,'b')
   hold on
   plot(tgpsV,ndot,'r')
   plot(tmoc,moc32,'g')
   xlabel('Time (s)')
   ylabel('Velocity (m/s)')
   title ([filename(1:end-4) ': GPS and ISAR Velocities'])
   legend('Edot','Ndot','ISAR')
   grid on; axis tight; hold off
   saveas(gcf,['GPS_Velociites_' filename(1:end-4)],'fig')

   Nf                 = 256;
   df                 = 1/(dtgps*Nf);

   [rspecf,f,ruse]    = tspect(Nf,VgpsR,tgpsV,dtgps,0);
   [cspecf,f,cuse]    = tspect(Nf,VgpsC,tgpsV,dtgps,0);

   nfmoc              = 48000;
   
   %==========================================================================================================

   down               = isarprf;         % Ratio of ISAR and GPS sampling rates

   % Downsample by pre-summing by 'down'
   mdown              = zeros(1,nfmoc/down);
   tdown              = zeros(1,nfmoc/down);
   for l = 1:nfmoc/down
      mdown(l) = mean(moc32(1+down*(l-1):down*l));
      tdown(l) = mean(tmoc(1+down*(l-1):down*l));
   end
   nfmocdown  = nfmoc/down;
   dtdown     = dtmoc*down;
   [mocspecf1,fm1,m1] = tspect(nfmocdown,mdown,tdown,dtdown,0);

   % Downsample by sampling by 'down'
   mdown2             = zeros(1,nfmoc/down);
   tdown2             = zeros(1,nfmoc/down);
   for l = 1:nfmoc/down
      mdown2(l) = moc32(nfmocdown/2+1+down*(l-1));
      tdown2(l) = tmoc(nfmocdown/2+1+down*(l-1));
   end
   [mocspecf2,fm2,m2] = tspect(nfmocdown,mdown2,tdown2,dtdown,0);

   %==========================================================================================================

   [mocspecf,fmoc,m0] = tspect(nfmoc,moc32,tmoc,dtmoc,0);

   fprintf(1,'\n Means of moc32 and mdown:   %12.6f %12.6f\n',mean(moc32),mean(mdown));
   
   fprintf(1,'\n Input and Output Variances: %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n', ...
                 var(moc32),var(m0),mean(mocspecf.^2),mean(mocspecf1.^2),var(ruse),mean(rspecf.^2));

   figure(248+k);clf
   plot(fmoc(2+nfmoc/2:end),20*log10(mocspecf(2+nfmoc/2:end)),'b')
   hold on
   plot(f(2+Nf/2:end),20*log10(rspecf(2+Nf/2:end)),'r')
   legend('ISAR','GPS')
   title ([filename(1:end-4) ': GPS vs ISAR In-range Velocity Spectra'])
   xlabel('Frequency (Hz)')
   ylabel('Power (dB)')
   axis tight; grid on; hold off
   xlim([0 0.5])
   saveas(gcf,['GPSvsISAR_VelocitySpectra_' filename(1:end-4)],'fig')

   figure(256+k);clf
   plot(fmoc(2+nfmoc/2:end),20*log10(mocspecf(2+nfmoc/2:end)),'b')
   hold on
   plot(f(2+Nf/2:end),20*log10(rspecf(2+Nf/2:end)),'r')
   plot(f(2+Nf/2:end),20*log10(cspecf(2+Nf/2:end)),'g')
   plot(fm1(2+nfmocdown/2:end),20*log10(mocspecf1(2+nfmocdown/2:end)),'m')
   plot(fm2(2+nfmocdown/2:end),20*log10(mocspecf2(2+nfmocdown/2:end)),'k')
   legend('ISAR','GPS In-Range','GPS Cross-Range','ISAR200Presum','ISAR200Resamp')
   title ([filename(1:end-4) ': GPS vs ISAR In-range Velocity Spectra'])
   xlabel('Frequency (Hz)')
   ylabel('Power (dB)')
   axis tight; grid on; hold off
   xlim([0 0.5])
   saveas(gcf,['GPSvsISAR_VelocitySpectra_All_' filename(1:end-4)],'fig')

   figure(272+k);clf
   plot(tmoc,m0,'b')
   hold on
   plot(tgpsV,ruse,'r')
   legend('ISAR','GPS')
   title ([filename(1:end-4) ': GPS vs ISAR Detrended Velocity'])
   axis tight; grid on; hold off
   saveas(gcf,['GPSvsISAR_DetrendedVelocity_' filename(1:end-4)],'fig')

end

fprintf(1,'\n Mean and Std. of heading from GPS Data:                %9.4f %9.4f\n',HdgBarLL,std(HdgLL));

%---------------------------------------------------------------------------
