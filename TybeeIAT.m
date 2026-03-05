function [] = TybeeIAT(icase)

if nargin == 0; icase = 5; end

% Parameters for the extrapolated map-drift equation
Rmult      = 10;    % Multiplier for range
xrpv       = 50;    % Synthetic cross-range platform velocity (m/s)
Nmult      = 1;     % Multiplier for noise
Amult      = 1.0;   % Multiplier for acceleration noise

% % Parameters for the extrapolated map-drift equation
% Rmult      = 1;     % Multiplier for range
% xrpv       = 0;     % Synthetic cross-range platform velocity (m/s)
% Nmult      = 1;     % Multiplier for noise
% Amult      = 0;     % Multiplier for acceleration noise

% Frperleg: Frames per octagon leg
% FrTime:   Frame time in minutes
% plotMoco: Option to plot many motion compensation details

if     icase == 1
   Frperleg = 16;    
   FrTime   = 15;
   ofile    = '15T\';
elseif icase == 2
   Frperleg = 8;    
   FrTime   = 30;
   ofile    = '30T\';
elseif icase == 3
   Frperleg = 4;
   FrTime   = 60;
   ofile    = '60T\';
elseif icase == 4
   Frperleg = 2;
   FrTime   = 120;
   ofile    = '120T\';
elseif icase == 5
   Frperleg = 1;
   FrTime   = 240;
   ofile    = '240T\';
end

set(0,'DefaultAxesTitleFontWeight','bold')
set(0,'DefaultAxesFontWeight','bold')
set(0,'defaultLineLineWidth',2)

%------------------------------------------------------------------------------
%
%      The ISAR AutoTrack (IAT) algorithm uses the Map-Drift equation to
%      estimate the cross-range velocity for a moving radar with no
%      tracker:
%
%                            LHS = RHS
%
%      Mocomp_Acceleration*Range = V^2-2*V*u+u^2
%
%      Where                   V = Crossrange platform velocity (m/s)
%
%                              u = Crossrange target velocity (m/s)
%
%      The above equation is used for the Extrapolated solution to the
%      problem
% 
%      For a stationary radar this reduces to:
%
%      Mocomp_Acceleration*Range = u^2
%
%      This equation is used for the Model solution - the result when all
%      parameters are estimated from the GPS data
%
%------------------------------------------------------------------------------

plotMoco   = 1;                % Control for mocomp and advanced diagnostic plots

numfiles   = 8;                % SCATR Data: 8 each cases for USCG Tybee

toffsetI3  = 2;                % Time offset (s) to first frame for ISAR and 3D data

toffsetGPS = 60;               % Time offset (s) to first frame for GPS data

isarprf    = 200;              % Pulse repetition frequency for ISAR

pulseskip  = 32;               % Skip the acquisition mode pulses used by the ISAR processor

NPulses    = isarprf*FrTime;

mocfiles   = dir('sc15*.moc'); % Motion compensation data from the ISAR processor

%------------------------------------------------------------------------------

% Input guidance data from the short NRL spreadsheet - not used here
% because the data derived from the GPS data works as well

% Initial range, Tybee Cases
R0         = [ 4253, 5120, 6822, 8380, 9065, 8649, 7136, 5126];

% Initial radial velocity, Tybee Cases
Vr0        = [ 0.7, 4.8, 3.5, 2.9, -0.2, -4.4, -7.6, -5.0];

% Aspect angles, Tybee Cases
A0         = [ -78.27, -16.25, 36.10, 70.01, 101.90, 130.96, 173.63, -145.62];
   
%------------------------------------------------------------------------------

% Right and left hand side of the map-drift equation for a stationary
% platform with no tracker

RHSM       = zeros(Frperleg,numfiles);         % Model version of the IAT equation - right
LHSM       = zeros(Frperleg,numfiles);         % and left hand sides
AccM       = zeros(Frperleg,numfiles);
AccV       = zeros(Frperleg,numfiles);
dVM        = zeros(Frperleg,numfiles);
xrvM       = zeros(1,numfiles);
AspM       = zeros(1,numfiles);
AccG       = zeros(Frperleg,numfiles);

RHSE       = zeros(Frperleg,numfiles);         % Extrapolated version of the IAT equation
LHSE       = zeros(Frperleg,numfiles);  
dVE        = zeros(Frperleg,numfiles);         % Time-integrated acceleration (m/s)
AccE       = zeros(Frperleg,numfiles);         % Extrapolated acceleration
AspDif     = zeros(Frperleg,numfiles);         % Aspect angle difference
MocoAcc0   = zeros(Frperleg,numfiles);         % Acceleration before data preparation
MocoAcc    = zeros(Frperleg,numfiles);         % Acceleration after data preparation

% Mean range and crossrange velocities and aspect angles for the IAT output and GPS input
rvIAT      = zeros(Frperleg,numfiles);
rvGPS      = zeros(Frperleg,numfiles);
xrvIAT     = zeros(Frperleg,numfiles);
xrvGPS     = zeros(Frperleg,numfiles);
AspIAT     = zeros(Frperleg,numfiles);
AspGPS     = zeros(Frperleg,numfiles);

% Other variables for plotting
RdotBar    = zeros(Frperleg,numfiles);
AspBar     = zeros(Frperleg,numfiles);
HdgBar     = zeros(Frperleg,numfiles);
Rbar       = zeros(Frperleg,numfiles);
dRmaxA     = zeros(Frperleg,numfiles);
dRmaxA_E   = zeros(Frperleg,numfiles);

RHSD       = zeros(1,numfiles);                % Demo version of the IAT equation
LHSD       = zeros(1,numfiles);
dVD        = zeros(1,numfiles);
AccD       = zeros(1,numfiles);
xrvIATD    = zeros(1,numfiles);
AspIATD    = zeros(Frperleg,numfiles);
AspDifD    = zeros(Frperleg,numfiles);

RHSQ       = zeros(Frperleg,numfiles);         % Data version of the IAT equation in the IAT paper
LHSQ       = zeros(Frperleg,numfiles);
dVQ        = zeros(Frperleg,numfiles);

Phi3DAcc   = zeros(Frperleg,numfiles);
VoverR     = zeros(Frperleg,numfiles);

for k = 1:numfiles       % Each file is a leg of the octagonal path
   
   for kFr = 1:Frperleg  % Frames per leg, equal 1 when using the whole 4 minutes
      
      fprintf(1,['\n\n ' mocfiles(k).name '\n']);
      
      filename          = mocfiles(k).name;
   
      % Get the parameters from the GPS routine that apply to all methods
   
      FirstPulse        = 1+pulseskip+toffsetI3*isarprf+(kFr-1)*FrTime*isarprf;

      [RdotBar(kFr,k),Rbar(kFr,k),SpeedBar0,VgpsR,VgpsC,AspBar(kFr,k),HdgBar(kFr,k),tgpsV,moc32,tmoc,dtmoc] = ...
                          GetSCATR_GPS(filename,k,kFr,isarprf,FirstPulse,NPulses,toffsetGPS,FrTime,plotMoco);
   
      % Get 3-D Diagnostics - not related to the IAT
      if plotMoco > 1 && (icase == 5)
         [time3D,dt3D,phi3D,theta3D,Phi3DAcc(kFr,k)] = GetSCATR_3D(filename,isarprf,NPulses,dtmoc,toffsetI3);
      end
      
      % Override the average speed with speed of the average components
      SpeedBar          = sqrt(mean(VgpsR)^2+mean(VgpsC)^2);
      
      fprintf(1,'\n Average speed, speed of average components:            %9.4f %9.4f\n',SpeedBar0,SpeedBar);
      
      %---------------------------------------------------------------------------
      
      [~,~,AccG(kFr,k)] = LinReg(VgpsR,tgpsV);
      
      rvGPS(kFr,k)      = mean(VgpsR);
      xrvGPS(kFr,k)     = mean(VgpsC);
      VoverR(kFr,k)     = mean(VgpsC)/Rbar(kFr,k);

      AspGPS(kFr,k)     = atan2d(xrvGPS(kFr,k),rvGPS(kFr,k));
   
      %---------------------------------------------------------------------------
      
      % Map-Drift Equation (RHS and LHS) for a stationary radar at the correct
      % range - using only mean values of the GPS data
      
      fprintf(1,'\n                 Results for GPS Data Only:')
   
      % Replicate Figure 26 in the current draft (Mod 23) of the IAT paper
      
      [AccM(kFr,k),dVM(kFr,k),LHSM(kFr,k),RHSM(kFr,k),xrvM(kFr,k)] = IAT_Test(SpeedBar,AspBar(kFr,k),0,Rbar(kFr,k),tmoc,tgpsV,VgpsR,VgpsC,1,0,[],0);
   
      AspM(k)           = atan2d(xrvM(kFr,k),rvGPS(kFr,k));
   
      %---------------------------------------------------------------------------
   
      fprintf(1,'\n                 Results needed for ISAR Data Analysis Only:')
   
      % Prepare the ISAR motion compensation data for processing
      
      % Compute noise estimate for the motion compensation velocity
      [Lin,m1,a1]       = LinReg(moc32,tmoc);
      mocNoise          = moc32-Lin;
      MocoAcc0(kFr,k)   = a1;
   
      fprintf(1,'\n Mean Rdot, Acceleration, Acc*T, std(noise):            %9.4f %9.4f %9.4f %9.4f\n', ...
                    m1,a1,a1*tmoc(end),std(mocNoise));
   
      % Implement the ISAR AutoTrack data preparation algorithm
      [LinRegC,MocoMean,MocoAcc(kFr,k),Mocomp] = IAT_DataPrepSCATR(tmoc,dtmoc,moc32,mocfiles(k).name(1:end-4),plotMoco);
   
      AccV(kFr,k)       = mean(VgpsC)^2/Rbar(kFr,k);
      
      AccErr            = MocoAcc(kFr,k)-AccV(kFr,k);
      T                 = tmoc(end)-tmoc(1);
      dRmaxA(kFr,k)     = AccErr*T^2/8;
      
      fprintf(1,'\n Acceleration Error, dRmaxA:                            %9.4f %9.4f\n',AccErr,dRmaxA(kFr,k));
      
      % Estimate the mocomp noise by removing the linear model; but then add back
      % in a fraction of the acceleration error for testing.
   
      MocoNoise         = Mocomp-LinRegC+Amult*AccErr*(tmoc-mean(tmoc));
   
      %---------------------------------------------------------------------------
      %---------------------------------------------------------------------------
      
      % Replicate Figure 27 in the current draft (Mod 23) of the IAT paper
      
      LHSQ(kFr,k)       = MocoAcc(kFr,k)*Rbar(kFr,k);                % Left Hand Side: Acc*R
   
      RHSQ(kFr,k)       = SpeedBar^2-RdotBar(kFr,k)^2;               % Right Hand Side: Square of cross-range velocity
      RHSQ(kFr,k)       = max(RHSQ(kFr,k),0);
   
      dVQ(kFr,k)        = MocoAcc(kFr,k)*(tmoc(end)-tmoc(1));
   
      %---------------------------------------------------------------------------
      %---------------------------------------------------------------------------
   
      % Map-Drift Equation (RHS and LHS) for a moving radar at a different range - using only
      % selected mean values of the GPS data - speed, aspect angle and range - and using the
      % radar data. Use a bearing of zero so that the heading is equal to the aspect angle.
      % The reference value for the cross-range velocity is SpeedBar*sind(AspBar(kFr,k)). Ignore
      % the mocomp noise.
   
      fprintf(1,'\n                 Results for a moving radar at a different range (DEMO):')
   
      rvIAT(kFr,k)      = mean(moc32);
      
      % Calculate the map-drift model, ignoring the noise term
      [AccD(k),dVD(k),LHSD(k),RHSD(k),xrvIATD(kFr,k)] = IAT_Test(SpeedBar,AspBar(kFr,k),0,Rbar(kFr,k),tmoc,   ...
                                                                 tgpsV,VgpsR,VgpsC,Rmult,xrpv,MocoNoise,0);
   
      fprintf(1,'\n AccG(kFr,k), AccM(kFr,k), AccD(k):                     %9.4f %9.4f %9.4f\n',AccG(kFr,k),AccM(kFr,k),AccD(k));
      
      AspIATD(kFr,k)    = atan2d(xrvIATD(kFr,k),rvIAT(kFr,k));
      
      AspDifD(kFr,k)    = AspIATD(kFr,k)-AspGPS(kFr,k);
      if AspDifD(kFr,k) >  180; AspIATD(kFr,k) = AspIATD(kFr,k)-360; AspDifD(kFr,k) = AspDifD(kFr,k)-360; end
      if AspDifD(kFr,k) < -180; AspIATD(kFr,k) = AspIATD(kFr,k)+360; AspDifD(kFr,k) = AspDifD(kFr,k)+360; end
   
      fprintf(1,'\n AspGPS(kFr,k), AspIATD(kFr,k), AspDifD(kFr,k):         %9.1f %9.1f %9.1f\n',AspGPS(kFr,k),AspIATD(kFr,k),AspDifD(kFr,k));
      
      %---------------------------------------------------------------------------
      
      fprintf(1,'\n                 Results for a moving radar at a different range (Data):')
   
      % Map-Drift Equation (RHS and LHS) for a moving radar at a different range - using only
      % selected mean values of the GPS data - speed, aspect angle and range - and using the
      % radar data. Use a bearing of zero so that the heading is equal to the aspect angle.
      % The reference value for the cross-range velocity is SpeedBar*sind(AspBar(kFr,k)). Include
      % the mocomp noise.
   
      rvGPS(kFr,k)      = mean(VgpsR);
      xrvGPS(kFr,k)     = mean(VgpsC);
      
      % Calculate the map-drift model, using the noise term
      [AccE(kFr,k),dVE(kFr,k),LHSE(kFr,k),RHSE(kFr,k),xrvIAT(kFr,k)] = IAT_Test(SpeedBar,AspBar(kFr,k),0,Rbar(kFr,k),tmoc, ...
                                                                                tgpsV,VgpsR,VgpsC,Rmult,xrpv,MocoNoise,Nmult);
   
      fprintf(1,'\n AccG(kFr,k), AccM(kFr,k), AccE(kFr,k), MocoAcc(kFr,k): %9.4f %9.4f %9.4f %9.4f\n',AccG(kFr,k),AccM(kFr,k),AccE(kFr,k),MocoAcc(kFr,k));
      
      dRmaxA_E(kFr,k)   = AccE(k)*T^2/8;
      AspGPS(kFr,k)     = atan2d(xrvGPS(kFr,k),rvGPS(kFr,k));
      AspIAT(kFr,k)     = atan2d(xrvIAT(kFr,k),rvIAT(kFr,k));
      
      AspDif(kFr,k)     = AspIAT(kFr,k)-AspGPS(kFr,k);
      if AspDif(kFr,k) >  180; AspIAT(kFr,k) = AspIAT(kFr,k)-360; AspDif(kFr,k) = AspDif(kFr,k)-360; end
      if AspDif(kFr,k) < -180; AspIAT(kFr,k) = AspIAT(kFr,k)+360; AspDif(kFr,k) = AspDif(kFr,k)+360; end
   
      fprintf(1,'\n AspGPS(kFr,k), AspIAT(kFr,k), AspDif(kFr,k):           %9.1f %9.1f %9.1f\n',AspGPS(kFr,k),AspIAT(kFr,k),AspDif(kFr,k));
      
      %---------------------------------------------------------------------------
      
   end  % for kFr = 1:frperleg

   if plotMoco && (icase == 5) % Only done if the plotMoco flag set and using the whole 4 minutes data
      
      % Plot the various mocomp velocities for each file
      figure(300+k);clf
      plot(tmoc,moc32,'b')
      hold on
      plot(tmoc,Mocomp,'g')
      plot(tmoc,LinRegC,'r')
      plot(tgpsV,VgpsR,'m')
      plot(tmoc,MocoNoise,'k')
      grid on
      title(mocfiles(k).name(1:end-4));
      xlabel('Time(s)')
      ylabel('Mocomp Velocity (m/s)')
      legend('ISAR Original','Data Prep','Linear','GPS','Noise')
      axis tight; hold off
      saveas(gcf,['RadialVelocity_' mocfiles(k).name(1:end-4)],'fig')
   
      SpdGPS = sqrt(VgpsR.^2+VgpsC.^2);
      
      % Plot just the GPS versus the ISAR Mocomp
      figure(308+k);clf
      plot(tmoc,moc32-mean(moc32),'b')
      hold on
      plot(tgpsV,VgpsR-mean(VgpsR),'r')
      plot(tgpsV,SpdGPS-mean(SpdGPS),'k')
      grid on
      title(['In-range Velocity, ISAR vs GPS: ' mocfiles(k).name(1:end-4)]);
      xlabel('Time(s)')
      ylabel('In-range Velocity (m/s)')
      legend('ISAR','GPS','SpeedGPS')
      axis tight; hold off
      saveas(gcf,[ofile 'ISARvsGPS_' mocfiles(k).name(1:end-4)],'fig')

      % Plot just the GPS versus the ISAR Mocomp
      figure(316+k);clf
      plot(tmoc,moc32,'b')
      hold on
      plot(tgpsV,VgpsR,'r')
      grid on
      title(['In-range Velocity, ISAR vs GPS: ' mocfiles(k).name(1:end-4)]);
      xlabel('Time(s)')
      ylabel('In-range Velocity (m/s)')
      legend('ISAR','GPS')
      axis tight; hold off
      saveas(gcf,[ofile 'ISARvsGPS_InRangeVel_' mocfiles(k).name(1:end-4)],'fig')

      % Plot just the GPS versus the ISAR Mocomp
      figure(324+k);clf
      plot(tgpsV,VgpsR-mean(VgpsR),'b')
      hold on
      plot(tgpsV,-(SpdGPS-mean(SpdGPS)),'r')
      grid on
      title(['GPS In-range Velocity vs GPS Speed: ' mocfiles(k).name(1:end-4)]);
      xlabel('Time(s)')
      ylabel('Velocity (m/s)')
      legend('GPS Radial Velocity','GPS Speed')
      axis tight; hold off
      saveas(gcf,[ofile 'GPSVrvsGPSSpeed_' mocfiles(k).name(1:end-4)],'fig')

   end

end     % for k = 1:numfiles ... each file is a leg of the octagonal path

fprintf(1,'\n RMS AspDifD, AspDif:  %9.4f %9.4f\n',sqrt(mean(AspDifD(:))^2+var(AspDifD(:))),sqrt(mean(AspDif(:))^2+var(AspDif(:))));

for k = 1:numfiles
   fprintf(1,'\n RMS AspDifD, AspDif:  %9.4f %9.4f\n',sqrt(mean(AspDifD(:,k))^2+var(AspDifD(:,k))),sqrt(mean(AspDif(:,k))^2+var(AspDif(:,k))));
end

for kFr = 1:Frperleg
   fprintf(1,'\n RMS AspDifD, AspDif:  %9.4f %9.4f\n',sqrt(mean(AspDifD(kFr,:))^2+var(AspDifD(kFr,:))),sqrt(mean(AspDif(kFr,:))^2+var(AspDif(kFr,:))));
end

fprintf(1,'\n RMS Error for AccG, MocoAcc0:                          %9.5f %9.5f\n',rms(AccG(:)-AccV(:)),rms(MocoAcc0(:)-AccV(:)));
   
fprintf(1,'\n Correlation G/I, std(G), std(I), std(V):               %9.5f %9.5f %9.5f %9.5f\n',correlation(AccG(:),MocoAcc0(:)),std(AccG(:)),std(MocoAcc0(:)),std(AccV(:)));
   
figure(341);clf
scatter(LHSM(:),RHSM(:),'filled')
hold on
grid on
title ('Map-Drift Equation (Model) for 8 Tybee Cases')
xlabel('Acc*R')
ylabel('u^2')
xequal = min(RHSM(:)):0.1:max(RHSM(:));
plot(xequal,xequal,'r')
saveas(gcf,[ofile 'Figure1_IAT_Equation_Model'],'fig')

figure(342);clf
scatter(LHSD(:),RHSD(:),'filled')
hold on
grid on
title (['Map-Drift Equation (Data) for 8 Tybee Cases, Vp =' sprintf(' %3.0f ',xrpv) ' m/s'])
xlabel('Acc*R')
ylabel('V^2-2Vu+u^2')
xequal = min(RHSD(:)):0.1:max(RHSD(:));
plot(xequal,xequal,'r')
saveas(gcf,[ofile 'Figure2_IAT_Equation_Demo'],'fig')

figure(343);clf
scatter(LHSE(:),RHSE(:),'filled')
hold on
grid on
title (['Map-Drift Equation (Data) for 8 Tybee Cases, Vp =' sprintf(' %3.0f ',xrpv) ' m/s'])
xlabel('Acc*R')
ylabel('V^2-2Vu+u^2')
xequal = min(RHSE(:)):0.1:max(RHSE(:));
plot(xequal,xequal,'r')
saveas(gcf,[ofile 'Figure3_IAT_Equation_Extrapolated'],'fig')

figure(344);clf
scatter(LHSQ(:),RHSQ(:),'filled')
hold on
grid on
title ('Map-Drift Equation (IAT Paper) for 8 Tybee Cases')
xlabel('Acc*R')
ylabel('u^2')
xequal = min(RHSQ(:)):0.1:max(RHSQ(:));
plot(xequal,xequal,'r')
saveas(gcf,[ofile 'Figure4_IAT_Equation_Paper'],'fig')

figure(345);clf
scatter(Rbar(:)/1000,dVE(:),'b','filled')
hold on
scatter(Rbar(:)/1000,dVM(:),'r','filled')
grid on
title ('Integrated Acceleration for 8 USCG Tybee Cases')
legend('Extrapolated','Model')
xlabel('Range (km)')
ylabel('dV (m/s)')
saveas(gcf,[ofile 'Figure5_Integrated_Acceleration_for_Data and Model'],'fig')

figure(346);clf
scatter(AspBar(:),rvGPS(:),'b','filled')
hold on
grid on
scatter(AspBar(:),rvIAT(:),'r','filled')
title ('Rdot versus Aspect angle')
xlabel('GPS Aspect Angle (deg)')
ylabel('Mean Rdot (m/s)')
legend('GPS','ISAR')
xlim([-180 180])

figure(347);clf
scatter(MocoAcc0(:),MocoAcc(:),'b','filled')
hold on
grid on
title ('Map-Drift Equation (Data) for 8 Tybee Cases')
xlabel('Acc Before DataPrep')
ylabel('Acc After DataPrep')
xequal = min(MocoAcc(:)):0.00001:max(MocoAcc(:));
plot(xequal,xequal,'r')

figure(348);clf
scatter(Rbar(:),MocoAcc(:),'b','filled')
grid on
title ('Map-Drift Equation (Data) for 8 Tybee Cases')
xlabel('Range (m)')
ylabel('Acceleration After DataPrep')

figure(349);clf
scatter(xrvGPS(:),xrvM(:),'b','filled')
hold on
grid on
xequal = min(xrvM(:)):0.00001:max(xrvM(:));
plot(xequal,xequal,'r')
title ('Map-Drift Equation (Model) for 8 Tybee Cases')
xlabel('GPS Crossrange Velocity (m/s)')
ylabel('Model Crossrange Velocity (m/s)')

figure(350);clf
scatter(xrvGPS(:),xrvIATD(:),'b','filled')
hold on
grid on
xequal = min(xrvGPS(:)):0.00001:max(xrvGPS(:));
plot(xequal,xequal,'r')
title ('Map-Drift Equation (Demo) for 8 Tybee Cases')
xlabel('GPS Crossrange Velocity (m/s)')
ylabel('IAT Crossrange Velocity (m/s)')

figure(351);clf
scatter(AspGPS(:),AspIATD(:),'b','filled')
hold on
grid on
scatter(AspGPS(:),AspDifD(:),'g','filled')
xequal = -180:180;
plot(xequal,xequal,'r')
xlim([-180 180])
ylim([-180 180])
title (['Map-Drift Equation (Demo) for 8 Tybee Cases, Vp =' sprintf(' %3.0f ',xrpv) ' m/s'])
xlabel('GPS Aspect Angle (deg)')
ylabel('IAT Aspect Angle (deg)')

figure(352);clf
scatter(xrvGPS(:),xrvIAT(:),'b','filled')
hold on
grid on
xequal = min(xrvGPS(:)):0.00001:max(xrvGPS(:));
plot(xequal,xequal,'r')
title (['Map-Drift Equation (Data) for 8 Tybee Cases, Vp =' sprintf(' %3.0f ',xrpv) ' m/s'])
xlabel('GPS Crossrange Velocity (m/s)')
ylabel('IAT Crossrange Velocity (m/s)')

figure(353);clf
scatter(AspGPS(:),AspIAT(:),'b','filled')
hold on
grid on
scatter(AspGPS(:),AspDif(:),'g','filled')
xequal = -180:180;
plot(xequal,xequal,'r')
legend('IAT Aspect','IAT Error','Equal')
xlim([-180 180])
ylim([-180 180])
title (['Map-Drift Equation (Data) for 8 Tybee Cases, Vp =' sprintf(' %3.0f ',xrpv) ' m/s'])
legend('location','northwest')
xlabel('GPS Aspect Angle (deg)')
ylabel('IAT Aspect Angle (deg)')
saveas(gcf,[ofile 'Figure6_Aspect Error'],'fig')

figure(354);clf
scatter(Rbar(:)*Rmult/1000,dRmaxA(:),'b','filled')
hold on
grid on
scatter(Rbar(:)*Rmult/1000,dRmaxA_E(:),'r','filled')
title ('Integrated Acceleration Error (Data) for 8 Tybee Cases')
legend('ISAR after Data Prep','Extrapolated')
xlabel('Range (km)')
ylabel('dRmaxA (m)')
saveas(gcf,[ofile 'Figure7_Integrated Acceleration Error'],'fig')

figure(355);clf
scatter(AspGPS(:),dRmaxA(:),'b','filled')
hold on
grid on
scatter(AspGPS(:),dRmaxA_E(:),'r','filled')
title ('Integrated Acceleration Error (Data) for 8 Tybee Cases')
legend('Point Loma Fixed','Extrapolated')
xlabel('GPS Aspect Angle (deg)')
ylabel('dRmaxA (m)')
saveas(gcf,[ofile 'Figure8_Integrated Acceleration Error'],'fig')

figure(356);clf
scatter(AccG(:),MocoAcc0(:),'b','filled')
hold on
xequal = min(AccG(:)):0.0001:max(AccG(:));
plot(xequal,xequal,'r')
title ('Raw Data Acceleration: GPS vs ISAR')
grid on
xlabel('GPS Acceleration (m/s^2)')
ylabel('ISAR Acceleration (m/s^2)')
saveas(gcf,[ofile 'Figure9_Raw Data Acceleration_GPS vs ISAR'],'fig')

figure(357);clf
scatter(AccG(:),AccM(:),'b','filled')
hold on
xequal = min(AccG(:)):0.0001:max(AccG(:));
plot(xequal,xequal,'r')
title ('Acceleration: Model vs GPS')
grid on
xlabel('GPS Acceleration (m/s^2)')
ylabel('Model Acceleration (m/s^2)')

figure(358);clf
scatter(MocoAcc0(:),AccM(:),'b','filled')
hold on
xequal = min(AccG(:)):0.0001:max(AccG(:));
plot(xequal,xequal,'r')
title ('Acceleration: Model vs ISAR')
grid on
xlabel('ISAR Acceleration (m/s^2)')
ylabel('Model Acceleration (m/s^2)')

figure(359);clf
scatter(Rbar(:),MocoAcc0(:)-AccM(:),'b','filled')
hold on
scatter(Rbar(:),AccG(:)-AccM(:),'r','filled')
legend('ISAR','GPS')
title ('Acceleration Error: ISAR vs GPS')
xlabel('Range (m)')
ylabel('Acceleration (m/s^2)')
grid on; hold off

figure(360);clf
scatter(AspGPS(:),MocoAcc0(:)-AccM(:),'b','filled')
hold on
scatter(AspGPS(:),AccG(:)-AccM(:),'r','filled')
legend('ISAR','GPS')
title ('Acceleration Error: ISAR vs GPS')
xlabel('GPS Aspect (deg)')
ylabel('Acceleration (m/s^2)')
grid on; hold off

figure(361);clf
scatter(AccM(:),AccV(:),'b','filled')
hold on
xequal = min(AccM(:)):0.0001:max(AccM(:));
plot(xequal,xequal,'r')
title ('Acceleration Error: Model vs Vc')
xlabel('Model Acceleration (m/s^2)')
ylabel('Vc Acceleration (m/s^2)')
grid on; hold off

figure(362);clf
scatter(Rbar(:)/1000,abs(AspDif(:)),'b','filled')
title ('IAT Aspect Error vs Range')
xlabel('GPS Range (km)')
ylabel('IAT Aspect Error (deg)')
grid on

figure(363);clf
scatter(AspGPS(:),abs(AspDif(:)),'b','filled')
title ('IAT Aspect Error vs Aspect')
xlabel('GPS Aspect Angle (deg)')
ylabel('IAT Aspect Error (deg)')
xlim([-180 180])
grid on

figure(364);clf
scatter(HdgBar(:),abs(AspDif(:)),'b','filled')
title ('IAT Aspect Error vs Heading')
xlabel('GPS Heading (deg)')
ylabel('IAT Aspect Error (deg)')
xlim([-180 180])
grid on

figure(365);clf
scatter(AccV(:),AccG(:),'b','filled')
hold on
scatter(AccV(:),MocoAcc0(:),'g','filled')
xequal = min(AccV(:)):0.0001:max(AccV(:));
plot(xequal,xequal,'r')
title ('GPS and ISAR Accelerations vs True Value')
xlabel('True Acceleration (m/s^2)')
ylabel('Acceleration (m/s^2)')
legend('GPS','ISAR','Equal','location','southeast')
grid on

figure(366);clf
scatter(AccV(:),MocoAcc0(:),'b','filled')
hold on
xequal = min(AccM(:)):0.0001:max(AccM(:));
plot(xequal,xequal,'r')
title ('ISAR Acceleration vs True Value')
xlabel('True Acceleration (m/s^2)')
ylabel('ISAR Acceleration (m/s^2)')
grid on
saveas(gcf,[ofile 'Figure10_ISAR Acceleration vs True'],'fig')

% Plot 3-D Diagnostics - not related to the IAT
if plotMoco > 1 && (icase == 5)

   figure(367);clf
   scatter(-VoverR(:)*180/pi,Phi3DAcc(:)*180/pi,'b','filled')
   hold on
   xequal = min(-VoverR(:)*180/pi):0.0001:max(-VoverR(:)*180/pi);
   plot(xequal,xequal,'r')
   title ('Rotation Rates')
   xlabel('Rotation Rate GPS (deg/s)')
   ylabel('Rotation Rate 3D (deg/s)')
   grid on
   saveas(gcf,'RotationRates_GPS_3D','fig')

end

