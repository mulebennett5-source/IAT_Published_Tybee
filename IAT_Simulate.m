
function IAT_Simulate(SimOption,IAT_Method,AIS_LatLon,AIS_CourseSpd)

outdir        = ['Si0' filesep];      
fidIAT        = fopen([outdir 'IAT_Si0.txt'],'wt');
S             = fopen([outdir 'IAT_Si0_S.txt'],'wt');

set(0,'DefaultAxesTitleFontWeight','bold')
set(0,'DefaultAxesFontWeight','bold')
set(0,'defaultLineLineWidth',2)

% 22 NotAISAided cases
CasesNAA = char(            ...
   'elandraspruce1.'      , ...
   'transsibbridge2.'     , ...
   'atlanticsail1.'       , ...
   'maerskdenver1.'       , ...
   'maerskshams2.'        , ...
   'colomboexpress2.'     , ...
   'stoltcreativity2.'    , ...
   'msckrystal4.'         , ...
   'shouchenshan1.'       , ...
   'independenthorizon2.' , ...
   'overseaschinook3.'    , ...
   'seamelody4.'          , ...
   'seawaysredwood1.'     , ...
   'seawaysredwood2.'     , ...
   'coscojasmine1.'       , ...
   'mscheidi3.'           , ...
   'sagaadventure3.'      , ...
   'mediatlantico6.'      , ...
   'angeles4.'            , ...
   'sagaadventure4.'      , ...
   'sagaadventure5.'      , ...
   'georgmaersk1.'        );

% For all cases assume that the track estimate of the ship is at range 100 km,
% azimuth 45 degrees, speed 6 m/s
tgtCueRngKm   = 100;
tgtCueAzDeg   = 45;

% True speed
SpeedAIS      = 6;

prf           = 512;
dt            = 1/prf;

Npulse        = 32*prf+101-52;

pulseTov      = (0:Npulse)*dt-Npulse*dt/2;
pulseTov1     = (0:Npulse+1)*dt-Npulse*dt/2;

% Target position relative to platform at Cue Time
tgtE0Cue      = tgtCueRngKm*cosd(tgtCueAzDeg);        % Target East  (km)
tgtN0Cue      = tgtCueRngKm*sind(tgtCueAzDeg);        % Target North (km)

% Platform track - heading East at 100 m/s and assumed zero offset at first pulse
platHdg       = 90;
platSpd       = 100;
platEvel      = platSpd*sind(platHdg);
platNvel      = platSpd*cosd(platHdg);

platEtrack    = platEvel*pulseTov1;
platNtrack    = platNvel*pulseTov1;
   
% Arrays for output plot
Asp           = zeros(3,17);
Asp13         = zeros(3,17);

HdgErr        = 0;
SpdErr        = 0;

NDifKm        = 0;
EDifKm        = 0;

for kErr = 1:3
   
   if     SimOption == 1
      
      HdgErr      = 25*(kErr-2);
   
   elseif SimOption == 2
      
      CDifKm      = 1*(kErr-2);
      EDifKm      = CDifKm*abs(cosd(tgtCueAzDeg));
      NDifKm      = CDifKm*abs(sind(tgtCueAzDeg));
   
   elseif SimOption == 3

      SpdErr      = 2*(kErr-2);
   
   elseif SimOption == 4
      
      HdgErr      = 25;
      CDifKm      = 1;
      
      tgtCueAzDeg = 85-40*(kErr-1);  % Azimuth = 85,45,5 Degrees
      
      EDifKm      = CDifKm*abs(sind(tgtCueAzDeg));
      NDifKm      = CDifKm*abs(cosd(tgtCueAzDeg));
      
      SpdErr      = 2;
   
   end
   
   % Loop over true ship heading, 5 to 85 degrees
   for ihdg = 1:17
   
      SpeedAux       = SpeedAIS+SpdErr;
      
      % True heading and aspect
      tgtHdgAISDeg   = 5*ihdg+tgtCueAzDeg;
      Asp(kErr,ihdg) = tgtHdgAISDeg-tgtCueAzDeg;
      
      % East and North target velocities for AIS
      tgtEvelAIS     = SpeedAIS*sind(tgtHdgAISDeg);
      tgtNvelAIS     = SpeedAIS*cosd(tgtHdgAISDeg);
      
      % Tracker heading
      tgtHdgAuxDeg   = tgtHdgAISDeg+HdgErr;
   
      % East and North target velocities for Aux
      tgtEvelAux     = SpeedAux*sind(tgtHdgAuxDeg);
      tgtNvelAux     = SpeedAux*cosd(tgtHdgAuxDeg);
   
      % Loop over noise file
      for icase = 13:13
      
         if icase == 0
            VelNoise        = zeros(1,Npulse+1);
            wavePeriod      = Npulse*dt/4.25;
            for k = 1:Npulse+1
               VelNoise(k) = 1.0*sin(2*pi*k*dt/wavePeriod);
            end
         else
            NoiseFile       = fopen(['VelNoise\' dotlim(CasesNAA(icase,:)) '_VelNoise.dat'],'rb');
            VelNoise        = fread(NoiseFile,16434,'float32')';
            fclose(NoiseFile);
         end
         
         % Target East and North coordinates at the initial pulse time
         tgtE0           = tgtE0Cue;
         tgtN0           = tgtN0Cue;
         
         % Set the AIS initial coordinates to the cue values
         tgtE0AIS        = tgtE0Cue+EDifKm;
         tgtN0AIS        = tgtN0Cue+NDifKm;
         
         % Base file name for output
         outfile         = [outdir sprintf('%2.0f',10*tgtHdgAISDeg)];
      
         % True range to target versus time
         mocompPointAIS  = sqrt((tgtE0AIS*1000+tgtEvelAIS*pulseTov1-platEtrack).^2 ...
                               +(tgtN0AIS*1000+tgtNvelAIS*pulseTov1-platNtrack).^2);
      
         mocompSpdAIS    = diff(mocompPointAIS)/dt;
         
         % Tracker range to target versus time
         mocompPointAux  = sqrt((tgtE0*1000+tgtEvelAux*pulseTov1-platEtrack).^2 ...
                               +(tgtN0*1000+tgtNvelAux*pulseTov1-platNtrack).^2);
      
         mocompDetSpdAux = diff(mocompPointAux)/dt;
      
         % Adaptive mocomp value versus time for Aux
         mocompCoSpdAux  = mocompSpdAIS-mocompDetSpdAux+VelNoise;
   
         % Adaptive mocomp value versus time for AIS
         mocompCoSpdAIS  = mocompSpdAIS-mocompDetSpdAux+VelNoise;
   
         mocompPointG    = dt*cumsum(mocompCoSpdAux+mocompDetSpdAux);
   
         [~,m1,a1]       = LinReg(mocompCoSpdAux,pulseTov);
         fprintf(1,'\n    Mean and Linear Terms for Aux, Acc*Ti: %8.4f %8.4f %8.4f\n',m1,a1,a1*pulseTov(end));
         
         [~,m2,a2]       = LinReg(mocompCoSpdAIS,pulseTov);
         fprintf(1,  '    Mean and Linear Terms for AIS, Acc*Ti: %8.4f %8.4f %8.4f\n',m2,a2,a2*pulseTov(end));

         % Call the ISAR AutoTrack routine to process the data
         
         IATHdg          = IAT_Est(tgtE0,tgtN0,tgtHdgAuxDeg,SpeedAux,tgtHdgAISDeg,SpeedAIS, ...
                                   pulseTov,dt,platEtrack(1:end-1),platNtrack(1:end-1),     ...
                                   mocompPointG,mocompDetSpdAux,mocompCoSpdAux,fidIAT,      ...
                                   outfile,tgtE0AIS,tgtN0AIS,0,S,IAT_Method,AIS_LatLon,     ...
                                   AIS_CourseSpd);
   
         if icase == 13 || icase == 0; Asp13(kErr,ihdg) = IATHdg-tgtCueAzDeg; end  
   
      end  % for icase = 0:0 (simulated radial velocity noise)

   end     % for ihdg  = 1:17

end        % for kErr  = 1:3

figure(2001);clf
plot(Asp(2,:),Asp13(2,:),'g')
hold on
plot(Asp(2,:),Asp13(1,:),'k')
plot(Asp(2,:),Asp13(3,:),'m')
title ('Aspect Accuracy for the ISAR AutoTrack Algorithm')
xlabel('True Aspect (Deg)')
ylabel('IAT Estimate (Deg)')
if     SimOption == 1
   legend('Correct Aspect','Aspect -25 Deg','Aspect +25 Deg')
   legend('Location','southeast')
elseif SimOption == 2
   legend('Correct CrossRange','CrossRange -1 km','CrossRange +1 km')
   legend('Location','southeast')
elseif SimOption == 3
   legend('Correct Speed','Speed -2 m/s','Speed +2 m/s')
   legend('Location','southeast')
elseif SimOption == 4
   legend('Azimuth 45','Azimuth 85','Azimuth 5')
   legend('Location','southeast')
end
axis tight; grid on

if     SimOption == 1
   saveas(gcf,'IAT_AspectError1','png')
elseif SimOption == 2
   saveas(gcf,'IAT_CrossRangeError1','png')
elseif SimOption == 3
   saveas(gcf,'IAT_SpeedError1','png')
elseif SimOption == 4
   saveas(gcf,'IAT_Azimuth1','png')
end

figure(2002);clf
hold on
plot(Asp(2,:),Asp13(2,:)-Asp(2,:),'g')
plot(Asp(2,:),Asp13(1,:)-Asp(2,:),'k')
plot(Asp(2,:),Asp13(3,:)-Asp(2,:),'m')
title ('Aspect Accuracy for the ISAR AutoTrack Algorithm')
xlabel('True Aspect (Deg)')
ylabel('IAT Aspect Error (Deg)')
if     SimOption == 1
   legend('Correct Aspect','Aspect -25 Deg','Aspect +25 Deg')
   legend('Location','southwest')
elseif SimOption == 2
   legend('Correct CrossRange','CrossRange -1 km','CrossRange +1 km')
   legend('Location','northeast')
elseif SimOption == 3
   legend('Correct Speed','Speed -2 m/s','Speed +2 m/s')
   legend('Location','northeast')
elseif SimOption == 4
   legend('Azimuth 45','Azimuth 85','Azimuth 5')
   legend('Location','southwest')
end
axis tight; grid on

if     SimOption == 1
   saveas(gcf,'IAT_AspectError2','png')
elseif SimOption == 2
   saveas(gcf,'IAT_CrossRangeError2','png')
elseif SimOption == 3
   saveas(gcf,'IAT_SpeedError2','png')
elseif SimOption == 4
   saveas(gcf,'IAT_Azimuth2','png')
end

figure(2003);clf
hold on
plot(Asp(2,:),(Asp13(2,:)-Asp(2,:))*(pi/180)*100.*tand(Asp(2,:)),'g')
plot(Asp(2,:),(Asp13(1,:)-Asp(2,:))*(pi/180)*100.*tand(Asp(2,:)),'k')
plot(Asp(2,:),(Asp13(3,:)-Asp(2,:))*(pi/180)*100.*tand(Asp(2,:)),'m')
title ('LOA Accuracy for the ISAR AutoTrack Algorithm')
xlabel('True Aspect (Deg)')
ylabel('IAT LOA Error (%)')
if     SimOption == 1
   legend('Correct Aspect','Aspect -25 Deg','Aspect +25 Deg')
   legend('Location','southeast')
elseif SimOption == 2
   legend('Correct CrossRange','CrossRange -1 km','CrossRange +1 km')
   legend('Location','northwest')
elseif SimOption == 3
   legend('Correct Speed','Speed -2 m/s','Speed +2 m/s')
   legend('Location','southeast')
elseif SimOption == 4
   legend('Azimuth 45','Azimuth 85','Azimuth 5')
   legend('Location','southwest')
end
axis tight; grid on

if     SimOption == 1
   saveas(gcf,'IAT_AspectError3','png')
elseif SimOption == 2
   saveas(gcf,'IAT_CrossRangeError3','png')
elseif SimOption == 3
   saveas(gcf,'IAT_SpeedError3','png')
elseif SimOption == 4
   saveas(gcf,'IAT_Azimuth3','png')
end

fclose(fidIAT);