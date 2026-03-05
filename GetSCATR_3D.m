function [time3D,dt3D,phi3D,theta3D,Phi3DAcc] = GetSCATR_3D(filename,isarprf,NPulses,dtmoc,toffsetI3)

%-------------------------------------------------------------------------------

% Load the 3-D output from Focus3D - 3-D parameters, Moments. and values for
% LOA estimation

load(['Data3D\' filename(1:end-4) '_3d.3Doutput.mat'], 'time','theta','omega','omegadot','phi','phidot','phidotdot');

load(['Data3D\' filename(1:end-4) '_0_Moments.mat'],   'CovRF_data','CovFF_data','CovRA_data', ...
                                                       'CovFA_data','CovAA_data','D_data');
  
load(['Data3D\' filename(1:end-4) '_EstLOA_Input.mat'],'time','nFrame','NA','rsR','phi0_Sign','theta0_User', ...
                                                       'UsePGAtgts','outfile','filenumber','EstLOA');

%-------------------------------------------------------------------------------

% First and last times for 3D information, allowing buffers at each end for
% smoothing
FirstPulse3    = toffsetI3*isarprf;
tfirst         = (FirstPulse3-1)*dtmoc;
tlast          = tfirst+(NPulses-1)*dtmoc;

time3D         = time    ((time   >= tfirst));
time3D         = time3D  ((time3D <= tlast));
phi3D          = phi     ((time   >= tfirst));
phi3D          = phi3D   ((time3D <= tlast));
theta3D        = theta   ((time   >= tfirst));
theta3D        = theta3D ((time3D <= tlast));

dt3D           = mean(diff(time3D));
time3D         = time3D-time3D(1)+dt3D/2;

[~,~,Phi3DAcc] = LinReg(phi3D,time3D);

fprintf(1,'\n Means CovRF, phi-term, theta-term:  %12.6f %12.6f %12.6f\n',mean(CovRF_data),mean(tan(phi).*phidot),mean(theta.*omega));
fprintf(1,'\n Stds. CovRF, phi-term, theta-term:  %12.6f %12.6f %12.6f\n',std(CovRF_data),std(tan(phi).*phidot),std(theta.*omega));

figure(201);clf
plot(time,CovRF_data,'b')
hold on
plot(time,tan(phi).*phidot,'r')
plot(time,tan(theta).*omega,'g')
xlabel('Time (s)')
ylabel('CovRF or Equivalent')
legend('CovRF','tan(phi).*phidot','tan(theta).*omega')
title ([filename(1:end-4) ': Focus3D RF-Covariance'])
axis tight; grid on; hold off
saveas(gcf,['Focus3D_RF-Covariance_' filename(1:end-4)],'fig')

figure(202);clf
plot(time,CovRF_data+tan(theta).*omega,'b')
hold on
plot(time,-tan(phi).*phidot,'r')
xlabel('Time (s)')
ylabel('CovRF or Equivalent')
legend('CovRF Mod','tan(phi).*phidot')
title ([filename(1:end-4) ': Focus3D RF-Covariance'])
axis tight; grid on; hold off
saveas(gcf,['Focus3D_RF-Covariance_2' filename(1:end-4)],'fig')

figure(203);clf
plot(time,theta,'b')
hold on
plot(time,omega,'r')
plot(time,theta.*omega,'g')

% 
% EstLOA                   = EstimateLOA(fidall,time,nFrame,NA,rsR,phi0_User, ...
%                                        theta0_User,UsePGAtgts,outfile,      ...
%                                        filenumber,figType,[]);
% 
%-------------------------------------------------------------------------------
