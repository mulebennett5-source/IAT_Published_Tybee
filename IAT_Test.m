function [Acc,dV,LHS,RHS,xrv] = IAT_Test(SpeedBar,Heading,Bearing,Rbar,tmoc,tgpsV,VgpsR,VgpsC,Rmult,xrpv,MocoNoise,Nmult)

%------------------------------------------------------------------------------
%
% Evaluate the solution of the map-drift equation
%
% SpeedBar  - Mean ship speed (m/s)
%
% Heading   - Ship heading (deg)
%
% Bearing   - Bearing to ship (deg)
%
% Rbar      - Mean range to ship before multiplying by Rmult (m)
%
% tmoc      - Times of pulses (sec)
%
% tgpsV     - Times for the GPS data
%
% VgpsR     - In-range velocity of the ship from the GPS data (m/s)
%
% VgpsC     - Cross-range velocity of the ship from the GPS data (m/s)
%
% Rmult     - Multiple of range for the extrapolated data
%
% xrpv      - Synthetic cross-range platform velocity (m/s)
%
% MocoNoise - Adaptive motion compensation noise from the ISAR processor (m/s)
%
% Nmult     - Multiplier for motion compensation noise
%
%------------------------------------------------------------------------------

% East and North velocities of the target
Edotbar        = SpeedBar*sind(Heading);
Ndotbar        = SpeedBar*cosd(Heading);

% xrpv is the synthetic cross-range platform velocity, assumed to be
% perpendicular to the mean bearing
EdotPlat       = xrpv*sind(Bearing+90);
NdotPlat       = xrpv*cosd(Bearing+90);

% Adjust time to have center zero
tprime         = tmoc-mean(tmoc);
dtmoc          = mean(diff(tmoc));

% East and North displacement from the radar for the ship and platform
E              = Rbar*Rmult*sind(Bearing)+(Edotbar-EdotPlat)*tprime;
N              = Rbar*Rmult*cosd(Bearing)+(Ndotbar-NdotPlat)*tprime;

% Range and range velocity versus time
R              = sqrt(E.^2+N.^2);
Rdot           = diff(R)/dtmoc;

if ~isempty(MocoNoise) && Nmult > 0
   Rdot = Rdot+Nmult*MocoNoise(1:end-1);
end

% Linear model of the range velocity
[~,~,Acc]      = LinReg(Rdot,tprime(1:end-1));

dV             = Acc*(tmoc(end)-tmoc(1));

if dV < 0
   fprintf(1,'\n Negative acceleration: %8.4f',dV);
end

% Left-hand side of the map-drift equation
LHS            = max(0,Rbar*Rmult*Acc);

% Cross-range velocity estimates from the two roots of the quadratic equation
xrvPlus        = xrpv+sqrt(LHS);
xrvMinus       = xrpv-sqrt(LHS);

% Guidance value of the cross-range velocity, equal to half the correct value
% This is only used to choose the correct root of the quadratic equation
VgpsC_Half     = mean(VgpsC)/2;

% Choose the root that is closest to the guidance value
if abs(xrvPlus-VgpsC_Half) < abs(xrvMinus-VgpsC_Half)
   xrv = xrvPlus;
else
   xrv = xrvMinus;
end

% Right-hand side of the map-drift equation
RHS            = xrpv^2-2*xrpv*xrv+xrv^2;

fprintf(1,'\n Crossrange velocity for IAT, GPS, Plus and Minus:      %9.4f %9.4f %9.4f %9.4f',xrv,mean(VgpsC),xrvPlus,xrvMinus);
fprintf(1,'\n LHS, RHS, xrv:                                         %9.4f %9.4f %9.4f\n',LHS,RHS,xrv);
