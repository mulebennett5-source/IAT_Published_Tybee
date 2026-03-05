function [Aspect,Bearing,Heading,VIR,VCR] = GetAspect(edot,ndot,east,north)

Bearing = atan2d(east,north);
Heading = atan2d(edot,ndot);

Aspect  = Bearing-Heading;
for k = 1:numel(Aspect)
   if Aspect(k) >  180; Aspect(k) = Aspect(k)-360; end
   if Aspect(k) < -180; Aspect(k) = Aspect(k)+360; end
end

%--------------------------------------------------------------------------
% In-range and cross-range velocities
Xb      = sind(Bearing);
Yb      = cosd(Bearing);

VIR     =  Yb.*ndot+Xb.*edot;
VCR     = -Yb.*edot+Xb.*ndot;
%--------------------------------------------------------------------------

% Check on the algorithm for aspect angle
% Aspect2 = atan2d(VCR,VIR);  % The same as Aspect to machine accuracy

%--------------------------------------------------------------------------
% In-track and cross-track velocities - not returned or used
% Xh      = sind(Heading);
% Yh      = cosd(Heading);
% 
% VIT     =  Yh.*ndot+Xh.*edot;
% VCT     = -Yh.*edot+Xh.*ndot;
%--------------------------------------------------------------------------
