function RF_Win = addMFwin(RF,time,Period,W,train1,train2,apply1,apply2,Sync,PSync)

% Get coefficients of the match filter of the original data
[cosft,sinft] = MFwin(RF,train1,train2,time,Period);

if ~isempty(Sync)

   % Get match filter of the sync function
   [cosftS,sinftS] = MFwin(Sync,train1,train2,time,Period);
   
   sinftn          = sinft /sqrt(sinft^2+cosft^2);
   cosftn          = cosft /sqrt(sinft^2+cosft^2);
   
   sinftSn         = sinftS/sqrt(sinftS^2+cosftS^2);
   cosftSn         = cosftS/sqrt(sinftS^2+cosftS^2);
   
   % Compute weight for original and sync data
   delta           = min(1,1-(sinftn*sinftSn+cosftn*cosftSn));
   
%    fprintf(1,' Original and Sync angles, Delta: %10.4f %10.4f %10.4f\n', ...
%              atan2(sinft,cosft)*180/pi,atan2(sinftS,cosftS)*180/pi,delta);
   
   % Rotate sync terms through offset angle PSync
   costmp          = cosftS;
   sintmp          = sinftS;
   cosftS          = cos(PSync)*costmp-sin(PSync)*sintmp;
   sinftS          = sin(PSync)*costmp+cos(PSync)*sintmp;
   
   % Correct for amplitude of sync function
   ampcor          = sqrt((cosft^2+sinft^2)/(cosftS^2+sinftS^2));
   
   % Blend in phase of the sync function and amplitude of the original data
   cosft           = (1-delta)*cosft+delta*cosftS*ampcor;
   sinft           = (1-delta)*sinft+delta*sinftS*ampcor;

end

% fprintf(1,' cosft, sinft: %10.2f %10.2f\n',1e6*cosft,1e6*sinft);

RF_Win        = zeros(1,apply2+1-apply1);
for k = apply1:apply2
   RF_Win(k+1-apply1) = W(k+1-apply1)*(cosft*cos(2*pi*time(k)/Period)+ ...
                                       sinft*sin(2*pi*time(k)/Period));
end
