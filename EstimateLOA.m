function EstLOA = EstimateLOA(fidall,time,nFrame,N,rsR,phi,theta,UsePGAtgts, ...
                              outfile,filenumber,figType,BadFit)

lxyzUMin = zeros(1,nFrame);
lxyzUMax = zeros(1,nFrame);
LOA      = zeros(1,nFrame);

coscos   = abs(cos(phi)).*abs(cos(theta));

if numel(coscos) == 1
   fprintf(fidall,'\n Length analysis using mean angles');
   for k = 1:nFrame
      lxyzUMin(k) = min(rsR(1:N(k),k))/coscos;
      lxyzUMax(k) = max(rsR(1:N(k),k))/coscos;
      LOA(k)      = lxyzUMax(k)-lxyzUMin(k);
   end
else
   fprintf(fidall,'\n Length analysis using time-dependent angles');
   for k = 1:nFrame
      lxyzUMin(k) = min(rsR(1:N(k),k))/coscos(k);
      lxyzUMax(k) = max(rsR(1:N(k),k))/coscos(k);
      LOA(k)      = lxyzUMax(k)-lxyzUMin(k);
   end
end

if ~isempty(BadFit) 
   LOA_BF = LOA(BadFit<2);
else
   LOA_BF = LOA;
end

EstLOA   = median(LOA_BF(LOA_BF>median(LOA_BF)));

% Rule of Thumb estimate of ship beam versus Length Overall (LOA)

% Some literature on relation of beam to length
% https://en.wikipedia.org/wiki/Beam_(nautical)#Rule_of_thumb_-_formula
% https://en.wikipedia.org/wiki/Panamax
% https://en.wikipedia.org/wiki/Seawaymax

LOAft    = EstLOA/0.3048;               % LOA in feet
Beamft   = 1+LOAft^(2/3);               % Rule of Thumb beam in feet
BeamROT  = Beamft*0.3048;               % Beam in meters
L_Corr   = (BeamROT/2)*tan(abs(mean(phi)));
  
EstLOA   = EstLOA-L_Corr;

figure(167);clf
plot(time,LOA,'b')
hold on
plot(time,lxyzUMax,'r')
plot(time,lxyzUMin,'g')
plot(time,median(LOA)*ones(1,nFrame),'--b')
plot(time,EstLOA*ones(1,nFrame),'y')
plot(time,median(lxyzUMax)*ones(1,nFrame),'--r')
plot(time,median(lxyzUMin)*ones(1,nFrame),'--g')
legend('ThinShip','TS\_Max','TS\_Min','med(LOA)','EstLOA')
xlabel('Time (s)')
ylabel('Length (m)')
axis tight; grid on; hold off
saveas(gcf,[outfile '_rminmax'],figType)
saveas(gcf,['Length\' sprintf('%i',filenumber) '_' sprintf('%i',UsePGAtgts) '_RMinMax'],'fig')
   
fprintf(fidall,'\n Width Correction:                              %8.2f',        L_Corr);
fprintf(fidall,'\n TS Length Double Median:                       %8.2f\n',      EstLOA);
fprintf(fidall,  ' Std. Dev. of Rmin, Rmax:                       %8.2f %8.2f\n',std(lxyzUMin),std(lxyzUMax));
