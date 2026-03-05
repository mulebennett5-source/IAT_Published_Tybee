function [LinRegC,MocoMean,MocoAcc,mocompCoSpdC] = IAT_DataPrepSCATR(t,dt,moc32,FileName,plotMoco)

% Data preparation stage - used for both the search and analytical methods

[L,m,a]        = LinReg(moc32,t);
OrigMean       = m;
OrigNoise      = std(moc32-L);
OrigAcc        = a;

% Smooth the adaptive mocomp radial velocity - modified for a PRF of 200 Hz
moc32S1         = medianfilter(moc32,10);       % Median filter
moc32S2         = c1smooth(moc32S1,20);         % 3-point smoother

% Compute the linear regression model of the adaptive motion compensation velocity
[L,m,a]         = LinReg(moc32S2,t);

% Create a measure of the radial velocity noise
VelNoise        = moc32S2-L;
StdNoise        = std(VelNoise);
   
fprintf(1,'\n OrigNoise, SmoothNoise (m/s) in IAT_DataPrepSCATR:     %9.4f %9.4f\n',OrigNoise,StdNoise);

[wavePeriod,A]  = getWavePeriod(VelNoise,t,0.05,0.5,0,'Adaptive Mocomp');

RMSmin          = std(VelNoise);  % Initial value of score

wavePeriodUse   = wavePeriod;     % Initial guess of wave period

mocompCoSpdC    = VelNoise;       % Start value of noise array

N_Waves         = round((length(t))*dt/wavePeriod);

if N_Waves > 3

   for kk = -5:5
   
      wavePeriodU   = wavePeriod*(1+0.02*kk);
   
      VNoiseChapeau = chapeau(VelNoise,wavePeriodU,t,dt,round((length(t))*dt/wavePeriodU),[],0);
   
      mocompCoSpdCT = VelNoise-VNoiseChapeau;
      mocompRMS     = std(mocompCoSpdCT);
      
      if mocompRMS < RMSmin
   
         RMSmin        = mocompRMS;
         wavePeriodUse = wavePeriodU;
         mocompCoSpdC  = mocompCoSpdCT;
         fprintf(1,'\n wavePeriod, wavePeriodUse, std(mocomp):                %9.4f %9.4f %9.4f',wavePeriod,wavePeriodUse,std(mocompCoSpdC));
   
      end
   
   end

else

   mocompCoSpdC  = A-L;

end

% Add back the linear term
mocompCoSpdC    = mocompCoSpdC+L;

% Compute the linear regression model of the adaptive motion compensation velocity
[LinRegC,MocoMean,MocoAcc] = LinReg(mocompCoSpdC,t);

fprintf(1,'\n Original Mean, Output Mean (m/s):                      %9.4f %9.4f\n',OrigMean,MocoMean);

fprintf(1,'\n Original Acc, Output Acc (m/s^2):                      %9.4f %9.4f\n',OrigAcc,MocoAcc);

% Radial velocity noise estimate
rVelNoiseC      = mocompCoSpdC;

StdNoiseC       = std(rVelNoiseC);

fprintf(1,'\n Chapeau estimate of wave period, std(V,C), Chapeau:    %9.4f %9.4f %9.4f %9.4f %9.4f\n', ...
            wavePeriodUse,StdNoise,StdNoiseC,std(VelNoise-rVelNoiseC),correlation(VelNoise,rVelNoiseC));

if plotMoco > 1
   
   figure(501);clf
   plot(t,moc32,'b')
   hold on
   plot(t,mocompCoSpdC,'r')
   legend('Original','Chapeau')
   xlabel('Time (s)')
   ylabel('Mocomp Velocity (m/s)')
   grid on; hold off
   saveas(gcf,['Chapeau_Model_' FileName],'fig')
   
   figure(502);clf
   n               = numel(moc32);             % Number of time steps
   df              = 1/(n*dt);                    % Frequency step for normal FFT
   f               = (0:n-1)*df;
   spect = 20*log10(abs(fft(moc32)));
   plot(f(2:end),spect(2:end));
   grid on
   xlim([df 1])
   saveas(gcf,['Chapeau_Spectrum_' FileName],'fig')

end
