function [RF_Out,twidth] = chapeau(RF_In,Period,time,dtbar,N,Sync,PSync)

% This is a bandpass filter that works in the time domain

M             = numel(time);

nmiddle       = N-2;

if isempty(PSync); PSync = 0; end

% Don't do the algorithm if there are not enough blocks of data
if nmiddle < 1
   twindow = ones(M,1);
   RF_Out  = addMFwin(RF_In,time,Period,twindow,1,M,1,M,Sync,PSync);
   twidth  = M*dtbar;
%    RF_Out = RF_In;
%    twidth = 0;
   return
end

nguard        = round(Period/dtbar);
while mod(M-2*nguard,nmiddle) ~= 0 && nguard > 0
   nguard = nguard-1;
end

if nguard == 0
   nguard = round(Period/dtbar);
   while mod(M-2*nguard,nmiddle) ~= 0
      nguard = nguard+1;
   end
end

half          = round((M-2*nguard)/nmiddle);

if half < 4
   RF_Out = RF_In;
   twidth = [];
   return
end

twindow       = ones(2*half,1);
for k = 1:half
   twindow(k)          = (k-0.5)/half;
   twindow(2*half+1-k) = twindow(k);
end

twidth        = 2*half*dtbar;

%twindow       = (1-cos(pi*twindow))/2;  % Cosine on a pedestal weight

lwindow       = ones(1,half+nguard);
rwindow       = ones(1,half+nguard);
for k = 1:half
   lwindow(k+nguard) = twindow(k+half);
   rwindow(k)        = twindow(k);
end

RF_Out        = zeros(1,M);

first         = 1:half+nguard;
last          = M+1-nguard-half:M;

RF_Out(first) = RF_Out(first)     +addMFwin(RF_In,time,Period,lwindow,            ...
                                            nguard+1,half+nguard,                 ...
                                            1,half+nguard,Sync,PSync);

for kmiddle = 1:nmiddle-1
   midoff         = (kmiddle-1)*half;
   middle         = nguard+midoff+1:nguard+midoff+2*half;
   RF_Out(middle) = RF_Out(middle)+addMFwin(RF_In,time,Period,twindow,            ...
                                            nguard+midoff+1,nguard+midoff+2*half, ...
                                            nguard+midoff+1,nguard+midoff+2*half,Sync,PSync);
end

RF_Out(last)  = RF_Out(last)      +addMFwin(RF_In,time,Period,rwindow,            ...
                                            M+1-nguard-half,M-nguard,             ...
                                            M+1-nguard-half,M,Sync,PSync);
