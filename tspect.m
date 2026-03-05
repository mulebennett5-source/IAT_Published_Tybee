function [xspecf,f,xdotuse] = tspect(Nf,xdot,t,dt,nsmooth)

xdotuse              = xdot;
N                    = numel(xdotuse);

if mod(N,2)
   N       = N-1;
   xdotuse = xdotuse(1:end-1);
   t       = t(1:end-1);
end
noff                 = floor((Nf-N)/2);

[xLin,~,~]           = LinReg(xdotuse,t);

% Remove mean and trend
xdotuse              = xdotuse-xLin;

xspec                = zeros(1,Nf);

xspec(1+noff:noff+N) = xdotuse(1:N);

xspecf               = abs(fftshift(fft(fftshift(xspec))))*dt;

fprintf(1,'\n Input and Output Variances: %12.6f %12.6f\n',var(xdotuse),mean(xspecf.^2));

if nsmooth; xspecf = c1smooth(xspecf,nsmooth); end

f                    = (-Nf/2:-1+Nf/2)/(Nf*dt);
