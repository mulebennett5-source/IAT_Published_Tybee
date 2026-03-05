function [cosft,sinft] = MFwin(RF,train1,train2,time,Period)

cosft   = 0;
sinft   = 0;
cosftsq = 0;
sinftsq = 0;

for k = train1:train2
   cosft   = cosft  +RF(k)*cos(2*pi*time(k)/Period);
   cosftsq = cosftsq+      cos(2*pi*time(k)/Period)^2;
   sinft   = sinft  +RF(k)*sin(2*pi*time(k)/Period);
   sinftsq = sinftsq+      sin(2*pi*time(k)/Period)^2;
end

cosft   = cosft/cosftsq;
sinft   = sinft/sinftsq;
