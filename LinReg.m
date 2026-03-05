function [rdotLin,m,a] = LinReg(rdot,t0)

t       = t0-mean(t0);

m       = mean(rdot);

a       = mean((rdot-m).*t)/var(t);

rdotLin = m+a*t;

end
