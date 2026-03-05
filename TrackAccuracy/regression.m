function regress = regression(u,v)

u       = u-mean(u);
v       = v-mean(v);
regress = sum(u.*v)/sum(u.*u);