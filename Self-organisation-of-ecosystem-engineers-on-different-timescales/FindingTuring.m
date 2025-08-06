% Script to find value of Turing and Saddle-Node bifurcation
clc
aSN = 2*b*m + 2*m*sqrt(b^2+1);
aT = fzero(@(a) minqmu(a,m,b,D),guess);


function minq = minqmu(a,m,b,D)
vp = (a + sqrt(a^2-4*m*(a*b+m)))/2/(a*b+m);
d = m*b*vp/(1-b*vp);
mu_min_squared = -1/2/D * (d+(1+vp^2)*D-m);
mu_min = sqrt(mu_min_squared);
minq = qmu(mu_min,a,m,b,D);
end