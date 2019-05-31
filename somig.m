function [gamma]=somig(phi)
% To compute the normal gravity (in mGal) using Somigliana's formula.       
k=0.001931851353; %normal gravity constant
e1=0.081819191042811; %First eccentricity
gamma_e=978032.67715;%normal gravity value on the equator (m/sec^2)
a=sind(phi);
e2=e1*e1;
gamma=gamma_e*(1+k*(a*a))./(sqrt(1-e2*(a*a)));
