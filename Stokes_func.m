% This function computes the geoid using original Stokes kernel 
% Inputs:
%         fm : The latitudes truncated within the computation window
%         lm : The longitudes truncated within the computation window
% dfi & dlam : The resolution of the grid in degrees
%         c  : The truncated computation window
%         R  : The Earth's mean radius
%         ng : The normal gravity anomaly computed by Somigliana's formula
%         
%     
% Output: N1: High-frequency geoid height
%                        
%                      
%
%                            Ahmed Abdalla
%                     Louisiana State University
%                              May 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N1] =Stokes_func(fm,lm,dfi,dlam,c,R,ng)
fi=c(:,2);
lam=c(:,1);
s1=acos(sind(fm).*sind(fi)+cosd(fm).*cosd(fi).*cosd(lam-lm));
A1=2*dlam*pi/180*sind(dfi/2)*cosd(fi);
Sepsi1=(1./(sin(s1./2)))-6.*sin(s1./2)+1-5.*cos(s1)-3.*cos(s1).*log(sin(s1./2)+(sin(s1./2)).^2);
DF=(c(:,3).*Sepsi1.*A1);
T1=sum(DF);
N1=T1*(R/(4*pi*ng));
end

