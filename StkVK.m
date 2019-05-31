% This function computes the geoid using vanicek-kleusberg modified kernel 
% 
%
%Input:
%         fm : The latitudes truncated within the computation window
%         lm : The longitudes truncated within the computation window
% dfi & dlam : The resolution of the grid in degrees
%         c  : The truncated computation window
%         R  : The Earth's mean radius
%         ng : The normal gravity anomaly computed by Somigliana's formula
%         L  : The maximum modification degree
%         tk : The file of the modification parameters 
%     
% Output: NVK: High-frequency geoid height
%                        
%                      
%
%                            Ahmed Abdalla
%                     Louisiana State University
%                              May 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NVK] = StkVK(fm,lm,dfi,dlam,c,R,ng,L,tk)
fi=c(:,2);
lam=c(:,1);
% spherical distance
s1=acos(sind(fm).*sind(fi)+cosd(fm).*cosd(fi).*cosd(lam-lm));
% computing the block area A
A1=2*dlam*pi/180*sind(dfi/2)*cosd(fi);
% Stokes kernel
Stk=(1./(sin(s1./2)))-6.*sin(s1./2)+1-5.*cos(s1)-3.*cos(s1).*log(sin(s1./2)+(sin(s1./2)).^2);
t=cos(s1)';
tk=[0;tk];
n1=[1:L+1]';
[P,wgf, lsf] = lgpoly(t,L);
yy1=(wgf(3:length(wgf),1).*Pnn(3:length(wgf),1:length(c)));
yy2=(lsf(3:length(lsf),1).*tk(3:length(lsf)).*P(3:length(lsf),1:length(c)));
SVK1=sum(yy1(:,1:length(c)));
SVK2=sum(yy2(:,1:length(c)));
SVK=Stk-SVK1';
SVKM=SVK-SVK2';
DF=(c(:,3).* SVKM.*A1);
T1=sum(DF);
NVK=T1*(R/(4*pi*ng));
end


% end

