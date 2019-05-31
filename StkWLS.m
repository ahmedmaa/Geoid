% This function computes the geoid using least-squares modified Stokes kernel 
% 
%
%
% Inputs:
%         fm : The latitudes truncated within the computation window
%         lm : The longitudes truncated within the computation window
% dfi & dlam : The resolution of the grid in degrees
%         c  : The truncated computation window
%         R  : The Earth's mean radius
%         ng : The normal gravity anomaly computed by Somigliana's formula
%         L  : The maximum modification degree
%         Sn : The file of the modification parameters 
%     
% Output: NW: High-frequency geoid height
%                        
%                      
%
%                            Ahmed Abdalla
%                     Louisiana State University
%                              May 2019
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NW] = StkWLS(fm,lm,dfi,dlam,c,R,ng,L,Sn)
fi=c(:,2);
lam=c(:,1);
% spherical distance
s1=acos(sind(fm).*sind(fi)+cosd(fm).*cosd(fi).*cosd(lam-lm));
% computing the block area A
A1=2*dlam*pi/180*sind(dfi/2)*cosd(fi);
% Stokes kernel
Stk=(1./(sin(s1./2)))-6.*sin(s1./2)+1-5.*cos(s1)-3.*cos(s1).*log(sin(s1./2)+(sin(s1./2)).^2);
t=cos(s1)';
Sn=[0;Sn];
Pn1=zeros(L+1,length(c));
Pn1(1,1:length(c))=1;
Pn1(2,1:length(c))=t';
n1=[1:L+1]';
[P,lsf] = lgpoly(t,L);
yy=(lsf(3:length(lsf),1).*Sn(3:length(lsf)).*P(3:length(lsf),1:length(c)));
HG=sum(yy(:,1:length(c)));
STLS=Stk-HG';
DF=(c(:,3).* STLS.*A1);
T1=sum(DF);
NW=T1*(R/(4*pi*ng));
end


