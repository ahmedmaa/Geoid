% This function computes legendre polynomial
% 
%
%
% Inputs:
%         t : cos(theta)
%         L : the maximum degree of expansion
% 
%
% Output: 
%         Pn: The Legendre Polynomial
%        wgf: factorial for WG Kernel
%        lsf: factorial for VK and LSM knernels
%                        
%                      
%
%                            Ahmed Abdalla
%                     Louisiana State University
%                              May 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,wgf, lsf] = lgpoly(t,L)
P=zeros(L+1,length(t));
wgf=zeros(L+1,1);
lsf=zeros(L+1,1);
P(1,1:length(t)) = 1;    
P(2,1:length(t)) = t; 
for i = 2:L 
    P(i+1,:) = -(i-1)./i.*P(i-1)+(2.*i-1)./i.*t.*P(i);
    wgf(i+1,1) = (2.*i+1)./(i-1);
    lsf(i+1,1) = (2.*i+1)./2;
end
