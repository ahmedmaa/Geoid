%**************************************************************************
%****************************** Geoid Software ****************************
%**************************************************************************
%
% This program script uses Stokes's Kernel to evaluate the terrestrial
% gravity measurements and then uses Stokeks kernel to compute the geoid
% model
%
%                               Ahmed Abdalla
%                        Assistant Professor Research
%                         Louisiana State University
%                          Center for Geoinformatics
%                             aabdalla1@lsu.edu
%                                 May 2019
%                          
%                        Last update: April 2022
%
%
%**************************************************************************
%
%    User's guide
%    
%      Data Entry:
% 
%         1- Gravity data should be in the following format [lambda phi dg]
%
%         2- Define the outer offset to determine the size of the computation
%            area(inner grid)in arc-degrees. 
%
%         3- Insert the resolution of the gravity grid in phi & lambda
%
%         4- This software hadnles the original Stokes and other 5 methods
%            of the modiefied Stokes kernels. You will have an option to
%            select the number of the method you prefer in your
%            calculations [1 to 6].
%            *When selecting the modification methods from 3 to 6 you will
%            need to include the modification parameter files which are
%            described as follows:
%            a- For Vanicek-Kleusberg (no 3), the parameter file (tk.prn) is
%            obtained by Featherstone, W: (2003) Software for computing five 
%            existing types of deterministically modified integration kernel 
%            for gravimetric geoid determination
%            b- For Least-squares modification of Stokes, the paramteters
%            files (Sn_bias.prn, Bn_unb.prn and Bn_opt.prn) are obtained
%            from Ellmann, A (2005): Computation of three stochastic 
%            modifications of Stokes's formula for regional geoid 
%            determination
%            
%         5- Enter the name of the result file e.g. Result.TXT
%
%         6- The reuslts file will be in the following [lambda phi dg]
%
%        
%**************************************************************************
clc
clear all
anom=input('Insert the gravity data file [lambda phi dg]: (example: KRTOUT3.dat) \n','s');
fprintf ('\n WAIT! Loading the gravity data grid in [lambda phi dg] \n');
anomaly=load(anom);
R= 6378137.0;            % m
e2=0.00669437999014;
% R=6370000; %m
%Specifying parameters of the gravity file based on their locations
%(colmuns)
dg=anomaly(:,3);
lambda=anomaly(:,1);
fi=anomaly(:,2);
% calculating the number of rows in the gravity file 
n=length(dg);
indp=[1:n]';
% indl=[1:n]';
d=input('\n Insert the outer offset [e.g. 1,2 or 3 arc-degree, example: 3] = ');
% d=3;
% determination of the minimum and maximum latitudes and longitudes to
% define the boundary of the computation points in the inner grid
f2_in=max(fi)-d;
f1_in=min(fi)+d;
l2_in=max(lambda)-d;
l1_in=min(lambda)+d;
fprintf ('\n Insert the resolution of the gravity grid dphi and dlambda [in arc-minutes] : ');
dF=input('\n dphi    (example: 3) = ');
dL=input('\n dlambda (example: 3) = ');
% dF=ii;
% dL=ii;
dfi=dF/60;
dlam=dL/60;
d1=input('\n Insert the prefered spherical cap size (psi0) [0.25,0.5,...,3 arc-degree] = ');
% Excerption of the inner zone grid data from the main gravity grid 
% index by the following relation operators
inx=indp(f1_in<=fi & f2_in>=fi & l1_in<=lambda & l2_in>=lambda);
% Retrieving inner grid values by indexes
anomaly2=anomaly(inx,:);
% Loading the inner zone data
nm=length(anomaly2);
hh=zeros(nm,1);
    %####################################################################################
fprintf ('\n The following are the modification methods exist in this software : ');
fprintf ('\n 1- No modification');
fprintf ('\n 2- Wong and Gore');
fprintf ('\n 3- VK (Vanicek-Kleusberg)');
fprintf ('\n 4- LS biased');
fprintf ('\n 5- LS_unbiased');
fprintf ('\n 6- LS Optimum');
Sn=input('\nSelect method [1..6] = ');
    %####################################################################################

% switch Sn
% case 'LS biased'
if (Sn==1)
    filename=input('Enter the name of results file \n','s') ;
    out=fopen(filename,'wt');
    for i=1:nm   
        fm=anomaly2(i,2);
        lm=anomaly2(i,1);
        fi1=fm-fi;
        lambda1=lm-lambda;
        % Ex
        ind=indp(fi1<=d1 & fi1>=-d1 & lambda1<=d1 & lambda1>=-d1);
        f=fi1(ind,:);
        gravity2=anomaly(ind,:);
        f2=gravity2(:,2);
        l2=gravity2(:,1);
        % Calculating the spherical distance (psi) from equation 1-90 
        % (Hofmann-Wellenhof & Moritz, 2006) page 22
        psi=acos(sind(fm).*sind(f2)+cosd(fm).*cosd(f2).*cosd(lm-l2));
        m=length(psi);
        indX=[1:m]';
        % Specifying the truncated spherical cap around the computation point and
        % excluding the any identical points that are equal or close to the
        % computation point to avoid kernel singularity.
        indx=indX(psi<=(d*pi/180) & psi>1.0e-3);
        psi2=psi(indx);
        % Calling back the index values of the points within the truncation
        % boundary based on the index values of their positions
        c=gravity2(indx,:);
        ff=c(:,2);
        ll=c(:,1);
        dg2=c(:,3);
        gamma=somig(fm);
        [N1] =Stokes_func(fm,lm,dfi,dlam,c,R,gamma);
        fprintf('\n%2.8f %2.8f %2.4f',lm,fm,N1);
        fprintf(out,'\n%2.8f %2.8f %2.4f',lm,fm,N1);
end
fclose('all');

elseif (Sn==2)
    L=input('\n Insert the modification degree [1..360,example: 221] = ');
    filename=input('Enter the name of results file \n','s') ;
    out=fopen(filename,'wt');
    for i=1:nm
        fm=anomaly2(i,2);
        lm=anomaly2(i,1);
        fi1=fm-fi;
        lambda1=lm-lambda;
        % Ex
        ind=indp(fi1<=d1 & fi1>=-d1 & lambda1<=d1 & lambda1>=-d1);
        f=fi1(ind,:);
        gravity2=anomaly(ind,:);
        f2=gravity2(:,2);
        l2=gravity2(:,1);
        % Calculating the spherical distance (psi) from equation 1-90 
        % (Hofmann-Wellenhof & Moritz, 2006) page 22
        psi=acos(sind(fm).*sind(f2)+cosd(fm).*cosd(f2).*cosd(lm-l2));
        % psi=acosd(cosd(fi(indp(i),1)).*cosd(f2)+sind(fi(indp(i),1)).*sind(f2).*cosd(lambda(indp(i),1)-l2));
        m=length(psi);
        indX=[1:m]';
        % Specifying the truncated spherical cap around the computation point and
        % excluding the any identical points that are equal or close to the
        % computation point to avoid kernel singularity.
        indx=indX(psi<=(d*pi/180) & psi>1.0e-3);
        psi2=psi(indx);
        % Calling back the index values of the points within the truncation
        % boundary based on the index values of their positions
        c=gravity2(indx,:);
        ff=c(:,2);
        ll=c(:,1);
        dg2=c(:,3);
        gamma=somig(fm);
        [N1] =StkWGALL(fm,lm,dfi,dlam,c,R,gamma,L);
        fprintf('\n%2.8f %2.8f %2.4f',lm,fm,N1);
        fprintf(out,'\n%2.8f %2.8f %2.4f',lm,fm,N1);
    end
    fclose('all');

elseif (Sn==3)
    Snd=input('Insert the file of the modification Parameters (example: tk.prn) \n','s') ;
    Sn1=load(Snd);
    L=input('\n Insert the modification degree [1..360,example: 221] = ');
    filename=input('Enter the name of results file \n','s') ;
    out=fopen(filename,'wt');
    for i=1:nm
        fm=anomaly2(i,2);
        lm=anomaly2(i,1);
        fi1=fm-fi;
        lambda1=lm-lambda;
        % Ex
        ind=indp(fi1<=d1 & fi1>=-d1 & lambda1<=d1 & lambda1>=-d1);
        f=fi1(ind,:);
        gravity2=anomaly(ind,:);
        f2=gravity2(:,2);
        l2=gravity2(:,1);
         % Calculating the spherical distance (psi) from equation 1-90 
        % (Hofmann-Wellenhof & Moritz, 2006) page 22
        psi=acos(sind(fm).*sind(f2)+cosd(fm).*cosd(f2).*cosd(lm-l2));
        % psi=acosd(cosd(fi(indp(i),1)).*cosd(f2)+sind(fi(indp(i),1)).*sind(f2).*cosd(lambda(indp(i),1)-l2));
        m=length(psi);
        indX=[1:m]';
        % Specifying the truncated spherical cap around the computation point and
        % excluding the any identical points that are equal or close to the
        % computation point to avoid kernel singularity.
        indx=indX(psi<=(d*pi/180) & psi>1.0e-3);
        psi2=psi(indx);
        % Calling back the index values of the points within the truncation
        % boundary based on the index values of their positions
        c=gravity2(indx,:);
        ff=c(:,2);
        ll=c(:,1);
        dg2=c(:,3);
        gamma=somig(fm);
        [NVK] = StkVK(fm,lm,dfi,dlam,c,R,gamma,L,Sn1);
        fprintf('\n%2.8f %2.8f %2.4f',lm,fm,NVK);
        fprintf(out,'\n%2.8f %2.8f %2.4f',lm,fm,NVK);
    end
elseif (Sn==4)
    Snd=input('Insert the file of the modification Parameters (example: Sn_bias.prn) \n','s') ;
    Sn1=load(Snd);
    L=input('\n Insert the modification degree [1..360,example: 221] = ');
    filename=input('Enter the name of results file \n','s') ;
    out=fopen(filename,'wt');
    for i=1:nm
        fm=anomaly2(i,2);
        lm=anomaly2(i,1);
        fi1=fm-fi;
        lambda1=lm-lambda;
        % Ex
        ind=indp(fi1<=d1 & fi1>=-d1 & lambda1<=d1 & lambda1>=-d1);
        f=fi1(ind,:);
        gravity2=anomaly(ind,:);
        f2=gravity2(:,2);
        l2=gravity2(:,1);
         % Calculating the spherical distance (psi) from equation 1-90 
        % (Hofmann-Wellenhof & Moritz, 2006) page 22
        psi=acos(sind(fm).*sind(f2)+cosd(fm).*cosd(f2).*cosd(lm-l2));
        % psi=acosd(cosd(fi(indp(i),1)).*cosd(f2)+sind(fi(indp(i),1)).*sind(f2).*cosd(lambda(indp(i),1)-l2));
        m=length(psi);
        indX=[1:m]';
        % Specifying the truncated spherical cap around the computation point and
        % excluding the any identical points that are equal or close to the
        % computation point to avoid kernel singularity.
        indx=indX(psi<=(d*pi/180) & psi>1.0e-3);
        psi2=psi(indx);
        % Calling back the index values of the points within the truncation
        % boundary based on the index values of their positions
        c=gravity2(indx,:);
        ff=c(:,2);
        ll=c(:,1);
        dg2=c(:,3);
        gamma=somig(fm);
        [N1] = StkWLS(fm,lm,dfi,dlam,c,R,gamma,L,Sn1);
        fprintf('\n%2.8f %2.8f %2.4f',lm,fm,N1);
        fprintf(out,'\n%2.8f %2.8f %2.4f',lm,fm,N1);
    end
    fclose('all');

elseif (Sn==5)
    Snd=input('Insert the file of the modification Parameters (example: Bn_unb.prn) \n','s') ;
    Bn=load(Snd);
    L=input('\n Insert the modification degree [1..360, Example: 221] = ');
    filename=input('Enter the name of results file \n','s') ;
    out=fopen(filename,'wt');
    for i=1:nm
        fm=anomaly2(i,2);
        lm=anomaly2(i,1);
        fi1=fm-fi;
        lambda1=lm-lambda;
        % Ex
        ind=indp(fi1<=d1 & fi1>=-d1 & lambda1<=d1 & lambda1>=-d1);
        f=fi1(ind,:);
        gravity2=anomaly(ind,:);
        f2=gravity2(:,2);
        l2=gravity2(:,1);
         % Calculating the spherical distance (psi) from equation 1-90 
        % (Hofmann-Wellenhof & Moritz, 2006) page 22
        psi=acos(sind(fm).*sind(f2)+cosd(fm).*cosd(f2).*cosd(lm-l2));
        % psi=acosd(cosd(fi(indp(i),1)).*cosd(f2)+sind(fi(indp(i),1)).*sind(f2).*cosd(lambda(indp(i),1)-l2));
        m=length(psi);
        indX=[1:m]';
        % Specifying the truncated spherical cap around the computation point and
        % excluding the any identical points that are equal or close to the
        % computation point to avoid kernel singularity.
        indx=indX(psi<=(d*pi/180) & psi>1.0e-3);
        psi2=psi(indx);
        % Calling back the index values of the points within the truncation
        % boundary based on the index values of their positions
        c=gravity2(indx,:);
        ff=c(:,2);
        ll=c(:,1);
        dg2=c(:,3);
        gamma=somig(fm);
        [N1] = StkWLS(fm,lm,dfi,dlam,c,R,gamma,L,Bn);
        fprintf('\n%2.8f %2.8f %2.4f',lm,fm,N1);
        fprintf(out,'\n%2.8f %2.8f %2.4f',lm,fm,N1);
end
    fclose('all');

elseif (Sn==6)
    Snd=input('Insert the file of the modification Parameters (example: Bn_opt.prn) \n','s') ;
    Bn=load(Snd);
    L=input('\n Insert the modification degree [1..360, example: 221] = ');
    filename=input('Enter the name of results file \n','s') ;
    out=fopen(filename,'wt');
    for i=1:nm
        fm=anomaly2(i,2);
        lm=anomaly2(i,1);
        fi1=fm-fi;
        lambda1=lm-lambda;
        % Ex
        ind=indp(fi1<=d1 & fi1>=-d1 & lambda1<=d1 & lambda1>=-d1);
        f=fi1(ind,:);
        gravity2=anomaly(ind,:);
        f2=gravity2(:,2);
        l2=gravity2(:,1);
         % Calculating the spherical distance (psi) from equation 1-90 
        % (Hofmann-Wellenhof & Moritz, 2006) page 22
        psi=acos(sind(fm).*sind(f2)+cosd(fm).*cosd(f2).*cosd(lm-l2));
        % psi=acosd(cosd(fi(indp(i),1)).*cosd(f2)+sind(fi(indp(i),1)).*sind(f2).*cosd(lambda(indp(i),1)-l2));
        m=length(psi);
        indX=[1:m]';
        % Specifying the truncated spherical cap around the computation point and
        % excluding the any identical points that are equal or close to the
        % computation point to avoid kernel singularity.
        indx=indX(psi<=(d*pi/180) & psi>1.0e-3);
        psi2=psi(indx);
        % Calling back the index values of the points within the truncation
        % boundary based on the index values of their positions
        c=gravity2(indx,:);
        ff=c(:,2);
        ll=c(:,1);
        dg2=c(:,3);
        gamma=somig(fm);
        [N1] = StkWLS(fm,lm,dfi,dlam,c,R,gamma,L,Bn);
        fprintf('\n%2.8f %2.8f %2.4f',lm,fm,N1);
        fprintf(out,'\n%2.8f %2.8f %2.4f',lm,fm,N1);
end
    fclose('all');    
end


