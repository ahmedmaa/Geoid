%**************************************************************************
%****************************** Geoid Software ****************************
%**************************************************************************
%
% This program script uses Stokes's Kernel to evaluate the terrestrial
% gravity measurements and then uses Stokeks kernel to compute the geoid
% model
%
%                                Ahmed Abdalla
%                        Assistant Professor Research
%                         Louisiana State University
%                          Center for Geoinformatics
%                             aabdalla1@lsu.edu
%                                 May 2019
%                          
%
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
GEOWARE

