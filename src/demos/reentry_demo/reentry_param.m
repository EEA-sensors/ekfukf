% The parameters for re-entry mechanics.
 
% Copyright (C) 2005-2006 Simo Särkkä
%               2007      Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.
  
b0 = -0.59783;
H0 = 13.406;
Gm0 = 3.9860e5;
R0 = 6374;

dt = 0.1;
sim_iter = 1;

xr = R0;
yr = 0;

vr = (1e-3).^2;
va = (.17e-3).^2;


true_m0 = [6500.4;
           349.14;
           -1.8093;
           -6.7967;
           0.6932];

true_P0 = diag([1e-6 1e-6 1e-6 1e-6 0]);

true_Qc = diag([2.4064e-5 2.4064e-5 0]) / 0.1;

L = [0 0 0;
     0 0 0;
     1 0 0;
     0 1 0;
     0 0 1];

m0 = [6500.4;
      349.14;
      -1.8093;
      -6.7967;
      0];

P0 = diag([1e-6 1e-6 1e-6 1e-6 1]);

Qc = diag([2.4064e-5 2.4064e-5 1e-6]) / 0.1;

% Parameters needed by the dynamic model as cell array
d_param = {dt,b0,H0,Gm0,R0,L};

% Parameters needed by the measurement model as cell array
h_param = {xr,yr};