% Jacobian of the state transition function in reentry demo.

% 
% Copyright (C) 2005-2006 Simo Särkkä
%               2007      Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function da = reentry_df_dx(x,param)

  dt  = param{1};
  b0  = param{2};
  H0  = param{3};
  Gm0 = param{4};
  R0  = param{5};

  R = sqrt(x(1,:).^2 + x(2,:).^2);
  V = sqrt(x(3,:).^2 + x(4,:).^2);
  b = b0 * exp(x(5,:));
  D = b .* exp((R0-R)/H0) * V;
  G = -Gm0 ./ R.^3;

  dR_dx1 = x(1) / R;
  dR_dx2 = x(2) / R;
  dV_dx3 = x(3) / V;
  dV_dx4 = x(4) / V;
  db_dx5 = b;
  dD_dx1 = b * (-dR_dx1/H0) * exp((R0-R)/H0) * V;
  dD_dx2 = b * (-dR_dx2/H0) * exp((R0-R)/H0) * V;
  dD_dx3 = b * exp((R0-R)/H0) * dV_dx3;
  dD_dx4 = b * exp((R0-R)/H0) * dV_dx4;
  dD_dx5 = db_dx5 * exp((R0-R)/H0) * V;
  dG_dx1 = -Gm0 * (-3 * dR_dx1 / R^4);
  dG_dx2 = -Gm0 * (-3 * dR_dx2 / R^4);

  df = zeros(5,5);
  df(1,3) = 1;
  df(2,4) = 1;
  df(3,1) = dD_dx1 * x(3) + dG_dx1 * x(1) + G;
  df(3,2) = dD_dx2 * x(3) + dG_dx2 * x(1);
  df(3,3) = dD_dx3 * x(3) + D;
  df(3,4) = dD_dx4 * x(3);
  df(3,5) = dD_dx5 * x(3);
  
  df(4,1) = dD_dx1 * x(4) + dG_dx1 * x(2); 
  df(4,2) = dD_dx2 * x(4) + dG_dx2 * x(2) + G;
  df(4,3) = dD_dx3 * x(4);
  df(4,4) = dD_dx4 * x(4) + D;
  df(4,5) = dD_dx5 * x(4);

  da = eye(size(df,1)) + dt * df;
  
