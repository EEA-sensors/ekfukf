% Dynamical model function for reentry problem.
% Discretization is done using a simple Euler
% time integration.

% 
% Copyright (C) 2005-2006 Simo Särkkä
%               2007      Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function x = reentry_f(xw,param)

  dt  = param{1};
  b0  = param{2};
  H0  = param{3};
  Gm0 = param{4};
  R0  = param{5};
  
  x = xw(1:5,:);

  R = sqrt(x(1,:).^2 + x(2,:).^2);
  V = sqrt(x(3,:).^2 + x(4,:).^2);
  b = b0 * exp(x(5,:));
  D = b .* exp((R0-R)/H0) .* V;
  G = -Gm0 ./ R.^3;
  dot_x = zeros(size(x));
  dot_x(1,:) = x(3,:);
  dot_x(2,:) = x(4,:);
  dot_x(3,:) = D .* x(3,:) + G .* x(1,:);
  dot_x(4,:) = D .* x(4,:) + G .* x(2,:);
  dot_x(5,:) = zeros(1,size(x,2));

  % Euler integration
  x = x + dt * dot_x;

  % Add process noise if the state is augmented 
  if size(xw,1) > 5 && length(param) > 5
    L   = param{6};  
    w = xw(size(L,1)+1:size(L,1)+size(L,2),:);
    x = x + L * w;
  end
  
