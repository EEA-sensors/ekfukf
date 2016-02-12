% Hessian of the measurement function in BOT-demo.

% Copyright (C) 2007 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function dY = bot_d2h_dx2(x,s)
  % Space for Hessians. Note that we need a Hessian for
  % each dimension in the measurement space, that is we need
  % a Hessian for each sensor in this case.     
  dY = zeros(size(s,2),size(x,1),size(x,1));
  
  % Loop through sensors.
  for i=1:size(s,2)
    % Derivative twice wrt. x
    dx2 = -2*(x(1)-s(1,i)) / ((x(1)-s(1,i))^2+(x(2)-s(2,i))^2)^2;
    % Derivative twice wrt. y    
    dy2 = -2*(x(2)-s(2,i)) / ((x(1)-s(1,i))^2+(x(2)-s(2,i))^2)^2;
    % Derivative wrt. x and y
    dxdy = ((x(2)-s(2,i))^2-(x(1)-s(1,i))^2) / ((x(1)-s(1,i))^2+(x(2)-s(2,i))^2)^2;
    dh = [dx2  dxdy 0 0;...
	  dxdy dy2  0 0;...
	  0    0    0 0;...
          0    0    0 0];
    dY(i,:,:) = dh;
  end