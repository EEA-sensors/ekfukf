% Jacobian of the measurement model function in BOT demo.
%
%  dh_dx = -(y-sy) / (x-sx)^2 * 1 / (1 + (y-sy)^2 / (x-sx)^2)
%        = -(y-sy) / ((x-sx)^2 + (y-sy)^2)
%  dh_dy = 1 / (x-sx) * 1 / (1 + (y-sy)^2 / (x-sx)^2)
%        = (x-sx) / ((x-sx)^2 + (y-sy)^2)
%

% Copyright (C) 2003 Simo Särkkä
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function dY = bot_dh_dx(x,s)
  % Reserve space for the Jacobian. 
  dY = zeros(size(s,2),size(x,1));
  
  % Loop through sensors
  for i=1:size(s,2)
    dh = [-(x(2)-s(2,i)) / ((x(1)-s(1,i))^2 + (x(2)-s(2,i))^2);...
	   (x(1)-s(1,i)) / ((x(1)-s(1,i))^2 + (x(2)-s(2,i))^2);...
	  zeros(size(x,1)-2,1)]';
    dY(i,:) = dh; 
  end

