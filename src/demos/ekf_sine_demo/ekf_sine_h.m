% Measurement model function for the random sine signal demo

% Copyright (C) 2007 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function Y = ekf_sine_h(x,param)
   f = x(1,:);
   a = x(3,:);

   Y = a.*sin(f);
   if size(x,1) == 7
      Y = Y + x(7,:);
   end
