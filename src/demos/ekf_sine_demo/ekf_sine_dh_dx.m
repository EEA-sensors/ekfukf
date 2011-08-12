% Jacobian of the measurement model function in the random sine signal demo

% Copyright (C) 2007 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function dY = ekf_sine_dh_dx(x,param)
f = x(1,:);
w = x(2,:);
a = x(3,:);
  
dY = [(a.*cos(f))' zeros(size(f,2),1) (sin(f))'];