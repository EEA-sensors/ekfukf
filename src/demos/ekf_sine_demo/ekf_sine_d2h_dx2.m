% Hessian of the measurement model function in the random sine signal demo

% Copyright (C) 2007 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function df = ekf_sine_d2h_dx2(x,param)
f = x(1);
a = x(3);

df = zeros(1,3,3);
df(1,:,:) = [-a*sin(f) 0 cos(f);
             0         0      0;
             cos(f)    0      0];
