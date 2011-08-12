% Dynamical model function for the random sine signal demo

% Copyright (C) 2007 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function x_n = ekf_sine_f(x,param)
    dt = param(1);
    A = [1 dt 0;0 1 0;0 0 1];
    x_n = A*x(1:3,:);
    if size(x,1) == 7 || size(x,1) == 6
        x_n(1:3,:) = x_n(1:3,:) + x(4:6,:);
    end
    