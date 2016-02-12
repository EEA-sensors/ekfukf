% Hessian of the measurement model function in UNGM-model.

% Copyright (C) 2007 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function d = ungm_d2h_dx2(x,param)
    d = zeros(1,1,1);
    d(1,1,1) = 1/10;