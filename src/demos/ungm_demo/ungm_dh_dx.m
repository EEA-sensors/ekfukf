% Jacobian of the measurement model function for the UNGM-model.
%
% Copyright (C) 2007 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function dh = ungm_dh_dx(x,param)
dh = x/10;