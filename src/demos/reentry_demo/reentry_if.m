% Inverse dynamical function for reentry demo.

% 
% Copyright (C) 2005-2006 Simo Särkkä
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function x = reentry_if(x,param)

  y = reentry_f(x,param);
  x = 2 * x(1:size(y,1)) - y;

