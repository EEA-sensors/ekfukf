% Measurement model function for the UNGM-model.
%
% Copyright (C) 2007 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function y_n = ungm_h(x_n,param)
y_n = x_n(1,:).*x_n(1,:) ./ 20;
if size(x_n,1) == 3
   y_n = y_n + x_n(3,:);
end
