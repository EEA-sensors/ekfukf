%AZ_H Azimuth measurement function for EKF/UKF.
%
%  h = atan((y-sy) / (x-sx))
%

% Copyright (C) 2003 Simo Särkkä
%
% $Id: az_h.m,v 1.1 2006/09/26 21:32:32 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function Y = az_h(x,s)
  Y = zeros(size(s,2),size(x,2));

  for i=1:size(s,2)
    h = atan2(x(2,:)-s(2,i),x(1,:)-s(1,i));
    np = find(h>0.5*pi);
    nn = find(h<-0.5*pi);
    if length(nn)>length(np)
      h(np) = h(np)-2*pi;
    else
      h(nn) = h(nn)+2*pi;
    end
    Y(i,:) = h;
  end


