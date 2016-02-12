%KF_LHOOD  Kalman Filter measurement likelihood
%
% Syntax:
%   LH = KF_LHOOD(X,P,Y,H,R)
%
% In:
%   X - Nx1 state mean
%   P - NxN state covariance
%   Y - Dx1 measurement vector.
%   H - Measurement matrix.
%   R - Measurement noise covariance.
%
% Out:
%   LH - Likelihood of measurement.
%   
% Description:
%   Calculate likelihood of measurement in Kalman filter.
%   If and X and P define the parameters of predictive
%   distribution (e.g. from KF_PREDICT)
%
%     p(x[k] | y[1:k-1]) = N(x[k] | m-[k], P-[k])
%
%   then this likelihood is the probability of measurement
%   in innovation distribution:
%
%     p(y[k] | y[1:k-1]) = N(y[k] | IM, IS)
%
% See also:
%   KF_PREDICT, KF_UPDATE

% History:
%   20.11.2002  The first official version.
%
% Copyright (C) 2002 Simo Särkkä
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function LH = kf_lhood(m,P,y,H,R)

  %
  % Check which arguments are there
  %
  if nargin < 5
    error('Too few arguments');
  end

  IM = H*m;
  IS = (H*P*H'+R);
  LH = gauss_pdf(y,IM,IS);
