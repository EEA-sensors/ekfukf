%KF_UPDATE  Kalman Filter update step
%
% Syntax:
%   [X,P,K,IM,IS,LH] = KF_UPDATE(X,P,Y,H,R)
%
% In:
%   X - Nx1 mean state estimate after prediction step
%   P - NxN state covariance after prediction step
%   Y - Dx1 measurement vector.
%   H - Measurement matrix.
%   R - Measurement noise covariance.
%
% Out:
%   X  - Updated state mean
%   P  - Updated state covariance
%   K  - Computed Kalman gain
%   IM - Mean of predictive distribution of Y
%   IS - Covariance or predictive mean of Y
%   LH - Predictive probability (likelihood) of measurement.
%   
% Description:
%   Kalman filter measurement update step. Kalman Filter
%   model is
%
%     x[k] = A*x[k-1] + B*u[k-1] + q,  q ~ N(0,Q)
%     y[k] = H*x[k]   + r,             r ~ N(0,R)
%
%   Prediction step of Kalman filter computes predicted
%   mean m-[k] and covariance P-[k] of state:
%
%     p(x[k] | y[1:k-1]) = N(x[k] | m-[k], P-[k])
%
%   See for instance KF_PREDICT how m-[k] and P-[k] are
%   calculated. 
%
%   Update step computes the posterior mean m[k] and
%   covariance P[k]  of state given new measurement:
%
%     p(x[k] | y[1:k]) = N(x[k] | m[k], P[k])
%
%   Innovation distribution is defined as
%
%     p(y[k] | y[1:k-1]) = N(y[k] | IM[k], IS[k])
%   
%   Updated mean x[k] and covarience P[k] are given by
%   the following equations (not the only possible ones):
%
%     v[k] = y[k] - H[k]*m-[k]
%     S[k] = H[k]*P-[k]*H[k]' + R[k]
%     K[k] = P-[k]*H[k]'*[S[k]]^(-1) 
%     m[k] = m-[k] + K[k]*v[k]
%     P[k] = P-[k] - K[k]*S[k]*K[k]'
%
% Example:
%   m = m0;
%   P = P0;
%   M = m0;
%   for i=1:size(Y,2)
%     [m,P] = kf_predict(m,P,A,Q);
%     [m,P] = kf_update(m,P,Y(:,i),H,R);
%     M = [M m];
%   end
%
% See also:
%   KF_PREDICT, EKF_UPDATE

% History:
%   26.02.2007 JH Added the equations for calculating the updated
%                 means and covariances to the description section.
%   12.01.2003 SS Symmetrized covariance update
%   20.11.2002 SS The first official version.
%
% Copyright (C) 2002, 2003 Simo Särkkä
% Copyright (C) 2007 Jouni Hartikainen
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [X,P,K,IM,IS,LH] = kf_update(X,P,y,H,R)

  %
  % Check which arguments are there
  %
  if nargin < 5
    error('Too few arguments');
  end

  %
  % update step
  %
  IM = H*X;
  IS = (R + H*P*H');
  K = P*H'/IS;
  X = X + K * (y-IM);
  P = P - K*IS*K';
  if nargout > 5
    LH = gauss_pdf(y,IM,IS);
  end
