
function [M,P,K,MU,S,LH] = ckf_update(M,P,Y,h,R,h_param)
% CKF_UPDATE - Cubature Kalman filter update step
%
% Syntax:
%   [M,P,K,MU,S,LH] = CKF_UPDATE(M,P,Y,h,R,param)
%
% In:
%   M  - Mean state estimate after prediction step
%   P  - State covariance after prediction step
%   Y  - Measurement vector.
%   h  - Measurement model function as a matrix H defining
%        linear function h(x) = H*x, inline function,
%        function handle or name of function in
%        form h(x,param)
%   R  - Measurement covariance.
%   h_param - Parameters of h.
%
% Out:
%   M  - Updated state mean
%   P  - Updated state covariance
%   K  - Computed Kalman gain
%   MU - Predictive mean of Y
%   S  - Predictive covariance Y
%   LH - Predictive probability (likelihood) of measurement.
%   
% Description:
%   Perform additive form spherical-radial cubature Kalman filter (CKF)
%   measurement update step. Assumes additive measurement noise.
%
%   Function h should be such that it can be given
%   DxN matrix of N sigma Dx1 points and it returns 
%   the corresponding measurements for each sigma
%   point. This function should also make sure that
%   the returned sigma points are compatible such that
%   there are no 2pi jumps in angles etc.
%
% Example:
%   h = inline('atan2(x(2,:)-s(2),x(1,:)-s(1))','x','s');
%   [M2,P2] = ckf_update(M1,P1,Y,h,R,S);
%
% See also:
%   CKF_PREDICT, CRTS_SMOOTH, CKF_TRANSFORM, SPHERICALRADIAL
% 
% References:
%   Arasaratnam and Haykin (2009). Cubature Kalman Filters.
%    IEEE Transactions on Automatic Control, vol. 54, no. 5, pp.1254-1269

% Copyright (c) 2010 Arno Solin
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.
%%

  %
  % Check that all arguments are there
  %
  if nargin < 5
    error('Too few arguments');
  end


  %
  % Do transform and make the update
  %
  if nargin < 6
    [MU,S,C,X] = ckf_transform(M,P,h);
  else  
    [MU,S,C,X] = ckf_transform(M,P,h,h_param);
  end
  
  S = S + R;
  K = C / S;
  M = M + K * (Y - MU);
  P = P - K * S * K';
  
  if nargout > 5
    LH = gauss_pdf(Y,MU,S);
  end
  
 