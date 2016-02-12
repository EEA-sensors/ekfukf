%UKF_UPDATE1 -  Additive form Unscented Kalman Filter update step
%
% Syntax:
%   [M,P,K,MU,S,LH] = UKF_UPDATE1(M,P,Y,h,R,param,alpha,beta,kappa,mat)
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
%   h_param - Parameters of h               (optional, default empty)
%   alpha - Transformation parameter      (optional)
%   beta  - Transformation parameter      (optional)
%   kappa - Transformation parameter      (optional)
%   mat   - If 1 uses matrix form         (optional, default 0)
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
%   Perform additive form Discrete Unscented Kalman Filter (UKF)
%   measurement update step. Assumes additive measurement
%   noise.
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
%   [M2,P2] = ukf_update(M1,P1,Y,h,R,S);
%
% See also:
%   UKF_PREDICT1, UKF_PREDICT2, UKF_UPDATE2, UKF_PREDICT3, UKF_UPDATE3,
%   UT_TRANSFORM, UT_WEIGHTS, UT_MWEIGHTS, UT_SIGMAS
%
% History:
%   08.02.2008 JH Fixed a typo in the syntax description. 
%   04.05.2007 JH Made corrections to the description.
%   02.05.2007 JH Fixed a bug in likelihood calculation and added
%              a "See also"-section. 
%   2002-2006  SS Initial version
%   
% 
% References:
%   [1] Wan, Merwe: The Unscented Kalman Filter

% Copyright (C) 2002-2006 Simo S�rkk�
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [M,P,K,MU,S,LH] = ukf_update1(M,P,Y,h,R,h_param,alpha,beta,kappa,mat)

  %
  % Check that all arguments are there
  %
  if nargin < 5
    error('Too few arguments');
  end
  if nargin < 6
    h_param = [];
  end
  if nargin < 7
    alpha = [];
  end
  if nargin < 8
    beta = [];
  end
  if nargin < 9
    kappa = [];
  end
  if nargin < 10
    mat = [];
  end

  %
  % Apply defaults
  %
  if isempty(mat)
    mat = 0;
  end

  %
  % Do transform and make the update
  %
  tr_param = {alpha beta kappa mat};
  [MU,S,C] = ut_transform(M,P,h,h_param,tr_param);
  
  S = S + R;
  K = C / S;
  M = M + K * (Y - MU);
  P = P - K * S * K';
  
  if nargout > 5
    LH = gauss_pdf(Y,MU,S);
  end
