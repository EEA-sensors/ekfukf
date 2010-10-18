%UKF_UPDATE2 - Augmented form Unscented Kalman Filter update step
%
% Syntax:
%   [M,P,K,MU,IS,LH] = UKF_UPDATE3(M,P,Y,h,R,X,w,h_param,alpha,beta,kappa,mat,sigmas)
%
% In:
%   M  - Mean state estimate after prediction step
%   P  - State covariance after prediction step
%   Y  - Measurement vector.
%   h  - Measurement model function as a matrix H defining
%        linear function h(x) = H*x+r, inline function,
%        function handle or name of function in
%        form h([x;r],param)
%   R  - Measurement covariance.
%   X - Sigma points of x
%   w - Weights as cell array {mean-weights,cov-weights,c}
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
%   Perform augmented form Discrete Unscented Kalman Filter (UKF)
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
%   [M2,P2] = ukf_update2(M1,P1,Y,h,R,S);
%
% See also:
%   UKF_PREDICT1, UKF_UPDATE1, UKF_PREDICT2, UKF_UPDATE2, UKF_PREDICT3
%   UT_TRANSFORM, UT_WEIGHTS, UT_MWEIGHTS, UT_SIGMAS
%

% History:
%   08.02.2008 JH Fixed a typo in the syntax description. 
%   04.05.2007 JH Initial version. Modified from ukf_update1.m
%              originally created by SS.
%   
% 
% References:
%   [1] Wan, Merwe: The Unscented Kalman Filter
%
% Copyright (C) 2007 Jouni Hartikainen, Simo S�rkk�
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [M,P,K,MU,S,LH] = ukf_update3(M,P,Y,h,R,X,w,h_param,alpha,beta,kappa,mat)

  %
  % Check that all arguments are there
  %
  if nargin < 5
    error('Too few arguments');
  end
  if nargin < 8
    h_param = [];
  end
  if nargin < 9
    alpha = [];
  end
  if nargin < 10
    beta = [];
  end
  if nargin < 11
    kappa = [];
  end
  if nargin < 12
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
  tr_param = {alpha beta kappa mat X w};
  [MU,S,C,X,Y_s] = ut_transform(M,P,h,h_param,tr_param);
  K = C / S;
  M = M + K * (Y - MU);
  P = P - K * S * K';
  if nargout > 5
    LH = gauss_pdf(Y,MU,S);
  end