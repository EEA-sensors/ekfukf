%UKF_PREDICT2  Augmented (state and process noise) UKF prediction step
%
% Syntax:
%   [M,P] = UKF_PREDICT2(M,P,a,Q,[param,alpha,beta,kappa])
%
% In:
%   M - Nx1 mean state estimate of previous step
%   P - NxN state covariance of previous step
%   f - Dynamic model function as inline function,
%       function handle or name of function in
%       form a([x;w],param)
%   Q - Non-singular covariance of process noise w
%   f_param - Parameters of f               (optional, default empty)
%   alpha - Transformation parameter      (optional)
%   beta  - Transformation parameter      (optional)
%   kappa - Transformation parameter      (optional)
%   mat   - If 1 uses matrix form         (optional, default 0)
%
% Out:
%   M - Updated state mean
%   P - Updated state covariance
%
% Description:
%   Perform augmented form Unscented Kalman Filter prediction step
%   for model
%
%    x[k+1] = a(x[k],w[k],param)
%
%   Function a should be such that it can be given
%   DxN matrix of N sigma Dx1 points and it returns 
%   the corresponding predictions for each sigma
%   point. 
%
% See also:
%   UKF_PREDICT1, UKF_UPDATE1, UKF_UPDATE2, UKF_PREDICT3, UKF_UPDATE3,
%   UT_TRANSFORM, UT_WEIGHTS, UT_MWEIGHTS, UT_SIGMAS

% Copyright (C) 2003-2006 Simo S�rkk�
%
% $Id$
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function [M,P] = ukf_predict2(M,P,f,Q,f_param,alpha,beta,kappa,mat)

  %
  % Check which arguments are there
  %
  if nargin < 2
    error('Too few arguments');
  end
  if nargin < 3
    f = [];
  end
  if nargin < 4
    Q = [];
  end
  if nargin < 5
    f_param = [];
  end
  if nargin < 6
    alpha = [];
  end
  if nargin < 7
    beta = [];
  end
  if nargin < 8
    kappa = [];
  end
  if nargin < 9
    mat = [];
  end

  %
  % Apply defaults
  %
  if isempty(mat)
    mat = 0;
  end

  %
  % Do transform
  % and add process noise
  %
  m = size(M,1);
  n = size(Q,1);
  MA = [M;zeros(size(Q,1),1)];
  PA = zeros(size(P,1)+size(Q,1));
  PA(1:size(P,1),1:size(P,1)) = P;
  PA(1+size(P,1):end,1+size(P,1):end) = Q;
  tr_param = {alpha beta kappa mat};
  [M,P] = ut_transform(MA,PA,f,f_param,tr_param);
