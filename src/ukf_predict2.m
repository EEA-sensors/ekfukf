%UKF_PREDICT2  Augmented form Unscented Kalman Filter prediction step
%
% Syntax:
%   [M,P] = UKF_PREDICT2(M,P,a,Q,[param,alpha,beta,kappa])
%
% In:
%   M - Nx1 mean state estimate of previous step
%   P - NxN state covariance of previous step
%   a - Dynamic model function as inline function,
%       function handle or name of function in
%       form a([x;w],param)
%   Q - Non-singular covariance of process noise w
%   param - Parameters of a               (optional, default empty)
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

% Copyright (C) 2003-2006 Simo Särkkä
%
% $Id: ukf_predict2.m,v 1.1 2006/09/29 13:43:35 ssarkka Exp $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function [M,P] = ukf_predict2(M,P,a,Q,param,alpha,beta,kappa,mat)

  %
  % Check which arguments are there
  %
  if nargin < 2
    error('Too few arguments');
  end
  if nargin < 3
    a = [];
  end
  if nargin < 4
    Q = [];
  end
  if nargin < 5
    param = [];
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
  [M,P] = ut_transform(MA,PA,a,param,alpha,beta,kappa,mat);
