%EKF_UPDATE2  2nd order Extended Kalman Filter update step
%
% Syntax:
%   [M,P,K,MU,S,LH] = EKF_UPDATE2(M,P,Y,H,H_xx,R,[h,V,param])
%
% In:
%   M  - Nx1 mean state estimate after prediction step
%   P  - NxN state covariance after prediction step
%   Y  - Dx1 measurement vector.
%   H  - Derivative of h() with respect to state as matrix,
%        inline function, function handle or name
%        of function in form H(x,param)
%   H_xx - DxNxN Hessian of h() with respect to state as matrix,
%          inline function, function handle or name of function
%          in form H_xx(x,param) 
%   R  - Measurement noise covariance.
%   h  - Mean prediction (measurement model) as vector,
%        inline function, function handle or name
%        of function in form h(x,param).                 (optional, default H(x)*X)
%   V  - Derivative of h() with respect to noise as matrix,
%        inline function, function handle or name
%        of function in form V(x,param).                 (optional, default identity)
%   param - Parameters of h                              (optional, default empty)
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
%   Extended Kalman Filter measurement update step.
%   EKF model is
%
%     y[k] = h(x[k],r),   r ~ N(0,R)
%
% See also:
%   EKF_PREDICT1, EKF_UPDATE1, EKF_PREDICT2, DER_CHECK, LTI_DISC, 
%   KF_UPDATE, KF_PREDICT

% Copyright (C) 2002-2006 Simo Särkkä
% Copyright (C) 2007 Jouni Hartikainen
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [M,P,K,IM,S,LH] = ekf_update2(M,P,y,H,H_xx,R,h,V,param)

  %
  % Check which arguments are there
  %
  if nargin < 6
    error('Too few arguments');
  end
  if nargin < 7
    h = [];
  end
  if nargin < 8
    V = [];
  end
  if nargin < 9
    param = [];
  end

  %
  % Apply defaults
  %
  if isempty(V)
    V = eye(size(R,1));
  end

  %
  % Evaluate matrices
  %
  if isnumeric(H)
    % nop
  elseif ischar(H) | strcmp(class(H),'function_handle')
    H = feval(H,M,param);
  else
    H = H(M,param);
  end
  
  if isnumeric(H_xx)
    % nop
  elseif ischar(H_xx) | strcmp(class(H_xx),'function_handle')
    H_xx = feval(H_xx,M,param);
  else
    H_xx = H_xx(M,param);
  end

  if isempty(h)
    MU = H*M;
  elseif isnumeric(h)
    MU = h;
  elseif ischar(h) | strcmp(class(h),'function_handle')
    MU = feval(h,M,param);
  else
    MU = h(M,param);
  end

  if isnumeric(V)
    % nop
  elseif ischar(V) | strcmp(class(V),'function_handle')
    V = feval(V,M,param);
  else
    V = V(M,param);
  end

  
  %
  % update step
  %
  v = y - MU;
  for i = 1:size(H_xx,1)
    H_i = squeeze(H_xx(i,:,:));
    v(i) = v(i) - 0.5*trace(H_i*P);
  end
  
  S = (V*R*V' + H*P*H');
  for i = 1:size(H_xx,1)
    for j = 1:size(H_xx,1)
      H_i = squeeze(H_xx(i,:,:));
      H_j = squeeze(H_xx(j,:,:));
      S(i,j) = S(i,j) + 0.5*trace(H_i*P*H_j*P);
    end
  end
  K = P*H'/S;
  M = M + K * v;
  P = P - K*S*K';

  if nargout > 5
    LH = gauss_pdf(y,MU,S);
  end
