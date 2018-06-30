%EKF_PREDICT1  1st order Extended Kalman Filter prediction step
%
% Syntax:
%   [M,P] = EKF_PREDICT1(M,P,[A,Q,a,W,param])
%
% In:
%   M - Nx1 mean state estimate of previous step
%   P - NxN state covariance of previous step
%   A - Derivative of a() with respect to state as
%       matrix, inline function, function handle or
%       name of function in form A(x,param)       (optional, default eye())
%   Q - Process noise of discrete model               (optional, default zero)
%   a - Mean prediction E[a(x[k-1],q=0)] as vector,
%       inline function, function handle or name
%       of function in form a(x,param)                (optional, default A(x)*X)
%   W - Derivative of a() with respect to noise q
%       as matrix, inline function, function handle
%       or name of function in form W(x,param)        (optional, default identity)
%   param - Parameters of a                           (optional, default empty)
%
% Out:
%   M - Updated state mean
%   P - Updated state covariance
%   
% Description:
%   Perform Extended Kalman Filter prediction step.
%
% See also:
%   EKF_UPDATE1, EKF_PREDICT2, EKF_UPDATE2, DER_CHECK,
%   LTI_DISC, KF_PREDICT, KF_UPDATE

% Copyright (C) 2002-2006 Simo Särkkä
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [M,P] = ekf_predict1(M,P,A,Q,a,W,param)

  %
  % Check arguments
  %
  if nargin < 3
    A = [];
  end
  if nargin < 4
    Q = [];
  end
   if nargin < 5
    a = [];
  end
  if nargin < 6
    W = [];
  end
  if nargin < 7
    param = [];
  end
  
  %
  % Apply defaults
  %
  if isempty(A)
    A = eye(size(M,1));
  end
  if isempty(Q)
    Q = zeros(size(M,1));
  end
  if isempty(W)
    W = eye(size(M,1),size(Q,2));
  end

  if isnumeric(A)
    % nop
  elseif ischar(A) | strcmp(class(A),'function_handle')
    A = feval(A,M,param);
  else
    A = A(M,param);
  end
  %
  % Perform prediction
  %

  if isempty(a)
    M = A*M;
  elseif isnumeric(a)
    M = a;
  elseif ischar(a) | strcmp(class(a),'function_handle')
    M = feval(a,M,param);
  else
    M = a(M,param);
  end



  if isnumeric(W)
    % nop
  elseif ischar(W) | strcmp(class(W),'function_handle')
    W = feval(W,M,param);
  else
    W = W(M,param);
  end

  P = A * P * A' + W * Q * W';
