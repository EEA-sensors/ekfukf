%EKF_PREDICT2  2nd order Extended Kalman Filter prediction step
%
% Syntax:
%   [M,P] = EKF_PREDICT2(M,P,[A,F,Q,a,W,param])
%
% In:
%   M - Nx1 mean state estimate of previous step
%   P - NxN state covariance of previous step
%   A - Derivative of a() with respect to state as
%       matrix, inline function, function handle or
%       name of function in form A(x,param)                 (optional, default identity)
%   F - NxNxN Hessian matrix of the state transition function
%       w.r.t. state variables as matrix, inline
%       function, function handle or name of function
%       in form F(x,param)                                  (optional, default identity)
%   Q - Process noise of discrete model                     (optional, default zero)
%   a - Mean prediction E[a(x[k-1],q=0)] as vector,
%       inline function, function handle or name
%       of function in form a(x,param)                      (optional, default A(x)*X)
%   W - Derivative of a() with respect to noise q
%       as matrix, inline function, function handle
%       or name of function in form W(x,k-1,param)          (optional, default identity)
%   param - Parameters of a                                 (optional, default empty)
%
%   
%
% Out:
%   M - Updated state mean
%   P - Updated state covariance
%   
% Description:
%   Perform Extended Kalman Filter prediction step.
%
% See also:
%   EKF_PREDICT1, EKF_UPDATE1, EKF_UPDATE2, DER_CHECK, LTI_DISC, 
%   KF_PREDICT, KF_UPDATE

% History:
%   22.5.07  JH Initial version. Modified from ekf_predict1.m
%            originally created by SS.
%
% Copyright (C) 2007 Jouni Hartikainen, Simo S�rkk�
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.
function [M,P] = ekf_predict2(M,P,A,F,Q,a,W,param)

  %
  % Check arguments
  %
  if nargin < 3
    A = [];
  end
  if nargin < 4
    F = []
  end
  if nargin < 5
    Q = [];
  end
   if nargin < 6
    a = [];
  end
  if nargin < 7
    W = [];
  end
  if nargin < 8
    param = [];
  end
  
  % Apply defaults
  if isempty(A)
    A = eye(size(M,1));
  end
  if isempty(F)
    F = eye(size(M,1));
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

  if isnumeric(F)
    % nop
  elseif ischar(F) | strcmp(class(F),'function_handle')
    F = feval(F,M,param);
  else
    F = F(M,param);
  end
  
  % Perform prediction
  if isempty(a)
    M = A*M;
  elseif isnumeric(a)
    M = a;
  elseif ischar(a) | strcmp(class(a),'function_handle')
    M = feval(a,M,param);
  else
    M = a(M,param);
  end
  
  for i=1:size(F,1)
    F_i = squeeze(F(i,:,:));  
    M(i) = M(i) + 0.5*trace(F_i*P);
  end


  if isnumeric(W)
    % nop
  elseif ischar(W) | strcmp(class(W),'function_handle')
    W = feval(W,M,param);
  else
    W = W(M,param);
  end

  P_new = A * P * A' + W * Q * W';
  for i = 1:size(F,1)
    for j = 1:size(F,1)
      F_i = squeeze(F(i,:,:));
      F_j = squeeze(F(j,:,:));
      P_new(i,j) = P_new(i,j)+0.5*trace(F_i*P*F_j*P);
    end
  end
  P = P_new;