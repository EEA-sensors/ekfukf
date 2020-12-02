%ERTS_SMOOTH1  Extended Rauch-Tung-Striebel smoother
%
% Syntax:
%   [M,P,D] = ERTS_SMOOTH1(M,P,A,Q,[a,W,param,same_p])
%
% In:
%   M - NxK matrix of K mean estimates from Unscented Kalman filter
%   P - NxNxK matrix of K state covariances from Unscented Kalman Filter
%   A - Derivative of a() with respect to state as
%       matrix, inline function, function handle or
%       name of function in form A(x,param)                 (optional, default eye())
%   Q - Process noise of discrete model                       (optional, default zero)
%   a - Mean prediction E[a(x[k-1],q=0)] as vector,
%       inline function, function handle or name
%       of function in form a(x,param)                        (optional, default A(x)*X)
%   W - Derivative of a() with respect to noise q
%       as matrix, inline function, function handle
%       or name of function in form W(x,param)                (optional, default identity)
%   param - Parameters of a. Parameters should be a single cell array, vector or a matrix
%           containing the same parameters for each step or if different parameters
%           are used on each step they must be a cell array of the format
%           { param_1, param_2, ...}, where param_x contains the parameters for
%           step x as a cell array, a vector or a matrix.     (optional, default empty)
%   same_p - 1 if the same parameters should be
%            used on every time step                          (optional, default 1)
%                                   
%                         
%
% Out:
%   K - Smoothed state mean sequence
%   P - Smoothed state covariance sequence
%   D - Smoother gain sequence
%   
% Description:
%   Extended Rauch-Tung-Striebel smoother algorithm. Calculate
%   "smoothed" sequence from given Kalman filter output sequence by
%   conditioning all steps to all measurements.
%
% Example:
%
% See also:
%   EKF_PREDICT1, EKF_UPDATE1

% History:
%   04.05.2007 JH Added the possibility to pass different parameters for a and h
%                 for each step.
%   2006       SS Initial version.  
%
% Copyright (C) 2006 Simo Särkkä
% Copyright (C) 2007 Jouni Hartikainen
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [M,P,D] = erts_smooth1(M,P,A,Q,a,W,param,same_p)

  %
  % Check which arguments are there
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
  if nargin < 8
    same_p = [];
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
  if isempty(same_p)
    same_p = 1;
  end

  %
  % Extend Q if NxN matrix
  %
  if size(Q,3)==1
    Q = repmat(Q,[1 1 size(M,2)]);
  end

  %
  % Run the smoother
  %
  D = zeros(size(M,1),size(M,1),size(M,2));
  for k=(size(M,2)-1):-1:1
    if isempty(param)
      params = [];  
    elseif same_p
        params = param;
    else
        params = param{k};
    end

    %
    % Perform prediction
    %
    if isempty(a)
      m_pred = A*M(:,k);
    elseif isnumeric(a)
      m_pred = a;
    elseif ischar(a) | strcmp(class(a),'function_handle')
      m_pred = feval(a,M(:,k),params);
    else
      m_pred = a(M(:,k),params);
    end
    
    if isnumeric(A)
      F = A;
    elseif ischar(A) | strcmp(class(A),'function_handle')
      F = feval(A,M(:,k),params);
    else
      F = A(M(:,k),params);
    end
    
    if isnumeric(W)
      B = W;
    elseif ischar(W) | strcmp(class(W),'function_handle')
      B = feval(W,M(:,k),params);
    else
      B = W(M(:,k),params);
    end

    P_pred = F * P(:,:,k) * F' + B * Q(:,:,k) * B';
    C = P(:,:,k) * F';
    
    D(:,:,k) = C / P_pred;
    M(:,k)   = M(:,k) + D(:,:,k) * (M(:,k+1) - m_pred);
    P(:,:,k) = P(:,:,k) + D(:,:,k) * (P(:,:,k+1) - P_pred) * D(:,:,k)';
  end
