%URTS_SMOOTH1  Additive form Unscented Rauch-Tung-Striebel smoother
%
% Syntax:
%   [M,P,D] = URTS_SMOOTH1(M,P,a,Q,[param,alpha,beta,kappa,mat,same_p])
%
% In:
%   M - NxK matrix of K mean estimates from Unscented Kalman filter
%   P - NxNxK matrix of K state covariances from Unscented Kalman Filter
%   a - Dynamic model function as a matrix A defining
%       linear function a(x) = A*x, inline function,
%       function handle or name of function in
%       form a(x,param)                   (optional, default eye())
%   Q - NxN process noise covariance matrix or NxNxK matrix
%       of K state process noise covariance matrices for each step.
%   param - Parameters of a. Parameters should be a single cell array,
%           vector or a matrix containing the same parameters for each
%           step, or if different parameters are used on each step they
%           must be a cell array of the format { param_1, param_2, ...},
%           where param_x contains the parameters for step x as a cell array,
%           a vector or a matrix.   (optional, default empty)
%   alpha - Transformation parameter      (optional)
%   beta  - Transformation parameter      (optional)
%   kappa - Transformation parameter      (optional)
%   mat   - If 1 uses matrix form         (optional, default 0)
%   same_p - If 1 uses the same parameters 
%            on every time step      (optional, default 1) 
%

% Out:
%   K - Smoothed state mean sequence
%   P - Smoothed state covariance sequence
%   D - Smoother gain sequence
%   
% Description:
%   Unscented Rauch-Tung-Striebel smoother algorithm. Calculate
%   "smoothed" sequence from given Kalman filter output sequence by
%   conditioning all steps to all measurements.
%
% Example:
%   m = m0;
%   P = P0;
%   MM = zeros(size(m,1),size(Y,2));
%   PP = zeros(size(m,1),size(m,1),size(Y,2));
%   for k=1:size(Y,2)
%     [m,P] = ukf_predict1(m,P,a,Q);
%     [m,P] = ukf_update1(m,P,Y(:,k),h,R);
%     MM(:,k) = m;
%     PP(:,:,k) = P;
%   end
%   [SM,SP] = urts_smooth(MM,PP,a,Q);
%
% See also:
%   UKF_PREDICT1, UKF_UPDATE1

% Copyright (C) 2006 Simo Särkkä
%
% $Id: urts_smooth1.m,v 1.3 2006/10/10 20:18:57 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [M,P,D] = urts_smooth1(M,P,a,Q,param,alpha,beta,kappa,mat,same_p)

  %
  % Check which arguments are there
  %
  if nargin < 4
    error('Too few arguments');
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
  if nargin < 10
    same_p = 1;
  end

  %
  % Apply defaults
  %
  if isempty(a)
    a = eye(size(M,1));
  end
  if isempty(Q)
    Q = zeros(size(M,1));
  end
  if isempty(mat)
    mat = 0;
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
    [m_pred,P_pred,C] = ...
	ut_transform(M(:,k),P(:,:,k),a,params,alpha,beta,kappa);
    P_pred = P_pred + Q(:,:,k);
    D(:,:,k) = C / P_pred;
    M(:,k)   = M(:,k) + D(:,:,k) * (M(:,k+1) - m_pred);
    P(:,:,k) = P(:,:,k) + D(:,:,k) * (P(:,:,k+1) - P_pred) * D(:,:,k)';
  end
