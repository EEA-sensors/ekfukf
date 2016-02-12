%URTS_SMOOTH2  Augmented form Unscented Rauch-Tung-Striebel smoother
%
% Syntax:
%   [M,P,S] = URTS_SMOOTH2(M,P,f,Q,[f_param,alpha,beta,kappa,mat,same_p])
%
% In:
%   M - NxK matrix of K mean estimates from Unscented Kalman filter
%   P - NxNxK matrix of K state covariances from Unscented Kalman Filter
%   f - Dynamic model function as inline function,
%       function handle or name of function in
%       form a([x;w],param)
%   Q - Non-singular covariance of process noise w
%   f_param - Parameters of a. Parameters should be a single cell array,
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
%   [...]
%
% See also:
%   URTS_SMOOTH1, UKF_PREDICT2, UKF_UPDATE2

% Copyright (C) 2006 Simo S�rkk�
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [M,P,D] = urts_smooth2(M,P,f,Q,f_param,alpha,beta,kappa,mat,same_p)

  %
  % Check which arguments are there
  %
  if nargin < 4
    error('Too few arguments');
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
  if nargin < 10
    same_p = 1;
  end

  %
  % Apply defaults
  %
  if isempty(mat)
    mat = 0;
  end
  
  %
  % Run the smoother
  %
  D = zeros(size(M,1),size(M,1),size(M,2));
  for k=(size(M,2)-1):-1:1
    if isempty(f_param)
        params = [];
    elseif same_p
        params = f_param;
    else
        params = f_param{k};
    end
      
    MA = [M(:,k);zeros(size(Q,1),1)];
    PA = zeros(size(P,1)+size(Q,1));
    PA(1:size(P,1),1:size(P,1)) = P(:,:,k);
    PA(1+size(P,1):end,1+size(P,1):end) = Q;
    
    tr_param = {alpha beta kappa mat};
    [m_pred,P_pred,C] = ...
	ut_transform(MA,PA,f,params,tr_param);
    C = C(1:size(M,1),:);
    
    D(:,:,k) = C / P_pred;
    M(:,k)   = M(:,k) + D(:,:,k) * (M(:,k+1) - m_pred);
    P(:,:,k) = P(:,:,k) + D(:,:,k) * (P(:,:,k+1) - P_pred) * D(:,:,k)';
  end
