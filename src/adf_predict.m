function [M,P,C] = adf_predict(M,P,f,Q,f_param,tr_method,tr_param)
% ADF_PREDICT - (Gaussian) Assumed Density Filter prediction step
%
% Syntax:
%   [M,P] = ADF_PREDICT(M,P,f,Q,f_param,tr_method,tr_param)
%
% In:
%   M - Nx1 mean state estimate of previous step
%   P - NxN state covariance of previous step
%   f - Dynamic model function as a matrix A defining
%       linear function f(x) = A*x, inline function,
%       function handle or name of function in
%       form f(x,param)                     (optional, default eye())
%   Q - Process noise of discrete model     (optional, default zero)
%   f_param - Parameters of f               (optional, default empty)
%   tr_method - Approximative method for estimating the joint Gaussian distribution
%               of x and y, where x ~ N(M,P) and y = g(x).
%   tr_param - Parameters of transformation 
% 
% Out:
%   M - Updated state mean
%   P - Updated state covariance
%   C - Covariance betwee x_{k-1} and x_k
%
% Description:
%   Perform a general additive form (Gaussian) Assumed Density Filter
%   prediction step.
%
%   Function f(.) should be such that it can be given a
%   DxN matrix of N sigma Dx1 points and it returns 
%   the corresponding predictions for each sigma
%   point. 
%
% See also:
%   ADF_UPDATE, ADRTS_SMOOTH, X_TRANSFORM

% History:
%   4 Aug, 2010 - First version

% Copyright (C) 2009 Hartikainen, Särkkä, Solin
%
% $Id: gh_predict.m,v 1.2 2009/07/01 06:34:40 ssarkka Exp $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.
%%

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
  if nargin < 7
     tr_param = []; 
  end

  %
  % Apply defaults
  %
  if isempty(f)
    f = eye(size(M,1));
  end
  if isempty(Q)
    Q = zeros(size(M,1));
  end

  %
  % Do transform and add process noise
  %
  [M,P,C] = feval(tr_method,M,P,f,f_param,tr_param);
  P = P + Q;


