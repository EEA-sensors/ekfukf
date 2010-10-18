%UT_TRANSFORM  Perform unscented transform
%
% Syntax:
%   [mu,S,C,X,Y,w] = UT_TRANSFORM(M,P,g,g_param,tr_param)
%
% In:
%   M - Random variable mean (Nx1 column vector)
%   P - Random variable covariance (NxN pos.def. matrix)
%   g - Transformation function of the form g(x,param) as
%       matrix, inline function, function name or function reference
%   g_param - Parameters of g               (optional, default empty)
%   tr_param - Parameters of the transformation as:       
%       alpha = tr_param{1} - Transformation parameter      (optional)
%       beta  = tr_param{2} - Transformation parameter      (optional)
%       kappa = tr_param{3} - Transformation parameter      (optional)
%       mat   = tr_param{4} - If 1 uses matrix form         (optional, default 0)
%       X     = tr_param{5} - Sigma points of x
%       w     = tr_param{6} - Weights as cell array {mean-weights,cov-weights,c}
%
% Out:
%   mu - Estimated mean of y
%    S - Estimated covariance of y
%    C - Estimated cross-covariance of x and y
%    X - Sigma points of x
%    Y - Sigma points of y
%    w - Weights as cell array {mean-weights,cov-weights,c}
%
% Description:
%   ...
%   For default values of parameters, see UT_WEIGHTS.
%
% See also
%   UT_WEIGHTS UT_MWEIGHTS UT_SIGMAS

% Copyright (C) 2006 Simo S�rkk�
%               2010 Jouni Hartikainen
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [mu,S,C,X,Y,w] = ut_transform(M,P,g,g_param,tr_param)

  if nargin < 4
     g_param = [];
  end
  
  if nargin < 5
     tr_param = []; 
  end

  %
  % Apply defaults
  %
  if isempty(tr_param) 
     alpha = [];
     beta = [];
     kappa = [];
     mat = [];
     X = [];
     w = [];
  else 
      alpha = tr_param{1}; 
      if length(tr_param) >= 2
          beta = tr_param{2};
      else 
          beta = [];
      end
      if length(tr_param) >= 3          
          kappa = tr_param{3};
      else
          kappa = [];
      end
      if length(tr_param) >= 4
          mat = tr_param{4};
      else
          mat = [];
      end
      if length(tr_param) >= 5
          X = tr_param{5};
      else
          X = [];
      end
      if length(tr_param) >= 6
          w = tr_param{6};
      else
          w = [];
      end
  end
    
  if isempty(mat)
     mat = 0;
  end  
  
  %
  % Calculate sigma points
  %
  if isempty(w) == 0
      WM = w{1};
      c  = w{3};
      if mat
        W  = w{2};
      else
        WC = w{2};
      end
  elseif mat
    [WM,W,c] = ut_mweights(size(M,1),alpha,beta,kappa);
    X = ut_sigmas(M,P,c);
    w = {WM,W,c};
  else
    [WM,WC,c] = ut_weights(size(M,1),alpha,beta,kappa);
    X = ut_sigmas(M,P,c);
    w = {WM,WC,c};
  end
  
  %
  % Propagate through the function
  %
  if isnumeric(g)
    Y = g*X;
  elseif ischar(g) | strcmp(class(g),'function_handle')
    Y = [];
    for i=1:size(X,2)
      Y = [Y feval(g,X(:,i),g_param)];
    end
  else
    Y = [];
    for i=1:size(X,2)
      Y = [Y g(X(:,i),g_param)];
    end
  end
  
  if mat
    mu = Y*WM;
    S  = Y*W*Y';
    C  = X*W*Y';
  else
    mu = zeros(size(Y,1),1);
    S  = zeros(size(Y,1),size(Y,1));
    C  = zeros(size(M,1),size(Y,1));
    for i=1:size(X,2)
      mu = mu + WM(i) * Y(:,i);
    end
    for i=1:size(X,2)
      S = S + WC(i) * (Y(:,i) - mu) * (Y(:,i) - mu)';
      C = C + WC(i) * (X(1:size(M,1),i) - M) * (Y(:,i) - mu)';
    end
  end

  
