
%UT_TRANSFORM  Perform unscented transform
%
% Syntax:
%   [mu,S,C,X,Y,w] = UT_TRANSFORM(M,P,g,[param,alpha,beta,kappa,mat],n,X,w)
%
% In:
%   M - Random variable mean (Nx1 column vector)
%   P - Random variable covariance (NxN pos.def. matrix)
%   g - Transformation function of the form g(x,param) as
%       matrix, inline function, function name or function reference
%   param - Parameters of g               (optional, default empty)
%   alpha - Transformation parameter      (optional)
%   beta  - Transformation parameter      (optional)
%   kappa - Transformation parameter      (optional)
%   mat   - If 1 uses matrix form         (optional, default 0)
%   X - Sigma points of x
%   w - Weights as cell array {mean-weights,cov-weights,c}
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
%   UT_WEIGHTS

% Copyright (C) 2006 Simo Särkkä
%
% $Id: ut_transform.m,v 1.2 2006/10/10 20:18:58 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [mu,S,C,X,Y,w] = ut_transform(M,P,g,param,alpha,beta,kappa,mat,X,w)

  if nargin < 4
    param = [];
  end
  if nargin < 5
    alpha = [];
  end
  if nargin < 6
    beta = [];
  end
  if nargin < 7
    kappa = [];
  end
  if nargin < 8
    mat = [];
  end
  if nargin < 10
    generate_sigmas = 1;
  else
    generate_sigmas = 0; 
  end
  

  %
  % Apply defaults
  %
  if isempty(mat)
    mat = 0;
  end
  
  %
  % Calculate sigma points
  %
  if generate_sigmas == 0
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
  elseif isstr(g) | strcmp(class(g),'function_handle')
    Y = [];
    for i=1:size(X,2)
      Y = [Y feval(g,X(:,i),param)];
    end
  else
    Y = [];
    for i=1:size(X,2)
      Y = [Y g(X(:,i),param)];
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

  
