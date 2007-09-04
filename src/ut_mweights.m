%UT_MWEIGHTS - Generate matrix form unscented transformation weights
%
% Syntax:
%   [WM,W,c] = ut_mweights(n,alpha,beta,kappa)
%
% In:
%   n     - Dimensionality of random variable
%   alpha - Transformation parameter  (optional, default 0.5)
%   beta  - Transformation parameter  (optional, default 2)
%   kappa - Transformation parameter  (optional, default 3-size(X,1))
%
% Out:
%   WM - Weight vector for mean calculation
%    W - Weight matrix for covariance calculation
%    c - Scaling constant
%
% Description:
%   Computes matrix form unscented transformation weights.
%

% Copyright (C) 2006 Simo Särkkä
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [WM,W,c] = ut_mweights(n,alpha,beta,kappa)

  %
  % Check which arguments are there
  %
  if nargin < 1
    error('At least dimensionality n required.');
  end
  if nargin < 2
    alpha = [];
  end
  if nargin < 3
    beta = [];
  end
  if nargin < 4
    kappa = [];
  end
  
  [WM,WC,c] = ut_weights(n,alpha,beta,kappa);

  W = eye(length(WC)) - repmat(WM,1,length(WM));
  W = W * diag(WC) * W';
