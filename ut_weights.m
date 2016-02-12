%UT_WEIGHTS - Generate unscented transformation weights
%
% Syntax:
%   [WM,WC,c] = ut_weights(n,alpha,beta,kappa)
%
% In:
%   n     - Dimensionality of random variable
%   alpha - Transformation parameter  (optional, default 0.5)
%   beta  - Transformation parameter  (optional, default 2)
%   kappa - Transformation parameter  (optional, default 3-n)
%
% Out:
%   WM - Weights for mean calculation
%   WC - Weights for covariance calculation
%    c - Scaling constant
%
% Description:
%   Computes unscented transformation weights.
%
% See also UT_MWEIGHTS UT_TRANSFORM UT_SIGMAS
% 

% Copyright (C) 2006 Simo S�rkk�
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [WM,WC,c] = ut_weights(n,alpha,beta,kappa)

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

%
% Apply default values
%
if isempty(alpha)
  alpha = 1;
end
if isempty(beta)
  beta = 0;
end
if isempty(kappa)
  kappa = 3 - n;
end
	  
%
% Compute the normal weights 
%
lambda = alpha^2 * (n + kappa) - n;
	  
WM = zeros(2*n+1,1);
WC = zeros(2*n+1,1);
for j=1:2*n+1
  if j==1
    wm = lambda / (n + lambda);
    wc = lambda / (n + lambda) + (1 - alpha^2 + beta);
  else
    wm = 1 / (2 * (n + lambda));
    wc = wm;
  end
  WM(j) = wm;
  WC(j) = wc;
end

c = n + lambda;
