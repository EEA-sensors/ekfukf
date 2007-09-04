%GAUSS_RND  Multivariate Gaussian random variables
%
% Syntax:
%   X = GAUSS_RND(M,S,N)
%
% In:
%   M - Dx1 mean of distibution or K values as DxK matrix.
%   S - DxD covariance matrix
%   N - Number of samples (optional, default 1)
%
% Out:
%   X - Dx(K*N) matrix of samples.
%   
% Description:
%   Draw N samples from multivariate Gaussian distribution
% 
%     X ~ N(M,S)
%
% See also:
%   GAUSS_PDF

% Copyright (C) 2002-2006 Simo Särkkä
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.


function X = gauss_rnd(M,S,N)

  if nargin < 3
    N = 1;
  end
  
  L = chol(S)';
  X = repmat(M,1,N) + L*randn(size(M,1),size(M,2)*N);
  
