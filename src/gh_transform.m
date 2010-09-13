function [mu,S,C,SX,W] = gh_transform(m,P,g,p,param)
% GH_TRANSFORM - Gauss-Hermite transform of random variables
%
% Syntax:
%   [mu,S,C,SX,W] = GH_TRANSFORM(M,P,g,p,param)
%
% In:
%   M - Random variable mean (Nx1 column vector)
%   P - Random variable covariance (NxN pos.def. matrix)
%   g - Transformation function of the form g(x,param) as
%       matrix, inline function, function name or function reference
%   p - Number of points in Gauss-Hermite integration
%   param - Parameters of g               (optional, default empty)
%
% Out:
%   mu - Estimated mean of y
%    S - Estimated covariance of y
%    C - Estimated cross-covariance of x and y
%   SX - Sigma points of x
%    W - Weights as cell array

% Copyright (c), 2009 Hartikainen, Särkkä, Solin
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.
%%

  % Estimate the mean of g
  if nargin < 5
    [mu,SX,x,W,F] = ngausshermi(g,p,m,P);
  else
    [mu,SX,x,W,F] = ngausshermi(g,p,m,P,param);
  end
  
  % Estimate the P and C
  if nargin < 5
    [pc,SX,x,W,FPC] = ngausshermi(@gh_packed_pc,p,m,P,{g,m,mu});
  else
    [pc,SX,x,W,FPC] = ngausshermi(@gh_packed_pc,p,m,P,{g,m,mu,param});
  end
  d = size(m,1);
  s = size(mu,1);
  S = reshape(pc(1:s^2),s,s);
  C = reshape(pc(s^2+(1:s*d)),d,s);

  