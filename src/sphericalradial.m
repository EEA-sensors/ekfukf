function [I,x,W,F] = sphericalradial(f,m,P,param)
% SPHERICALRADIAL - ND scaled spherical-radial cubature rule
%
% Syntax:
%   [I,x,W,F] = sphericalradial(f,m,P[,param])
%
% In:
%   f - Function f(x,param) as inline, name or reference
%   m - Mean of the d-dimensional Gaussian distribution
%   P - Covariance of the Gaussian distribution
%   param - Parameters for the function (optional)
%
% Out:
%   I - The integral
%   x - Evaluation points
%   W - Weights
%   F - Function values
%
% Description:
%   Apply the spherical-radial cubature rule to integrals of form:
%     int f(x) N(x | m,P) dx

% History:
%   Aug 5, 2010 - Renamed from 'cubature' to 'sphericalradial' (asolin)

% Copyright (c) 2010 Arno Solin
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

%% Spherical-radial cubature rule

  % The dimension of m is
  n = size(m,1);
  
  % Evaluation points (nx2n)
  x = [eye(n) -eye(n)];
  
  % Scaling
  x = sqrt(n)*x;
  x = chol(P)'*x + repmat(m,1,2*n);

  % Evaluate the function at the points
    if ischar(f) || strcmp(class(f),'function_handle')
      if nargin < 4
        F = feval(f,x);
      else
        F = feval(f,x,param);
      end
    elseif isnumeric(f)
         F = f*x;
    else
      if nargin < 4
         F = f(x);
      else
         F = f(x,param);
      end
    end
  
  % The weights are
  W = 1/(2*n);
  
  % Return integral value
  I = W*sum(F,2);
  
  % Return weights
  W = W*ones(1,2*n);
  
  