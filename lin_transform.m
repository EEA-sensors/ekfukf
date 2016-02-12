%UT_TRANSFORM  Perform linearization based transform of a Gaussian rv
%
%
% Syntax:
%   [mu,S,C,X,Y,w] = LIN_TRANSFORM(M,P,g,g_param,tr_param)
%
% In:
%   M - Random variable mean (Nx1 column vector)
%   P - Random variable covariance (NxN pos.def. matrix)
%   g - Transformation function of the form g(x,param) as
%       matrix, inline function, function name or function reference
%   g_param - Parameters of g               (optional, default empty)
%   tr_param - Parameters of linearization as a cell array with elements: 
%       tr_param{1} = derivative function gx of g
%
% Out:
%   mu - Estimated mean of y
%    S - Estimated covariance of y
%    C - Estimated cross-covariance of x and y
%
% Description:
%   ...
%   For default values of parameters, see UT_WEIGHTS.
%
% See also
%   UT_WEIGHTS UT_MWEIGHTS UT_SIGMAS

% Copyright (C) 2006 Simo S�rkk�
%
% $Id: ut_transform.m 111 2007-09-04 12:09:23Z ssarkka $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [mu,S,C] = lin_transform(M,P,g,g_param,tr_param)

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
      gx = [];
  else
      gx = tr_param{1};
  end
  
  if isempty(gx)
      gx = eye(size(M,1));
  end
  
  if isnumeric(gx)
    % nop
  elseif isstr(gx) | strcmp(class(gx),'function_handle')
    gx = feval(gx,M,g_param);
  else
    gx = gx(M,g_param);
  end
    
  if isempty(g)
    mu = g*M;
  elseif isnumeric(g)
    mu = g;
  elseif isstr(g) | strcmp(class(g),'function_handle')
    mu = feval(g,M,g_param);
  else
    mu = g(M,g_param);
  end
  
  S = gx*P*gx';
  C = P*gx';
  
  
 
      
