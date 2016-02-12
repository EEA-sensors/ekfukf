%UT_TRANSFORM  Perform quadratic approximation based transform of a Gaussian rv
%
%
% Syntax:
%   [mu,S,C,X,Y,w] = QUAD_TRANSFORM(M,P,g,g_param,tr_param)
%
% In:
%   M - Random variable mean (Nx1 column vector)
%   P - Random variable covariance (NxN pos.def. matrix)
%   g - Transformation function of the form g(x,param) as
%       matrix, inline function, function name or function reference
%   g_param - Parameters of g               (optional, default empty)
%   tr_param - Parameters of linearization as a cell array with elements: 
%       tr_param{1} = derivative function gx of g
%       tr_param{2} = 2nd derivative function gxx of g
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

function [mu,S,C] = quad_transform(M,P,g,g_param,tr_param)

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
      gxx = [];
  else
      gx = tr_param{1};

      if length(tr_param) >= 2
          gxx = tr_param{2};  
      else 
          gxx = [];
      end
  end
  if isempty(gx)
      gx = eye(size(M,1));
  end
  
  if isempty(gxx)
      gxx = eye(size(M,1));
  end
  
  if isnumeric(gx)
    % nop
  elseif isstr(gx) | strcmp(class(gx),'function_handle')
    gx = feval(gx,M,g_param);
  else
    gx = gx(M,g_param);
  end
  
  if isnumeric(gxx)
    % nop
  elseif isstr(gxx) | strcmp(class(gxx),'function_handle')
    gxx = feval(gxx,M,g_param);
  else
    gxx = gxx(M,g_param);
  end
    
  if isempty(g)
    mu = g*M;
  elseif isnumeric(g)
    mu = g;
  elseif ischar(g) | strcmp(class(g),'function_handle')
    mu = feval(g,M,g_param);
  else
    mu = g(M,g_param);
  end
  
  for i = 1:size(gxx,1)
    mu(i) = mu(i) + 0.5*trace(squeeze(gxx(i,:,:))*P);
  end
  
  S = gx*P*gx';
  
  for i = 1:size(gxx,1)
      gxx_i = squeeze(gxx(i,:,:));
      for j = 1:size(gxx,1)
        gxx_j = squeeze(gxx(j,:,:));
        S(i,j) = S(i,j) + 0.5*trace(gxx_i*P*gxx_j*P);
      end
  end
  C = P*gx';
  
  
 
      
