%DER_CHECK  Check derivatives using finite differences
%
% Syntax:
%   [D0,D1] = DER_CHECK(F,DF,INDEX,[P1,P2,P3,...])
%
% In:
%   F  - Name of actual function or inline function
%        in form F(P1,P2,...)
%   DF - Derivative value as matrix, name of derivative
%        function or inline function in form DF(P1,P2,...).
%   INDEX - Index of parameter of interest. DF should
%        Calculate the derivative with recpect to parameter
%        Pn, where n is the index.
%
% Out:
%   D0 - Actual derivative
%   D1 - Estimated derivative
%
% Description:
%   Evaluates function derivative analytically and
%   using finite differences. If no output arguments
%   are given, issues a warning if these two values
%   differ too much.
%
%   Function is intended to checking that derivatives
%   of transition and measurement equations of EKF are
%   bug free.
% 
% See also:
%   EKF_PREDICT1, EKF_UPDATE1, EKF_PREDICT2, EKF_UPDATE2

% History:
%   12.03.2003  SS  Support for function handles
%   27.11.2002  SS  The first official version.
%
% Copyright (C) 2002 Simo Särkkä
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [D0,D1] = der_check(f,df,index,varargin)
  %
  % Calculate function value and derivative
  %
  if ischar(f) | strcmp(class(f),'function_handle')
    y0 = feval(f,varargin{:});
  else
    y0 = f(varargin{:});
  end
  if isnumeric(df)
    D0 = df;
  elseif ischar(df) | strcmp(class(df),'function_handle')
    D0 = feval(df,varargin{:});
  else
    D0 = df(varargin{:});
  end

  %
  % Calculate numerical derivative
  %
  h = 0.0000001;
  D1 = [];
  X = varargin{index};
  for r=1:size(X,1)
    for c=1:size(X,2)
      H = zeros(size(X,1),size(X,2));
      H(r,c) = h;
      varargin{index} = X+H;
      if ischar(f) | strcmp(class(f),'function_handle')
        y1 = feval(f,varargin{:});
      else
        y1 = f(varargin{:});
      end
      d1 = (y1-y0)/h;
      if size(d1,1)>1
        D1 = [D1 d1];
      else
        D1(r,c) = d1;
      end
    end
  end

  if nargout == 0
    d = abs(D1(:)-D0(:));
    if max(d) > 0.001
      warning('Derivatives differ too much');
      D0
      D1
    else
      fprintf('Derivative check passed.\n');
    end
  end

