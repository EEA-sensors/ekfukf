%RK4  4th order Runge-Kutta integration
%
% Syntax:
%   [x,Y] = rk4(f,dt,x,[P1,P2,P3,Y])
%
% In:
%    f - Name of function in form f(x,P(:)) or
%        inline function taking the same parameters.
%        In chained case the function should be f(x,y,P(:)).
%   dt - Delta time as scalar.
%    x - Value of x from the previous time step.
%   P1 - Values of parameters of the function at initial time t
%        as a cell array (or single plain value). Defaults to empty
%        array (no parameters).
%   P2 - Values of parameters of the function at time t+dt/2 as
%        a cell array (or single plain value). Defaults to P1 and
%        each empty (or missing) value in the cell array is replaced
%        with the corresponding value in P1.
%   P3 - Values of parameters of the function at time t+dt.
%        Defaults to P2 similarly to above.
%    Y - Cell array of partial results y1,y2,y3,y4 in the RK algorithm
%        of the second parameter in the interated function. This can be
%        used for chaining the integrators. Defaults to empty.
% 
% Out:
%    x - Next value of X
%    Y - Cell array of partial results in Runge-Kutta algorithm.
%
% Description:
%   Perform one fourth order Runge-Kutta iteration step
%   for differential equation
%
%       dx/dt = f(x(t),P{:})
%
%   or in the chained case
%
%       dx/dt = f(x(t),y(t),P{:})
%       dy/dt = g(y(t),P{:})
%
%   - Example 1. Simple integration of model
%
%       dx/dt = tanh(x), x(0) = 1
%
%     can be done as follows:
%
%       X = [];
%       x = 1;
%       f = inline('tanh(x)','x');
%       for i=1:100
%         x = rk4(f,0.1,x);
%         X = [X x];
%       end
%
%   - Example 2. Chaining of integrators. Consider a
%     model of the form
%
%       dx/dt = x+y,     x(0)=1
%       dy/dt = tanh(y), y(0)=2
%
%     The equations can be now integrated as follows:
%
%       XY = [];
%       x = 1;
%       y = 2;
%       fx = inline('x+y','x','y');
%       fy = inline('tanh(y)','y');
%       for i=1:100
%         [y,YY] = rk4(fy,0.1,y);
%         x = rk4(fx,0.1,x,{},{},{},YY);
%         XY = [XY [x;y]];
%       end
%
%   which produces exactly the same result as
%
%      XY = [];
%      xy = [1;2];
%      fxy = inline('[xy(1)+xy(2);tanh(xy(2))]','xy');
%      for i=1:100
%         xy = rk4(fxy,0.1,xy);
%         XY = [XY xy];
%      end

% History:
%   14.10.2005  The first official version.
%
% Copyright (C) 2005 Simo Särkkä
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [x,Y] = rk4(f,dt,x,P1,P2,P3,Y)

  %
  % Apply default parameter values
  %
  if nargin < 4
    P1 = {};
  end
  if nargin < 5
    P2 = {};
  end
  if nargin < 6
    P3 = {};
  end
  if nargin < 7
    Y = {};
  end
  
  if ~isempty(P1) & ~iscell(P1)
    P1 = {P1};
  end
  if ~isempty(P2) & ~iscell(P2)
    P2 = {P2};
  end
  if ~isempty(P3) & ~iscell(P3)
    P3 = {P3};
  end
  
  if isempty(P2)
    P2 = P1;
  else
    for i=1:length(P1)
      if length(P2) >= i
	if isempty(P2{i})
	  P2{i} = P1{i};
	end
      else
	P2{i} = P1{i};
      end
    end
  end
  if isempty(P3)
    P3 = P2;
  else
    for i=1:length(P2)
      if length(P3) >= i
	if isempty(P3{i})
	  P3{i} = P2{i};
	end
      else
	P3{i} = P2{i};
      end
    end
  end

  %
  % Perform Runge-Kutta step
  %
  if ~isempty(Y)
    %
    % Chained integration
    %
    if ischar(f) | strcmp(class(f),'function_handle')
      x1  = x;
      dx1 = feval(f,x1,Y{1},P1{:}) * dt;
      x2  = x+0.5*dx1;
      dx2 = feval(f,x2,Y{2},P2{:}) * dt;
      x3  = x+0.5*dx2;
      dx3 = feval(f,x3,Y{3},P2{:}) * dt;
      x4  = x+dx3;
      dx4 = feval(f,x4,Y{4},P3{:}) * dt;
    else
      x1  = x;
      dx1 = f(x1,Y{1},P1{:}) * dt;
      x2  = x+0.5*dx1;
      dx2 = f(x2,Y{2},P2{:}) * dt;
      x3  = x+0.5*dx2;
      dx3 = f(x3,Y{3},P2{:}) * dt;
      x4  = x+dx3;
      dx4 = f(x4,Y{4},P3{:}) * dt;
    end  
  else
    if ischar(f) | strcmp(class(f),'function_handle')
      x1  = x;
      dx1 = feval(f,x1,P1{:}) * dt;
      x2  = x+0.5*dx1;
      dx2 = feval(f,x2,P2{:}) * dt;
      x3  = x+0.5*dx2;
      dx3 = feval(f,x3,P2{:}) * dt;
      x4  = x+dx3;
      dx4 = feval(f,x4,P3{:}) * dt;
    else
      x1  = x;
      dx1 = f(x1,P1{:}) * dt;
      x2  = x+0.5*dx1;
      dx2 = f(x2,P2{:}) * dt;
      x3  = x+0.5*dx2;
      dx3 = f(x3,P2{:}) * dt;
      x4  = x+dx3;
      dx4 = f(x4,P3{:}) * dt;
    end
  end

  Y = {x1,x2,x3,x4};
  x = x + 1/6 * (dx1 + 2*dx2 + 2*dx3 + dx4);
