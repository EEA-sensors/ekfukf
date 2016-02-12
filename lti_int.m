%LTI_INT  Integrate LTI ODE with Gaussian Noise
%
% Syntax:
%   [x,P,A] = lti_int(x,P,F,L,Q,T)
%
% Description:
%   Integrates LTI differential equation
%
%     x' = F*x + L*w  , w ~ N(0,Q)
%
%   from t0=0 to t1=T or over given steps t0,t1,t2,t3,...
%   Initial conditions can be in form
%
%     x(t0) = x0
%
%   or
%
%     x(t0) ~ N(x0,P0)
%
% See also
%   LTI_DISC

% History:
%   20.11.2002  The first official version.
%
% Copyright (C) 2002 Simo Särkkä
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [x,P,A] = lti_int(x,P,F,L,Q,T)

  %
  % Check number of arguments
  %
  if nargin < 3
    error('Too few arguments');
  end
  if nargin < 4
    L = [];
  end
  if nargin < 5
    Q = [];
  end
  if nargin < 6
    T = [];
  end

  %
  % Integrate one or many steps
  %
  if length(T) > 1
    %
    % Handle multiple time steps recursively
    %
    if isempty(P)
      if isempty(x)
        %
        % No x or P given, return only A's
        %
        tA = eye(size(x,1));
        A = zeros(size(x,1),size(x,1),length(T));
        A(:,:,1) = tA;
        for i=2:length(T)
          [tx,tP,tA] = lti_int([],[],F,L,Q,T(i)-T(i-1));
          A(:,:,i) = tA;
        end
      else
        %
        % Only x given, returns x's and A's
        %
        tx = x;
        tA = eye(size(x,1));
        x = zeros(size(x,1),length(T));
        A = zeros(size(x,1),size(x,1),length(T));
        x(:,1)   = tx;
        A(:,:,1) = tA;
        for i=2:length(T)
          [tx,tP,tA] = lti_int(tx,[],F,L,Q,T(i)-T(i-1));
          x(:,i)   = tx;
          A(:,:,i) = tA;
        end
      end
    else
      %
      % Both x and P given, return all
      %
      tx = x;
      tP = P;
      tA = eye(size(x,1));
      x = zeros(size(x,1),length(T));
      P = zeros(size(x,1),size(x,1),length(T));
      A = zeros(size(x,1),size(x,1),length(T));
      x(:,1)   = tx;
      P(:,:,1) = tP;
      A(:,:,1) = tA;
      for i=2:length(T)
        [tx,tP,tA] = lti_int(tx,tP,F,L,Q,T(i)-T(i-1));
        x(:,i)   = tx;
        P(:,:,i) = tP;
        A(:,:,i) = tA;
      end
    end

  else
    %
    % One step integration from 0 to T
    %
    if isempty(L)
      L = eye(size(x,1));
    end
    if isempty(Q)
      Q = zeros(size(x,1),size(x,1));
    end
    if isempty(T)
      T = 1;
    end

    %
    % Closed form integration of mean
    %
    A = expm(F*T);
    x = A*x;

    %
    % Runge-Kutta Integration of Covariance
    %
    if ~isempty(P) & (nargout > 1)
      f = inline('F*P+P*F''+L*Q*L','P','t','F','L','Q');
      P = rk(f,P,0,T,F,L,Q);
    end
  end

