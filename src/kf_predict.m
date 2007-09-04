%KF_PREDICT  Perform Kalman Filter prediction step
%
% Syntax:
%   [X,P] = KF_PREDICT(X,P,A,Q,B,U)
%
% In:
%   X - Nx1 mean state estimate of previous step
%   P - NxN state covariance of previous step
%   A - Transition matrix of discrete model (optional, default identity)
%   Q - Process noise of discrete model     (optional, default zero)
%   B - Input effect matrix                 (optional, default identity)
%   U - Constant input                      (optional, default empty)
%
% Out:
%   X - Predicted state mean
%   P - Predicted state covariance
%   
% Description:
%   Perform Kalman Filter prediction step. The model is
%
%     x[k] = A*x[k-1] + B*u[k-1] + q,  q ~ N(0,Q).
% 
%   The predicted state is distributed as follows:
%   
%     p(x[k] | x[k-1]) = N(x[k] | A*x[k-1] + B*u[k-1], Q[k-1])
%
%   The predicted mean x-[k] and covariance P-[k] are calculated
%   with the following equations:
%
%     m-[k] = A*x[k-1] + B*u[k-1]
%     P-[k] = A*P[k-1]*A' + Q.
%
%   If there is no input u present then the first equation reduces to
%     m-[k] = A*x[k-1]
%
% History:
%
%   26.2.2007 JH Added the distribution model for the predicted state
%                and equations for calculating the predicted state mean and
%                covariance to the description section.
%  
% See also:
%   KF_UPDATE, LTI_DISC, EKF_PREDICT, EKF_UPDATE

% Copyright (C) 2002-2006 Simo Särkkä
% Copyright (C) 2007 Jouni Hartikainen
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [x,P] = kf_predict(x,P,A,Q,B,u)

  %
  % Check arguments
  %
  if nargin < 3
    A = [];
  end
  if nargin < 4
    Q = [];
  end
  if nargin < 5
    B = [];
  end
  if nargin < 6
    u = [];
  end
  
  %
  % Apply defaults
  %
  if isempty(A)
    A = eye(size(x,1));
  end
  if isempty(Q)
    Q = zeros(size(x,1));
  end
  if isempty(B) & ~isempty(u)
    B = eye(size(x,1),size(u,1));
  end

  %
  % Perform prediction
  %
  if isempty(u)
    x = A * x;
    P = A * P * A' + Q;
  else
    x = A * x + B * u;
    P = A * P * A' + Q;
  end
