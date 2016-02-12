%KF_LOOP  Performs the prediction and update steps of the Kalman filter
%         for a set of measurements.   
%
% Syntax:
%   [MM,PP] = KF_LOOP(X,P,H,R,Y,A,Q)
% 
% In:
%   X - Nx1 initial estimate for the state mean 
%   P - NxN initial estimate for the state covariance
%   H - DxN measurement matrix
%   R - DxD measurement noise covariance
%   Y - DxM matrix containing all the measurements.
%   A - Transition matrix of the discrete model (optional, default identity)
%   Q - Process noise of the discrete model     (optional, default zero)
%   
% Out:
%   MM - Filtered state mean sequence
%   PP - Filtered state covariance sequence
%  
%  Description:
%    Calculates state estimates for a set measurements using the
%    Kalman filter. This function is for convience, as it basically consists
%    only of a space reservation for the estimates and of a for-loop which
%    calls the predict and update steps of the KF for each time step in
%    the measurements.  
%  
%  See also:
%    KF_PREDICT, KF_UPDATE

%  History:
%   
%    12.2.2007 JH Initial version.  
%
% Copyright (C) 2007 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [MM,PP] = kf_loop(X,P,H,R,Y,A,Q)

  % Check the input parameters.  
  if nargin < 5
     error('Too few arguments');
  end
  if nargin < 6
     A = [];
  end
  if nargin < 7
     Q = [];
  end

  % Apply the defaults
  if isempty(A)
    A = eye(size(X,1));
  end
  if isempty(Q)
    Q = zeros(size(X,1));
  end

  % Space for the estimates.
  MM = zeros(size(X,1), size(Y,2));
  PP = zeros(size(X,1), size(X,1), size(Y,2));

  % Filtering steps.
  for i = 1:size(Y,2)
     [X,P] = kf_predict(X,P,A,Q);
     [X,P] = kf_update(X,P,Y(:,i),H,R);
     MM(:,i) = X;
     PP(:,:,i) = P;
  end
