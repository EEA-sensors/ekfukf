%TF_SMOOTH  Two filter based Smoother
%
% Syntax:
%   [M,P] = TF_SMOOTH(M,P,Y,A,Q,H,R,[use_inf])
%
% In:
%   M - NxK matrix of K mean estimates from Kalman filter
%   P - NxNxK matrix of K state covariances from Kalman Filter
%   Y - Sequence of K measurement as DxK matrix
%   A - NxN state transition matrix.
%   Q - NxN process noise covariance matrix.
%   H - DxN Measurement matrix.
%   R - DxD Measurement noise covariance.
%   use_inf - If information filter should be used (default 1)
%
% Out:
%   M - Smoothed state mean sequence
%   P - Smoothed state covariance sequence
%   
% Description:
%   Two filter linear smoother algorithm. Calculate "smoothed"
%   sequence from given Kalman filter output sequence
%   by conditioning all steps to all measurements.
%
% Example:
%   m = m0;
%   P = P0;
%   MM = zeros(size(m,1),size(Y,2));
%   PP = zeros(size(m,1),size(m,1),size(Y,2));
%   for k=1:size(Y,2)
%     [m,P] = kf_predict(m,P,A,Q);
%     [m,P] = kf_update(m,P,Y(:,k),H,R);
%     MM(:,k) = m;
%     PP(:,:,k) = P;
%   end
%   [SM,SP] = tf_smooth(MM,PP,A,Q,H,R,Y);
%
% See also:
%   KF_PREDICT, KF_UPDATE

% History:
%   
%   02.8.2007 JH Changed the name to tf_smooth
%   26.3.2007 JH Fixed a bug in backward filter with observations having
%                having more than one dimension.
%             
% Copyright (C) 2006 Simo Särkkä
%               2007 Jouni Hartikainen
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.
%

function [M,P] = tf_smooth(M,P,Y,A,Q,H,R,use_inf)

  %
  % Check which arguments are there
  %
  if nargin < 4
    error('Too few arguments');
  end
  if nargin < 8
    use_inf = [];
  end
  
  if isempty(use_inf)
    use_inf = 1;
  end
  
  %
  % Run the backward filter
  %
  if use_inf
    zz = zeros(size(M));
    SS = zeros(size(P));
    IR = inv(R);
    IQ = inv(Q);
    z = zeros(size(M,1),1);
    S = zeros(size(M,1),size(M,1));
    for k=size(M,2):-1:1
      G = S / (S + IQ);
      S = A' * (eye(size(M,1)) - G) * S * A;
      z = A' * (eye(size(M,1)) - G) * z;
      zz(:,k)   = z;
      SS(:,:,k) = S;
      S = S + H'*IR*H;
      z = z + H'*IR*Y(:,k);
    end
  else
    BM = zeros(size(M));
    BP = zeros(size(P));
    IA = inv(A);
    IQ = IA*Q*IA';  
    fm = zeros(size(M,1),1);
    fP = 1e12*eye(size(M,1));
    BM(:,end) = fm;
    BP(:,:,end) = fP;
    for k=(size(M,2)-1):-1:1
      [fm,fP] = kf_update(fm,fP,Y(:,k+1),H,R);
      [fm,fP] = kf_predict(fm,fP,IA,IQ);
      BM(:,k) = fm;
      BP(:,:,k) = fP;
    end
  end

  %
  % Combine estimates
  %
  if use_inf
    for k=1:size(M,2)-1
      G = P(:,:,k) * SS(:,:,k) / (eye(size(M,1)) + P(:,:,k) * SS(:,:,k));
      P(:,:,k) = inv(inv(P(:,:,k)) + SS(:,:,k));
      M(:,k) = M(:,k) + P(:,:,k) * zz(:,k) - G * M(:,k);
    end
  else
    for k=1:size(M,2)-1
      tmp = inv(inv(P(:,:,k)) + inv(BP(:,:,k)));
      M(:,k) = tmp * (P(:,:,k)\M(:,k) + BP(:,:,k)\BM(:,k));
      P(:,:,k) = tmp;
    end
  end

