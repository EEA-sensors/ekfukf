%FP_SMOOTH  Fixed-point Rauch-Tung-Striebel smoother
%
% Syntax:
%   [MS,PS,B] = FP_SMOOTH(MS,PS,B,MP,PP,M,P,A,Q)
%
% In:
%   MS - Current fixed-point smoothing mean, initialize
%        to filtered mean in the beginning
%   PS - Current fixed-point smoothing covariance, initialize
%        to filtered covariance in the beginning
%    B - Current smoothing gain, set to empty matrix or
%        to identity matrix to initialize
%   OM - Previous filtered state mean
%   OP - Previous filtered covariance
%    M - Current filtered state mean
%    P - Updated filtered covariance
%    A - NxN state transition matrix of transition to current step
%    Q - NxN process noise covariance matrix of transition to
%        current step
%
% Out:
%   MS - Smoothed mean
%   PS - Smoothed covariance
%    B - New smoothing gain
%   
% Description:
%   Fixed-point Rauch-Tung-Striebel smoother algorithm.
%   Calculate smoothed estimate of single time point.
%   by conditioning it to all measurements.
%
% Example:
%   m = m0;
%   P = P0;
%   ms0 = m0;
%   Ps0 = P0;
%   B = [];
%   for k=1:size(Y,2)
%     old_m = m;
%     old_P = P;
%     [m,P] = kf_predict(m,P,A,Q);
%     [m,P] = kf_update(m,P,Y(:,k),H,R);
%     [ms0,Ps0,B] = fp_smooth(ms0,Ps0,B,old_m,old_P,m,P,A,Q)
%   end
%
% See also:
%   RTS_SMOOTH

% Copyright (C) 2006 Simo Särkkä
%
% $Id: fp_smooth.m,v 1.1 2006/11/12 12:10:31 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [ms,Ps,B] = fp_smooth(ms,Ps,B,om,oP,m,P,A,Q)

  if isempty(B)
    B = eye(size(ms,1));
  end

  mp = A*om;
  Pp = A*oP*A'+Q;

  B = B*oP*A'/Pp;

  ms = ms + B*(m-mp);
  Ps = Ps + B*(P-Pp)*B';
