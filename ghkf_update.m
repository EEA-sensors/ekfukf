function [M,P,K,MU,S,LH] = ghkf_update(M,P,Y,h,R,h_param,p)
% GHKF_UPDATE - Gauss-Hermite Kalman filter update step
%
% Syntax:
%   [M,P,K,MU,S,LH] = GHKF_UPDATE(M,P,Y,h,R,param,p)
%
% In:
%   M  - Mean state estimate after prediction step
%   P  - State covariance after prediction step
%   Y  - Measurement vector.
%   h  - Measurement model function as a matrix H defining
%        linear function h(x) = H*x, inline function,
%        function handle or name of function in
%        form h(x,param)
%   R  - Measurement covariance
%   h_param - Parameters of h
%   p  - Degree of approximation (number of quadrature points)
%
% Out:
%   M  - Updated state mean
%   P  - Updated state covariance
%   K  - Computed Kalman gain
%   MU - Predictive mean of Y
%   S  - Predictive covariance Y
%   LH - Predictive probability (likelihood) of measurement.
%   
% Description:
%   Perform additive form Gauss-Hermite Kalman filter (GHKF)
%   measurement update step. Assumes additive measurement
%   noise.
%
%   Function h(.) should be such that it can be given a
%   DxN matrix of N sigma Dx1 points and it returns 
%   the corresponding measurements for each sigma
%   point. This function should also make sure that
%   the returned sigma points are compatible such that
%   there are no 2pi jumps in angles etc.
%
% Example:
%   h = inline('atan2(x(2,:)-s(2),x(1,:)-s(1))','x','s');
%   [M2,P2] = ghkf_update(M1,P1,Y,h,R,S);
%
% See also:
%   GHKF_PREDICT, GHRTS_SMOOTH, GH_TRANSFORM

% History:
%   Jun 18, 2009 - Initial version  
%   May 24, 2010 - Fixed parameter input and description (asolin)
%   Aug 5,  2010 - Renamed from 'gh_update' to 'ghkf_update' (asolin)

% Copyright (C) 2009 Hartikainen, Särkkä, Solin
%
% $Id: gh_update.m,v 1.2 2009/07/01 06:34:41 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.
%%

  %
  % Check that all arguments are there
  %
  if nargin < 5
     error('Too few arguments');
  end
  if nargin < 6
     h_param = [];
  end
  if nargin < 7
     p = []; 
  end
  if isempty(p)
     p = 3;
  end

  %
  % Do the transform and make the update
  %
  tr_param = {p};
  [MU,S,C,X] = gh_transform(M,P,h,h_param,tr_param);
  S = S + R;
  K = C / S;
  M = M + K * (Y - MU);
  P = P - K * S * K';
  
  if nargout > 5
    LH = gauss_pdf(Y,MU,S);
  end
