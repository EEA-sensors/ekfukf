function [M,P,K,MU,S,LH] = adf_update(M,P,Y,h,R,h_param,tr_method,tr_param)
% ADF_UPDATE - (Gaussian) Assumed Density Kalman filter update step
%
% Syntax:
%   [M,P,K,MU,S,LH] = ADF_UPDATE(M,P,Y,h,R,h_param,tr_method,tr_param)
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
%   tr_method - Approximative method for estimating the joint Gaussian distribution
%               of x and y, where x ~ N(M,P) and y = g(x).
%   h_param - Parameters of h

%   tr_param - Parameters of the transformation
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
%   Perform a general (Gaussian) Assumed Density Filter (ADF)
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
%   [M2,P2] = adff_update(M1,P1,Y,h,R,S);
%
% See also:
%   ADF_PREDICT, ADRTS_SMOOTH

% History:
%   Oct 2, 2010 - Initial version  

% Copyright (C) 2010 Hartikainen, Särkkä, Solin
%
% $Id: adf_update.m,v 1.2 2009/07/01 06:34:41 ssarkka Exp $
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
     tr_param = []; 
  end

  %
  % Do the transform and make the update
  %
  [MU,S,C] = feval(tr_method,M,P,h,h_param,tr_param);

  S = S + R;
  K = C / S;
  M = M + K * (Y - MU);
  P = P - K * S * K';
  
  if nargout > 5
    LH = gauss_pdf(Y,MU,S);
  end
