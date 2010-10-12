function [I,SX,x,W,F] = ngausshermi(f,n,m,P,param)
% NGAUSSHERMI - ND scaled Gauss-Hermite quadrature (cubature) rule
%
% Syntax:
%   [I,x,W,F] = ngausshermi(f,p,m,P,param)
%
% In:
%   f - Function f(x,param) as inline, name or reference
%   n - Polynomial order
%   m - Mean of the d-dimensional Gaussian distribution
%   P - Covariance of the Gaussian distribution
%   param - Optional parameters for the function
%
% Out:
%   I - The integral value
%   x - Evaluation points
%   W - Weights
%   F - Function values
%
% Description:
%   Approximates a Gaussian integral using the Gauss-Hermite method
%   in multiple dimensions:
%     int f(x) N(x | m,P) dx

% History:
%   2009 - Initial version (ssarkka)
%   2010 - Partially rewritten (asolin)

% Copyright (c) 2010 Simo Sarkka, Arno Solin
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

%% THEORY OF OPERATION:
%
% Consider the multidimensional integral:
%
%  int f(x1,...,xn) exp(-x1^2-...-xn^2) dx1...dxn
%
% If we form quadrature for dimension 1, then
% for dimension 2 and so on, we get
%
% = sum w_i1 int f(x1^i1,...,xn) exp(-x2^2-...-xn^2) dx2...dxn
% = sum w_i1 w_i1 int f(x1^i1,x2^i2,...,xn) exp(-x3^2-...-xn^2) dx3...dxn
% = ...
% = sum w_i1 ... w_in f(x1^i1,...,xn^in) 
%
% Now let X = [x1,...,xn] and consider
%
%  1/(2pi)^{n/2}/sqrt(|P|) int f(X) exp(-1/2 (X-M)' iP (X-M)) dX
%
% Let P = L L' and change variable to
%
%   X = sqrt(2) L Y + M  <=>  Y = 1/sqrt(2) iL (X-M)
%
% then dX = sqrt(2)^n |L] dY = sqrt(2)^n sqrt(|P]) dY i.e. we get
%
%  1/(2pi)^{n/2}/sqrt(|P|) int f(X) exp(-1/2 (X-M)' iP (X-M)) dX
%   = 1/(pi)^{n/2} int f(sqrt(2) L Y + M) exp(-Y'Y) dY
%   = int g(Y) exp(-Y'Y) dY
%
% which is of the previous form if we define
%
%   g(Y) = 1/(pi)^{n/2} f(sqrt(2) L Y + M)
% 

%% The Gauss-Hermite cubature rule

  % The hermite polynomial of order n
  p = hermitepolynomial(n);
  
  % Evaluation points
  x = roots(p);
  
  % 1D coefficients
  Wc = pow2(n-1)*factorial(n)*sqrt(pi)/n^2;
  p2 = hermitepolynomial(n-1);
  W  = zeros(1,n);
  for i=1:n
     W(i)=Wc*polyval(p2,x(i)).^-2;
  end
    
  d = size(m,1);
  if d == 1
    x = x';
  end
  
  % Generate all n^d collections of indexes by
  % transforming numbers 0...n^d-1) into n-base system
  % and by adding 1 to each digit
  num = 0:(n^d-1);
  ind = zeros(d,n^d);
  for i=1:d
      ind(i,:) = rem(num,n)+1;
      num = floor(num / n);
  end
  
  % Form the sigma points and weights
  L = chol(P)';
  SX = sqrt(2)*L*x(ind)+repmat(m,1,size(ind,2));
  W = prod(W(ind),1); % ND weights
  
  % Evaluate the function at the sigma points
  if ischar(f) || strcmp(class(f),'function_handle')
      if nargin < 5
         F = feval(f,SX);
      else
         F = feval(f,SX,param);
      end
  elseif isnumeric(f)
         F = f*SX;
  else
      if nargin < 5
         F = f(SX);
      else
         F = f(SX,param);
      end
  end

  % Evaluate the integral
  I = 1/sqrt(pi)^d*sum(F.*repmat(W,size(F,1),1),2);
  