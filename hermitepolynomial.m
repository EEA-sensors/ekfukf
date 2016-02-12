function p=hermitepolynomial(n)
% HERMITEPOLYNOMIAL - Hermite polynomial
%
% Syntax:
%   p = hermitepolynomial(n)
%
% In:
%   n - Polynomial order
%
% Out:
%   p - Polynomial coefficients (starting from greatest order)
%
% Description:
%   Forms the Hermite polynomial of order n.
%
% See also:
%   POLYVAL, ROOTS

% History:
%   May 18, 2010 - Initial version (asolin)

% Copyright (c) 2010 Arno Solin
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

%% The "physicists' Hermite polynomials"
% To get the differently scaled "probabilists' Hermite polynomials"
% remove the coefficient *2 in (**).

  % Check the input argument values
  if (nargin ~= 1) error('Too few arguments.'); end;
  n = fix(max(n,0));
  
  % Allocate space for the polynomials and set H0(x) = -1
  H = zeros(n+1);
  r = 1:n;
  H(1) = -1;
  
  % Calculate the polynomials starting from H2(x)
  for i=2:n+1
      H(i,2:n+1) = H(i,2:n+1) + H(i-1,1:n)*(-1)*2; % (**)
      H(i,1:n)   = H(i,1:n)   + H(i-1,2:n+1).*r;
  end
  
  % Return results
  p=fliplr(H(n+1,:)).*(-1)^(n+1);
  
end
