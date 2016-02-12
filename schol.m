%SCHOL  Cholesky factorization for positive semidefinite matrices
%
% Syntax:
%   [L,def] = schol(A)
%
% In:
%   A - Symmetric pos.semi.def matrix to be factorized
%
% Out:
%   L   - Lower triangular matrix such that A=L*L' if def>=0.
%   def - Value 1,0,-1 denoting that A was positive definite,
%         positive semidefinite or negative definite, respectively.
%
% Description:
%   Compute lower triangular Cholesky factor L of symmetric positive
%   semidefinite matrix A such that
%
%     A = L*L'
%
% See also
%   CHOL

% Copyright (C) 2006 Simo Särkkä
%
% $Id$
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [L,def] = schol(A)

  L  = zeros(size(A));
  def = 1;
  for i=1:size(A,1)
    for j=1:i
      s = A(i,j);
      for k=1:j-1
	s = s - L(i,k)*L(j,k);
      end
      if j < i
	if L(j,j) > eps
          L(i,j) = s / L(j,j);
        else
          L(i,j) = 0;
        end
      else
	if (s < -eps)
	  s = 0;
	  def = -1;        
	elseif (s < eps)
	  s = 0;
	  def = min(0,def);
	end
	L(j,j) = sqrt(s);
      end
    end  
  end
  
  if (nargout < 2) & (def < 0)
    warning('Matrix is negative definite !!!!');
  end

