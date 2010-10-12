function [mu,S,C] = cd_transform(m,P,g,g_param,tr_param)
%CD_TRANSFORM  Central Difference transform of random variables
%
% Syntax:
%   [mu,S,C,SX,W] = CD_TRANSFORM(M,P,g,n,param)
%
% In:
%   M - Random variable mean (Nx1 column vector)
%   P - Random variable covariance (NxN pos.def. matrix)
%   g - Transformation function of the form g(x,param) as
%       matrix, inline function, function name or function reference
%   g_param - Parameters of g               (optional, default empty)
%   tr_param -  Parameters of the transformation in form:
%       h = tr_param{1}: Step length of central difference approximation 
%
% Out:
%   mu - Estimated mean of y
%    S - Estimated covariance of y
%    C - Estimated cross-covariance of x and y

% Copyright (c), 2009 Jouni Hartikainen
    if nargin < 4
        g_param = [];
    end
    
    if nargin < 5
        tr_param = [];
    end
    
    if isempty(tr_param)
        h = sqrt(3);
    else
        h = tr_param{1};
    end
    
    d = size(m,1);
    cholP = chol(P)';
    
    g_mu = feval(g,m,g_param);
    s = size(g_mu,1);
    
    % Compute the first and second order terms
    a = zeros(s,d);
    H = zeros(s,d);
    for i = 1:d
        e = zeros(d,1);
        e(i) = 1;
        f1 = feval(g,cholP*(h*e)+m,g_param);
        f2 = feval(g,cholP*(-h*e)+m,g_param);
        a(:,i) = (f1 - f2) / (2*h);
        H(:,i) = (f1 + f2 - 2*g_mu) / h^2;
    end
    
    % Transformed mean
    mu = g_mu + 0.5*sum(H,2);
    
    % Covariance of y
    S = zeros(s,s);
    for i = 1:d
        S = S + a(:,i)*a(:,i)' + 0.5*H(:,i)*H(:,i)';
    end
    
    % Cross-covariance of x and y
    C = cholP*a';
      