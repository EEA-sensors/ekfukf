% Kalman Filter demonstration with sine signal.
%
% History:
%    3.12.2002 SS  The first implementation
%
% Copyright (C) 2002 Simo Särkkä
%
% $Id: kf_demo.m,v 1.1 2005/08/02 11:10:46 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

  %
  % Create sine function
  %
  S1 = [0.2;1.0];
  S2 = [1.0;-0.2];
  sd = 0.1;
  dt = 0.1;
  w = 1;
  T = (0:dt:30);
  X = sin(w*T);
  Y = X + sd*randn(size(X));

  %
  % Initialize KF to values
  %
  %   x = 0
  %   dx/dt = 0
  %
  % with great uncertainty in derivative
  %
  M = [0;0];
  P = diag([0.1 2]);
  R = sd^2;
  H = [1 0];
  q = 0.1;
  F = [0 1;
       0 0];
  [A,Q] = lti_disc(F,[],diag([0 q]),dt);

  %
  % Track and animate
  %
  MM = zeros(size(M,1),size(Y,2));
  PP = zeros(size(M,1),size(M,1),size(Y,2));
  clf;
  for k=1:size(Y,2)
    %
    % Track with KF
    %
    [M,P] = kf_predict(M,P,A,Q);
    [M,P] = kf_update(M,P,Y(k),H,R);

    MM(:,k) = M;
    PP(:,:,k) = P;

    %
    % Animate
    %
    if rem(k,10)==1
      plot(T,X,'b--',...
           T,Y,'ro',...
           T(k),M(1),'k*',...
           T(1:k),MM(1,1:k),'k-');
      drawnow;
      pause;
    end
  end

  %
  % Apply Kalman smoother
  %
  SM = rts_smooth(MM,PP,A,Q);
  plot(T,X,'b--',...
       T,MM(1,:),'k-',...
       T,SM(1,:),'r-');

  %
  % Errors
  %
  fprintf('kf_rmse = %.3f\nks_rmse = %.3f\n',...
          sqrt(mean((MM(1,:)-X(1,:)).^2)),...
          sqrt(mean((SM(1,:)-X(1,:)).^2)));
