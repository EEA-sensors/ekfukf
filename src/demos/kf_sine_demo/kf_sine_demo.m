% Kalman Filter demonstration with sine signal.
%
% History:
%    3.12.2002 SS  The first implementation
%
% Copyright (C) 2002 Simo Särkkä
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
  clc;
  disp('In this demonstration we estimate a stationary sine signal from noisy measurements by using the classical Kalman filter.');
  disp(' ');
  disp('The filtering results are now displayed sequantially for 10 time step at a time.');
  disp(' ');
  disp('<push any key to proceed to next time steps>');
  
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
      legend('Real signal','Measurements','Latest estimate','Filtered estimate')
      title('Estimating a noisy sine signal with Kalman filter.');
      drawnow;
      
      pause;
    end
  end

  clc;
  disp('In this demonstration we estimate a stationary sine signal from noisy measurements by using the classical Kalman filter.');
  disp(' ');
  disp('The filtering results are now displayed sequantially for 10 time step at a time.');
  disp(' ');
  disp('<push any key to see the filtered and smoothed results together>')
  pause;  
  %
  % Apply Kalman smoother
  %
  SM = rts_smooth(MM,PP,A,Q);
  plot(T,X,'b--',...
       T,MM(1,:),'k-',...
       T,SM(1,:),'r-');
  legend('Real signal','Filtered estimate','Smoothed estimate') 
  title('Filtered and smoothed estimate of the original signal');
  
  clc;
  disp('The filtered and smoothed estimates of the signal are now displayed.')
  disp(' ');
  disp('RMS errors:');
  %
  % Errors
  %
  fprintf('KF = %.3f\nRTS = %.3f\n',...
          sqrt(mean((MM(1,:)-X(1,:)).^2)),...
          sqrt(mean((SM(1,:)-X(1,:)).^2)));
