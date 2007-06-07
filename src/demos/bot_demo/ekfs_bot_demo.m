%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Extended Kalman Filter / Smoother demonstration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2002 Simo Särkkä
%
% $Id: ekfs_demo1.m,v 1.6 2006/10/10 20:18:51 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

  keep_trajectory = 0;
  silent = 0;
  
  %
  % Measurement mean and derivative
  %
  h_func = @az_h;
  dh_dx_func = @az_dh_dx;
  d2h_dx2_func = @az_d2h_dx2;
  %
  % Create a bit curved trajectory and angle
  % measurements from two sensors
  %
  S1 = [-1;-2];
  S2 = [1;1];
  sd = 0.05;
  dt = 0.01;
  if ~keep_trajectory
    a = zeros(1,500);
    a(1,50:100)  = pi/2/51/dt + 0.01*randn(1,51);
    a(1,200:250) = pi/2/51/dt + 0.01*randn(1,51);
    a(1,350:400) = pi/2/51/dt + 0.01*randn(1,51);
    x = [0;0;1;0];
    t = 0;
    X = [];
    Y = [];
    T = [];
    for i=1:500
      F = [0 0  1    0;...
           0 0  0    1;...
           0 0  0   a(i);...
           0 0 -a(i) 0];
      x = expm(F*dt)*x;
      y1 = atan2(x(2)-S1(2), x(1)-S1(1)) + sd * randn;
      y2 = atan2(x(2)-S2(2), x(1)-S2(1)) + sd * randn;
      t  = t + dt;
      X = [X x];
      T = [T t];
      Y = [Y [y1;y2]];
    end
  end
  
  plot(1:size(Y,2),Y(1,:),'.',1:size(Y,2),Y(2,:),'.');
  legend('Sensor 1', 'Sensor 2');
  title('Measurements from sensors in radians');
  print -dpsc bot_demo_measurements.ps

  %
  % Initialize EKF to values
  %
  %   x = 0
  %   y = 0,
  %   dx/dt = 0
  %   dy/dt = 0
  %
  % with great uncertainty in velocity
  %
  M = [0;0;0;0];
  P = diag([4 4 4 4]);
  R = sd^2;
  if ~silent
    der_check(h_func, dh_dx_func, 1, M, S1);
    der_check(h_func, dh_dx_func, 1, M, S2);
  end
  qx = 0.1;
  qy = 0.1;
  F = [0 0 1 0;
       0 0 0 1;
       0 0 0 0;
       0 0 0 0];
  [A,Q] = lti_disc(F,[],diag([0 0 qx qy]),dt);

  %
  % Track and animate
  %
  MM = zeros(size(M,1),size(Y,2));
  PP = zeros(size(M,1),size(M,1),size(Y,2));
  ME = zeros(size(M,1),1);
  for k=1:size(Y,2)
    %
    % Track with EKF
    %
    [M,P] = ekf_predict1(M,P,A,Q);

    [M,P] = ekf_update1(M,P,Y(:,k),dh_dx_func,R*eye(2),h_func,[],[S1 S2]);
    %[M,P] = ekf_update2(M,P,Y(:,k),dh_dx_func,d2h_dx2_func,R*eye(2),h_func,[],[S1 S2]);
    
    MM(:,k)   = M;
    PP(:,:,k) = P;
    ME(k) = P(1,1) + P(2,2);
    
    %
    % Confidence ellipse
    %
    tt = (0:0.01:1)*2*pi;
    cc = repmat(M(1:2),1,length(tt)) + ...
	 2*chol(P(1:2,1:2))'*[cos(tt);sin(tt)];

    %
    % Animate
    %
    if ~silent
      if rem(k,10) == 1
        len = 2.5;
        dx1 = len*cos(Y(1,k));
        dy1 = len*sin(Y(1,k));
        dx2 = len*cos(Y(2,k));
        dy2 = len*sin(Y(2,k));
        plot(X(1,:),X(2,:),'r-',...
             M(1),M(2),'bo',...
             MM(1,1:k),MM(2,1:k),'b--',...
    	     cc(1,:),cc(2,:),'g-',...
             [S1(1);S1(1)+dx1],[S1(2);S1(2)+dy1],'k--',...
             [S2(1);S2(1)+dx2],[S2(2);S2(2)+dy2],'k--',...
             S1(1),S1(2),'k^',S2(1),S2(2),'k^');
        axis([-1.5 1.5 -2.5 1]);
        drawnow;
        pause;
      end
    end
  end

 
  
  %
  % Calculate RMSE
  %
  ekf_rmse = sqrt(mean((X(1,:)-MM(1,:)).^2+(X(2,:)-MM(2,:)).^2));
  fprintf('EKF-RMSE  = %.3f  [%.3f]\n',ekf_rmse,sqrt(mean(ME)));

  %
  % Smoother 1
  %
  [SM1,SP1] = erts_smooth1(MM,PP,A,Q);
  eks_rmse1 = sqrt(mean((X(1,:)-SM1(1,:)).^2+(X(2,:)-SM1(2,:)).^2));
  ME1 = squeeze(SP1(1,1,:)+SP1(2,2,:));
  fprintf('EKS-RMSE1 = %.4f [%.4f]\n',eks_rmse1,sqrt(mean(ME1)));

  %
  % Smoother 2
  %
  [SM2,SP2] = efbf_smooth1(MM,PP,Y,A,Q,[],[],[],...
		     dh_dx_func,R*eye(2),h_func,[],[S1 S2]);
  eks_rmse2 = sqrt(mean((X(1,:)-SM2(1,:)).^2+(X(2,:)-SM2(2,:)).^2));
  ME2 = squeeze(SP2(1,1,:)+SP2(2,2,:));
  fprintf('EKS-RMSE2 = %.4f [%.4f]\n',eks_rmse2,sqrt(mean(ME2)));

  EST = [];
  for k=size(Y,2):-1:1
    M = SM1(:,k);  
    P = SP1(:,:,k);
    EST = [EST M];

    % Confidence ellipse
    tt = (0:0.01:1)*2*pi;
    cc = repmat(M(1:2),1,length(tt)) + ...
	 2*chol(P(1:2,1:2))'*[cos(tt);sin(tt)];

    % Animate
    len = 1.5;
    if ~silent
        if rem(k,10) == 1
            plot(X(1,:),X(2,:),'r-',...
            M(1),M(2),'o',...
            EST(1,:),EST(2,:),'--',...
	    cc(1,:),cc(2,:),'g-',...
            MM(1,:),MM(2,:),'b--',...
            S1(1),S1(2),'k^',S2(1),S2(2),'k^');
            %Y(1,:),Y(2,:),'.',...         
            drawnow;
            pause;
        end
     end
  end
  

  %
  % Plot
  %
  if ~silent
    plot(X(1,:),X(2,:),'k-',...
         MM(1,:),MM(2,:),'b--',...
         SM1(1,:),SM1(2,:),'r-.',...
         SM2(1,:),SM2(2,:),'g-.',...
         S1(1),S1(2),'k^',S2(1),S2(2),'k^');
    axis([-1.5 1.5 -2.5 1]);
    legend('Real trajectory',...
           'EKF1 estimate',...
           'ERTS estimate',...
           'EFBF estimate',...
           'Positions of sensors',...
           'Location', 'NorthWest');
    title('Filtering and smoothing result with 1st order EKF');
    print -dpsc bot_demo_ekf1.ps
  end
  
  EMM = MM;
  ESM1 = SM1;
  ESM2 = SM2;

  
