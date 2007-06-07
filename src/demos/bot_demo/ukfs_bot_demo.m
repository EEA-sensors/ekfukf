% Bearings Only Tracking (BOT) demonstration with UKF. 

% Copyright (C) 2002, 2003 Simo Särkkä
%
% $Id: ukfs_demo1.m,v 1.5 2006/10/10 20:18:54 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

  keep_trajectory = 0;
  silent = 0;

  %
  % Measurement mean and derivative
  %
  %  h = atan((y-sy) / (x-sx))
  %
  h_func = @az_h;

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

  %
  % Initialize UKF to values
  %
  %   x = 0
  %   y = 0,
  %   dx/dt = 0
  %   dy/dt = 0
  %
  % with great uncertainty in velocity
  %
  M = [0;0;0;0];
  P = diag([0.1 0.1 10 10]);
  R = sd^2;
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
    % Track with UKF
    %
    [M,P] = ukf_predict1(M,P,A,Q);

%    [M,P] = ukf_update1(M,P,Y(1,k),h_func,R,S1);
%    [M,P] = ukf_update1(M,P,Y(2,k),h_func,R,S2);

    [M,P] = ukf_update1(M,P,Y(:,k),h_func,R*eye(2),[S1 S2]);

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
             [S2(1);S2(1)+dx2],[S2(2);S2(2)+dy2],'k--');
        axis([-1.5 1.5 -2.5 1]);
        drawnow;
      end
    end
  end

  %
  % Calculate RMSE
  %
  ukf_rmse = sqrt(mean((X(1,:)-MM(1,:)).^2+(X(2,:)-MM(2,:)).^2));
  fprintf('UKF-RMSE  = %.3f  [%.3f]\n',ukf_rmse,sqrt(mean(ME)));

  %
  % Smoother 1
  %
  [SM1,SP1] = urts_smooth1(MM,PP,A,Q);
  uks_rmse1 = sqrt(mean((X(1,:)-SM1(1,:)).^2+(X(2,:)-SM1(2,:)).^2));
  ME1 = squeeze(SP1(1,1,:)+SP1(2,2,:));
  fprintf('UKS-RMSE1 = %.4f [%.4f]\n',uks_rmse1,sqrt(mean(ME1)));  
  
  %
  % Smoother 2
  %
  IAW = inv(A)*[eye(size(A,1)) eye(size(A,1))];
  [SM2,SP2] = ufbf_smooth1(MM,PP,Y,IAW,Q,[],...
		     h_func,R*eye(2),[S1 S2]);

  uks_rmse2 = sqrt(mean((X(1,:)-SM2(1,:)).^2+(X(2,:)-SM2(2,:)).^2));
  ME2 = squeeze(SP2(1,1,:)+SP2(2,2,:));
  fprintf('UKS-RMSE2 = %.4f [%.4f]\n',uks_rmse2,sqrt(mean(ME2)));

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
    print -dpsc bot_demo_ukf.ps
  end
  
  UMM = MM;
  USM1 = SM1;
  USM2 = SM2;

  
