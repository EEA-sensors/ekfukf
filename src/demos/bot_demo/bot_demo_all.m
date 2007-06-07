% Bearings Only Tracking (BOT demonstration with EKF1,EKF2 and UKF

% Copyright (C) 2002, 2003 Simo Särkkä
%                     2007 Jouni Hartikainen
%
% $Id: ukfs_demo1.m,v 1.5 2006/10/10 20:18:54 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

  silent = 0;
  save_plots = 1;
  
  % Measurement mean and derivative
  %
  %  h = atan((y-sy) / (x-sx))
  h_func = @az_h;
  dh_dx_func = @az_dh_dx;
  d2h_dx2_func = @az_d2h_dx2;
  
  % Create a bit curved trajectory and angle
  % measurements from two sensors
  S1 = [-1;-2];
  S2 = [1;1];
  sd = 0.05;
  dt = 0.01;

  % Acceleration for the target to have a curved trajectory  
  a = zeros(1,500);
  a(1,50:100)  = pi/2/51/dt + 0.01*randn(1,51);
  a(1,200:250) = pi/2/51/dt + 0.01*randn(1,51);
  a(1,350:400) = pi/2/51/dt + 0.01*randn(1,51);
    
  % Starting state
  x = [0;0;1;0];
  t = 0;
  X = [];
  Y = [];
  T = [];
  for i=1:500
    % Dynamics in continous case
    F = [0 0  1    0;...
         0 0  0    1;...
         0 0  0   a(i);...
         0 0 -a(i) 0];
    % Discretization for the data generation
    x = expm(F*dt)*x;
    % Angular measurements for both sensors.
    y1 = atan2(x(2)-S1(2), x(1)-S1(1)) + sd * randn;
    y2 = atan2(x(2)-S2(2), x(1)-S2(1)) + sd * randn;
      
    % Save the data
    t  = t + dt;
    X = [X x];
    T = [T t];
    Y = [Y [y1;y2]];
  end

  % Prior for position and velocity
  M_0 = [0;0;0;0];
  P_0 = diag([0.1 0.1 10 10]);

  % Discretize the continous model
  qx = 0.1;
  qy = 0.1;
  F = [0 0 1 0;
       0 0 0 1;
       0 0 0 0;
       0 0 0 0];
  [A,Q] = lti_disc(F,[],diag([0 0 qx qy]),dt);


  % 
  % Initialize EKF1
  %
  M = M_0;
  P = P_0;
  R = sd^2;
  if ~silent
    der_check(h_func, dh_dx_func, 1, M, S1);
    der_check(h_func, dh_dx_func, 1, M, S2);
  end
  
  MM_EKF1 = zeros(size(M,1),size(Y,2));
  PP_EKF1 = zeros(size(M,1),size(M,1),size(Y,2));
  ME_EKF1 = zeros(size(M,1),1);

  % Filter with EKF1
  for k = 1:size(Y,2) 
    [M,P] = ekf_predict1(M,P,A,Q);
    [M,P] = ekf_update1(M,P,Y(:,k),dh_dx_func,R*eye(2),h_func,[],[S1 S2]);
    MM_EKF1(:,k)   = M;
    PP_EKF1(:,:,k) = P;
    ME_EKF1(k) = P(1,1) + P(2,2);
  end
  
  % RMSE for EKF1
  ekf1_rmse = sqrt(mean((X(1,:)-MM_EKF1(1,:)).^2+(X(2,:)-MM_EKF1(2,:)).^2));
  fprintf('EKF1-RMSE  = %.3f  [%.3f]\n',ekf1_rmse,sqrt(mean(ME_EKF1)));

  % ERTS smoother
  [SM1_EKF1,SP1_EKF1] = erts_smooth1(MM_EKF1,PP_EKF1,A,Q);
  eks1_rmse1 = sqrt(mean((X(1,:)-SM1_EKF1(1,:)).^2+(X(2,:)-SM1_EKF1(2,:)).^2));
  ME1_EKF1 = squeeze(SP1_EKF1(1,1,:)+SP1_EKF1(2,2,:));
  fprintf('EKS1-RMSE1 = %.4f [%.4f]\n',eks1_rmse1,sqrt(mean(ME1_EKF1)));

  % EFBF smoother
  [SM2_EKF1,SP2_EKF1] = efbf_smooth1(MM_EKF1,PP_EKF1,Y,A,Q,[],[],[],...
		        dh_dx_func,R*eye(2),h_func,[],[S1 S2]);
  eks1_rmse2 = sqrt(mean((X(1,:)-SM2_EKF1(1,:)).^2+(X(2,:)-SM2_EKF1(2,:)).^2));
  ME2_EKF1 = squeeze(SP2_EKF1(1,1,:)+SP2_EKF1(2,2,:));
  fprintf('EKS-RMSE2 = %.4f [%.4f]\n',eks1_rmse2,sqrt(mean(ME2_EKF1)));
  
  % Plot the results 
  if ~silent
    plot(X(1,:),X(2,:),'k-',...
         MM_EKF1(1,:),MM_EKF1(2,:),'b--',...
         SM1_EKF1(1,:),SM1_EKF1(2,:),'r-.',...
         SM2_EKF1(1,:),SM2_EKF1(2,:),'g-.',...
         S1(1),S1(2),'k^',S2(1),S2(2),'k^');
    axis([-1.5 1.5 -2.5 1]);
    legend('Real trajectory',...
           'EKF1 estimate',...
           'ERTS estimate',...
           'EFBF estimate',...
           'Positions of sensors',...
           'Location', 'NorthWest');
    title('Filtering and smoothing result with 1st order EKF');
    if save_plots
      print -dpsc bot_demo_ekf1.ps
    end
  end
  
  pause
  
  %
  % Initialize EKF2
  %
  M = M_0; 
  P = P_0;
  R = sd^2;
  if ~silent
      %der_check(dh_dx_func, d2h_dx2_func, 1, M, S1);
      %der_check(dh_dx_func, d2h_dx2_func, 1, M, S2);
  end
  
  MM_EKF2 = zeros(size(M,1),size(Y,2));
  PP_EKF2 = zeros(size(M,1),size(M,1),size(Y,2));
  ME_EKF2 = zeros(size(M,1),1);
  
  % Filter with EKF2
  for k = 1:size(Y,2) 
    [M,P] = ekf_predict1(M,P,A,Q);
    [M,P] = ekf_update2(M,P,Y(:,k),dh_dx_func,d2h_dx2_func,R*eye(2),h_func,[],[S1 S2]);
    MM_EKF2(:,k)   = M;
    PP_EKF2(:,:,k) = P;
    ME_EKF2(k) = P(1,1) + P(2,2);
  end
  
  ekf2_rmse = sqrt(mean((X(1,:)-MM_EKF2(1,:)).^2+(X(2,:)-MM_EKF2(2,:)).^2));
  fprintf('EKF2-RMSE  = %.3f  [%.3f]\n',ekf2_rmse,sqrt(mean(ME_EKF2)));

  %
  % Smoother 1
  %
  [SM1_EKF2,SP1_EKF2] = erts_smooth1(MM_EKF2,PP_EKF2,A,Q);
  eks2_rmse1 = sqrt(mean((X(1,:)-SM1_EKF2(1,:)).^2+(X(2,:)-SM1_EKF2(2,:)).^2));
  ME1_EKF2 = squeeze(SP1_EKF2(1,1,:)+SP1_EKF2(2,2,:));
  fprintf('EKS2-RMSE1 = %.4f [%.4f]\n',eks2_rmse1,sqrt(mean(ME1_EKF2)));

  %
  % Smoother 2
  %
  [SM2_EKF2,SP2_EKF2] = efbf_smooth1(MM_EKF2,PP_EKF2,Y,A,Q,[],[],[],...
		        dh_dx_func,R*eye(2),h_func,[],[S1 S2]);
  eks2_rmse2 = sqrt(mean((X(1,:)-SM2_EKF2(1,:)).^2+(X(2,:)-SM2_EKF2(2,:)).^2));
  ME2_EKF2 = squeeze(SP2_EKF2(1,1,:)+SP2_EKF2(2,2,:));
  fprintf('EKS2-RMSE2 = %.4f [%.4f]\n',eks2_rmse2,sqrt(mean(ME2_EKF2)));
  
  if ~silent
    plot(X(1,:),X(2,:),'k-',...
         MM_EKF2(1,:),MM_EKF2(2,:),'b--',...
         SM1_EKF2(1,:),SM1_EKF2(2,:),'r-.',...
         SM2_EKF2(1,:),SM2_EKF2(2,:),'g-.',...
         S1(1),S1(2),'k^',S2(1),S2(2),'k^');
    axis([-1.5 1.5 -2.5 1]);
    legend('Real trajectory',...
           'EKF2 estimate',...
           'ERTS estimate',...
           'EFBF estimate',...
           'Positions of sensors',...
           'Location', 'NorthWest');
    title('Filtering and smoothing result with 2nd order EKF');
    if save_plots
      print -dpsc bot_demo_ekf2.ps
    end
  end
  
  pause
  
  % Initialize UKF
  M = M_0;
  P = P_0;
  MM_UKF = zeros(size(M,1),size(Y,2));
  PP_UKF = zeros(size(M,1),size(M,1),size(Y,2));
  ME_UKF = zeros(size(M,1),1);
  
  % Filter with UKF
  for k=1:size(Y,2)
    [M,P] = ukf_predict1(M,P,A,Q);
    [M,P] = ukf_update1(M,P,Y(:,k),h_func,R*eye(2),[S1 S2]);
    MM_UKF(:,k)   = M;
    PP_UKF(:,:,k) = P;
    ME_UKF(k) = P(1,1) + P(2,2);
  end
  
  % Calculate RMSE
  ukf_rmse = sqrt(mean((X(1,:)-MM_UKF(1,:)).^2+(X(2,:)-MM_UKF(2,:)).^2));
  fprintf('UKF-RMSE  = %.3f  [%.3f]\n',ukf_rmse,sqrt(mean(ME_UKF)));

  % URTS Smoother   
  [SM1_UKF,SP1_UKF] = urts_smooth1(MM_UKF,PP_UKF,A,Q);
  uks_rmse1 = sqrt(mean((X(1,:)-SM1_UKF(1,:)).^2+(X(2,:)-SM1_UKF(2,:)).^2));
  ME1_UKF = squeeze(SP1_UKF(1,1,:)+SP1_UKF(2,2,:));
  fprintf('UKS-RMSE1 = %.4f [%.4f]\n',uks_rmse1,sqrt(mean(ME1_UKF)));  
  
  % UFBF Smoother
  IAW = inv(A)*[eye(size(A,1)) eye(size(A,1))];
  [SM2_UKF,SP2_UKF] = ufbf_smooth1(MM_UKF,PP_UKF,Y,IAW,Q,[],...
		                   h_func,R*eye(2),[S1 S2]);

  uks_rmse2 = sqrt(mean((X(1,:)-SM2_UKF(1,:)).^2+(X(2,:)-SM2_UKF(2,:)).^2));
  ME2_UKF = squeeze(SP2_UKF(1,1,:)+SP2_UKF(2,2,:));
  fprintf('UKS-RMSE2 = %.4f [%.4f]\n',uks_rmse2,sqrt(mean(ME2_UKF)));
  
  % Plot
  if ~silent
    plot(X(1,:),X(2,:),'k-',...
         MM_UKF(1,:),MM_UKF(2,:),'b--',...
         SM1_UKF(1,:),SM1_UKF(2,:),'r-.',...
         SM2_UKF(1,:),SM2_UKF(2,:),'g-.',...
         S1(1),S1(2),'k^',S2(1),S2(2),'k^');
    axis([-1.5 1.5 -2.5 1]);
    legend('Real trajectory',...
           'UKF estimate',...
           'URTS estimate',...
           'UFBF estimate',...
           'Positions of sensors',...
           'Location', 'NorthWest');
    title('Filtering and smoothing result with UKF');
    if save_plots
      print -dpsc bot_demo_ukf.ps
    end
  end