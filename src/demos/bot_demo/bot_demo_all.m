%% Bearings Only Tracking (BOT) demonstration with various filters
%
%  Description:
%    In this example various Kalman filters and Rauch-Tung-Striebel
%    smoothers are used to estimate the position and velocity of a
%    moving object on a plane. Two sensors track the position of the 
%    object by returning noisy measurements of angular direction of the 
%    target. The methods used are:
%      * Extended Kalman filter (1st and 2nd degree)
%      * Unscented Kalman filter
%      * Gauss-Hermite Kalman filter (degree 3)
%      * Cubature Kalman filter
%    Additionally, the corresponding Rauch-Tung-Striebel smoother results 
%    are also presented.
%
%  References:
%    Refer to the Toolbox documentation for details on the model.
%
%  See also:
%    ukf_predict1, ukf_update1, urts_smooth1,
%    ekf_predict1, ekf_update1, erts_smooth1, ekf_predict2, ekf_update2,
%    ghkf_predict, ghkf_update, ghrts_smooth, 
%    ckf_predict, ckf_update, crts_smooth
%
%  Author:
%    Copyright (C) 2002, 2003 Simo Särkkä
%                  2007       Jouni Hartikainen
%                  2010       Arno Solin
%
%  Licence:
%    This software is distributed under the GNU General Public 
%    Licence (version 2 or later); please refer to the file 
%    Licence.txt, included with the software, for details.

%% Set parameters

  silent = 0;
  save_plots = 0;
  
  
%% Simulate trajectory
  
  % Measurement mean and derivative
  %
  %  h = atan((y-sy) / (x-sx))
  h_func = @bot_h;
  dh_dx_func = @bot_dh_dx;
  d2h_dx2_func = @bot_d2h_dx2;
  
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

  
%% EKF1 and EKF2
  
  clc;  clf;
  disp(['In this demonstration we track a moving object with two sensors, ',...
        'which gives only bearings of the object with respect to sensors position. ',...
       'The state of the system is estimated with 1st and 2nd order EKF, and UKF.'])
  disp(' ');
  fprintf('Running 1st order EKF...')
  % 
  % Initialize EKF1
  %
  M = M_0;
  P = P_0;
  R = sd^2;
  %if ~silent
  %  der_check(h_func, dh_dx_func, 1, M, S1);
  %  der_check(h_func, dh_dx_func, 1, M, S2);
  %end
  
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

  fprintf('Done!\n')
  fprintf('Running smoothers...');
  % ERTS smoother
  [SM1_EKF1,SP1_EKF1] = erts_smooth1(MM_EKF1,PP_EKF1,A,Q);
  eks1_rmse1 = sqrt(mean((X(1,:)-SM1_EKF1(1,:)).^2+(X(2,:)-SM1_EKF1(2,:)).^2));
  ME1_EKF1 = squeeze(SP1_EKF1(1,1,:)+SP1_EKF1(2,2,:));

  % ETF smoother
  [SM2_EKF1,SP2_EKF1] = etf_smooth1(MM_EKF1,PP_EKF1,Y,A,Q,[],[],[],...
		        dh_dx_func,R*eye(2),h_func,[],[S1 S2]);
  eks1_rmse2 = sqrt(mean((X(1,:)-SM2_EKF1(1,:)).^2+(X(2,:)-SM2_EKF1(2,:)).^2));
  ME2_EKF1 = squeeze(SP2_EKF1(1,1,:)+SP2_EKF1(2,2,:));
  fprintf('Done!\n');
  
  % Plot the results 
  if ~silent
    plot(X(1,:),X(2,:),'k-',...
         MM_EKF1(1,:),MM_EKF1(2,:),'b--',...
         SM1_EKF1(1,:),SM1_EKF1(2,:),'r-.',...
         SM2_EKF1(1,:),SM2_EKF1(2,:),'g-.',...
         S1(1),S1(2),'k^',S2(1),S2(2),'k^');
    axis([-1.5 1.5 -2.5 1.5]);

    legend('Real trajectory',...
           'EKF1 estimate',...
           'ERTS estimate',...
           'ETF estimate',...
           'Positions of sensors',...
           'Location', 'NorthWest');
    title('Filtering and smoothing result with 1st order EKF');
    if save_plots
      print -dpsc bot_demo_ekf1.ps
    end
  end
  
  disp(' ');
  disp('1st order filtering and smoothing results are now displayed');
  disp(' ');
  disp('<push any key to filter and smooth with 2nd order EKF>');
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
  
  clc; 
  fprintf('Running 2nd order EKF...');
  % Filter with EKF2
  for k = 1:size(Y,2) 
    [M,P] = ekf_predict1(M,P,A,Q);
    [M,P] = ekf_update2(M,P,Y(:,k),dh_dx_func,d2h_dx2_func,R*eye(2),h_func,[],[S1 S2]);
    MM_EKF2(:,k)   = M;
    PP_EKF2(:,:,k) = P;
    ME_EKF2(k) = P(1,1) + P(2,2);
  end
  ekf2_rmse = sqrt(mean((X(1,:)-MM_EKF2(1,:)).^2+(X(2,:)-MM_EKF2(2,:)).^2));
  fprintf('Done!\n');
  
  fprintf('Running smoothers...');
  %
  % Smoother 1
  %
  [SM1_EKF2,SP1_EKF2] = erts_smooth1(MM_EKF2,PP_EKF2,A,Q);
  eks2_rmse1 = sqrt(mean((X(1,:)-SM1_EKF2(1,:)).^2+(X(2,:)-SM1_EKF2(2,:)).^2));
  ME1_EKF2 = squeeze(SP1_EKF2(1,1,:)+SP1_EKF2(2,2,:));

  %
  % Smoother 2
  %
  [SM2_EKF2,SP2_EKF2] = etf_smooth1(MM_EKF2,PP_EKF2,Y,A,Q,[],[],[],...
		        dh_dx_func,R*eye(2),h_func,[],[S1 S2]);
  eks2_rmse2 = sqrt(mean((X(1,:)-SM2_EKF2(1,:)).^2+(X(2,:)-SM2_EKF2(2,:)).^2));
  ME2_EKF2 = squeeze(SP2_EKF2(1,1,:)+SP2_EKF2(2,2,:));
  fprintf('Done!\n');
  
  if ~silent
    plot(X(1,:),X(2,:),'k-',...
         MM_EKF2(1,:),MM_EKF2(2,:),'b--',...
         SM1_EKF2(1,:),SM1_EKF2(2,:),'r-.',...
         SM2_EKF2(1,:),SM2_EKF2(2,:),'g-.',...
         S1(1),S1(2),'k^',S2(1),S2(2),'k^');
    axis([-1.5 1.5 -2.5 1.5]);
    legend('Real trajectory',...
           'EKF2 estimate',...
           'ERTS estimate',...
           'ETF estimate',...
           'Positions of sensors',...
           'Location', 'NorthWest');
    title('Filtering and smoothing result with 2nd order EKF');
    if save_plots
      print -dpsc bot_demo_ekf2.ps
    end
  end
  
  disp(' ');
  disp('2nd order filtering and smoothing results are now displayed.');
  disp(' ');

  
%% UKF
  
  disp('<push any key to filter and smooth with unscented filter and smoothers>')
  pause
  
  clc;
  fprintf('Running UKF...');
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
  
  % Calculate RMSE of UKF
  ukf_rmse = sqrt(mean((X(1,:)-MM_UKF(1,:)).^2+(X(2,:)-MM_UKF(2,:)).^2));

  fprintf('Done!\n');

  fprintf('Running smoothers...');  
  % URTS Smoother   
  [SM1_UKF,SP1_UKF] = urts_smooth1(MM_UKF,PP_UKF,A,Q);
  uks_rmse1 = sqrt(mean((X(1,:)-SM1_UKF(1,:)).^2+(X(2,:)-SM1_UKF(2,:)).^2));
  ME1_UKF = squeeze(SP1_UKF(1,1,:)+SP1_UKF(2,2,:));
  
  % UTF Smoother
  IAW = inv(A)*[eye(size(A,1)) eye(size(A,1))];
  [SM2_UKF,SP2_UKF] = utf_smooth1(MM_UKF,PP_UKF,Y,IAW,Q,[],...
		                   h_func,R*eye(2),[S1 S2]);
  uks_rmse2 = sqrt(mean((X(1,:)-SM2_UKF(1,:)).^2+(X(2,:)-SM2_UKF(2,:)).^2));
  ME2_UKF = squeeze(SP2_UKF(1,1,:)+SP2_UKF(2,2,:));

  fprintf('Done!\n\n');

  % Plot UKF, URTS, UTF
  if ~silent
    plot(X(1,:),X(2,:),'k-',...
         MM_UKF(1,:),MM_UKF(2,:),'b--',...
         SM1_UKF(1,:),SM1_UKF(2,:),'r-.',...
         SM2_UKF(1,:),SM2_UKF(2,:),'g-.',...
         S1(1),S1(2),'k^',S2(1),S2(2),'k^');
    axis([-1.5 1.5 -2.5 1.5]);
    legend('Real trajectory',...
           'UKF estimate',...
           'URTS estimate',...
           'UTF estimate',...
           'Positions of sensors',...
           'Location', 'NorthWest');
    title('Filtering and smoothing result with UKF');
    if save_plots
      print -dpsc bot_demo_ukf.ps
    end
  end
  
  disp(' ');
  disp('Unscented filtering and smoothing results are now displayed.')
  disp(' ');
  
  
  %% GHKF
  
  disp('<push any key to filter and smooth with Gauss-Hermite Kalman filter and smoother>')
  pause
  
  clc;
  fprintf('Running GHKF...');
  % Initialize GHKF
  M = M_0;
  P = P_0;
  MM_GHKF = zeros(size(M,1),size(Y,2));
  PP_GHKF = zeros(size(M,1),size(M,1),size(Y,2));
  ME_GHKF = zeros(size(M,1),1);
  
  % Filter with GHKF
  for k=1:size(Y,2)
    [M,P] = ghkf_predict(M,P,A,Q,[],3);
    [M,P] = ghkf_update(M,P,Y(:,k),h_func,R*eye(2),[S1 S2],3);
    MM_GHKF(:,k)   = M;
    PP_GHKF(:,:,k) = P;
    ME_GHKF(k) = P(1,1) + P(2,2);
  end
  
  % Calculate RMSE of GHKF
  ghkf_rmse = sqrt(mean((X(1,:)-MM_GHKF(1,:)).^2+(X(2,:)-MM_GHKF(2,:)).^2));

  fprintf('Done!\n');

  fprintf('Running smoothers...');  
  % GHRTS Smoother   
  [SM1_GHKF,SP1_GHKF] = ghrts_smooth(MM_GHKF,PP_GHKF,A,Q,[],3);
  ghks_rmse1 = sqrt(mean((X(1,:)-SM1_GHKF(1,:)).^2+(X(2,:)-SM1_GHKF(2,:)).^2));
  ME1_GHKF = squeeze(SP1_GHKF(1,1,:)+SP1_GHKF(2,2,:));
  
  fprintf('Done!\n\n');
  
  % Plot GHKF, GHRTS
  if ~silent
    plot(X(1,:),X(2,:),'k-',...
         MM_GHKF(1,:),MM_GHKF(2,:),'b--',...
         SM1_GHKF(1,:),SM1_GHKF(2,:),'r-.',...
         S1(1),S1(2),'k^',S2(1),S2(2),'k^');
    axis([-1.5 1.5 -2.5 1.5]);
    legend('Real trajectory',...
           'GHKF estimate',...
           'GHRTS estimate',...
           'Positions of sensors',...
           'Location', 'NorthWest');
    title('Filtering and smoothing result with GHKF');
    if save_plots
      print -dpsc bot_demo_ghkf.ps
    end
  end
  
  disp(' ');
  disp('Gauss-Hermite filtering and smoothing results are now displayed.')
  disp(' ');
  
  

  %% CKF
  
  disp('<push any key to filter and smooth with Cubature Kalman filter and smoother>')
  pause
  
  clc;
  fprintf('Running CKF...');
  % Initialize CKF
  M = M_0;
  P = P_0;
  MM_CKF = zeros(size(M,1),size(Y,2));
  PP_CKF = zeros(size(M,1),size(M,1),size(Y,2));
  ME_CKF = zeros(size(M,1),1);
  
  % Filter with CKF
  for k=1:size(Y,2)
    [M,P] = ckf_predict(M,P,A,Q,[]);
    [M,P] = ckf_update(M,P,Y(:,k),h_func,R*eye(2),[S1 S2]);
    MM_CKF(:,k)   = M;
    PP_CKF(:,:,k) = P;
    ME_CKF(k) = P(1,1) + P(2,2);
  end
  
  % Calculate RMSE of CKF
  ckf_rmse = sqrt(mean((X(1,:)-MM_CKF(1,:)).^2+(X(2,:)-MM_CKF(2,:)).^2));

  fprintf('Done!\n');

  fprintf('Running smoothers...');  
  % CRTS Smoother   
  [SM1_CKF,SP1_CKF] = crts_smooth(MM_CKF,PP_CKF,A,Q,[]);
  cks_rmse1 = sqrt(mean((X(1,:)-SM1_CKF(1,:)).^2+(X(2,:)-SM1_CKF(2,:)).^2));
  ME1_CKF = squeeze(SP1_CKF(1,1,:)+SP1_CKF(2,2,:));
  
  fprintf('Done!\n\n');

  % Plot CKF, CRTS
  if ~silent
    plot(X(1,:),X(2,:),'k-',...
         MM_CKF(1,:),MM_CKF(2,:),'b--',...
         SM1_CKF(1,:),SM1_CKF(2,:),'r-.',...
         S1(1),S1(2),'k^',S2(1),S2(2),'k^');
    axis([-1.5 1.5 -2.5 1.5]);
    legend('Real trajectory',...
           'CKF estimate',...
           'CRTS estimate',...
           'Positions of sensors',...
           'Location', 'NorthWest');
    title('Filtering and smoothing result with CKF');
    if save_plots
      print -dpsc bot_demo_ckf.ps
    end
  end
  
  disp(' ');
  disp('Cubature filtering and smoothing results are now displayed.')
  disp(' ');
  
  
%% Show RMSE values for all methods
  
  % Print errors
  disp('RMS errors:')
  fprintf('EKF1-RMSE  = %.4f [%.4f]\n',ekf1_rmse,  sqrt(mean(ME_EKF1)));
  fprintf('ERTS1-RMSE = %.4f [%.4f]\n',eks1_rmse1, sqrt(mean(ME1_EKF1)));
  fprintf('ETF1-RMSE  = %.4f [%.4f]\n',eks1_rmse2, sqrt(mean(ME2_EKF1)));
  fprintf('EKF2-RMSE  = %.4f [%.4f]\n',ekf2_rmse,  sqrt(mean(ME_EKF2)));
  fprintf('ERTS2-RMSE = %.4f [%.4f]\n',eks2_rmse1, sqrt(mean(ME1_EKF2)));
  fprintf('ETF2-RMSE  = %.4f [%.4f]\n',eks2_rmse2, sqrt(mean(ME2_EKF2)));  
  fprintf('UKF-RMSE   = %.4f [%.4f]\n',ukf_rmse,   sqrt(mean(ME_UKF)));
  fprintf('URTS-RMSE  = %.4f [%.4f]\n',uks_rmse1,  sqrt(mean(ME1_UKF)));    
  fprintf('UTF-RMSE   = %.4f [%.4f]\n',uks_rmse2,  sqrt(mean(ME2_UKF)));
  fprintf('GHKF-RMSE  = %.4f [%.4f]\n',ghkf_rmse,   sqrt(mean(ME_GHKF)));
  fprintf('GHRTS-RMSE = %.4f [%.4f]\n',ghks_rmse1,  sqrt(mean(ME1_GHKF)));    
  fprintf('CKF-RMSE   = %.4f [%.4f]\n',ckf_rmse,   sqrt(mean(ME_CKF)));
  fprintf('CRTS-RMSE  = %.4f [%.4f]\n',cks_rmse1,  sqrt(mean(ME1_CKF))); 
  