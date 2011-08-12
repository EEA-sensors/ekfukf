%% Discrete-time reentry dynamics demonstration with non-linear filters
%
%  Description:
%    In this example various different non-linear filters and smoothers are
%    applied to reentry tracking problem . The filters used in this 
%    demonstration are:
%      * Extended Kalman filter
%      * Unscented Kalman filter
%      * Gauss-Hermite Kalman filter (degree 3)
%      * Cubature Kalman filter
%    Additionally, the corresponding smoother results are also presented.
%
%  References:
%    Refer to the Toolbox documentation for details on the model.
%
%  See also:
%    ukf_predict1, ukf_update1, urts_smooth1,
%    ekf_predict1, ekf_update1, erts_smooth1, ghkf_predict, ghkf_update, 
%    ghrts_smooth, ckf_predict, ckf_update, crts_smooth
%
%  Author:
%    Copyright (C) 2006 Simo Särkkä
%                  2007 Jouni Hartikainen
%                  2010 Arno Solin
%
%  Licence:
%    This software is distributed under the GNU General Public 
%    Licence (version 2 or later); please refer to the file 
%    Licence.txt, included with the software, for details.

%% Set up model parameters

  reentry_param;
  make_reentry_data;

  silent = 0;

  % Handles to dynamic and measurement model functions,
  % and to their derivatives.
  func_f  = @reentry_f;
  func_if = @reentry_if;
  func_df = @reentry_df_dx;
  func_h  = @reentry_h;
  func_dh = @reentry_dh_dx;
  
  % Initial values and space for EKF
  m = m0;
  P = P0;
  Q = L*Qc*L'*dt;
  
  MM_EKF = zeros(size(m,1),size(Y,2));
  PP_EKF = zeros(size(m,1),size(m,1),size(Y,2));
  VV_EKF = zeros(size(m,1),size(Y,2));
  EE_EKF = zeros(size(m,1),size(Y,2));
  
  % Check derivatives (should be OK)
  %der_check(func_a, func_da, 1, m0, {dt,b0,H0,Gm0,R0});
  %der_check(func_h, func_dh, 1, m0, {xr,yr});

  clf; clc; 
  disp(['This is a demonstration for tracking a reentry vehicle ',...
        'using 1st order EKF and augmented UKF.'])
  disp(' ');
  
  
%% Extended Kalman filter
  
  fprintf('Running EKF...'); 
  % Filtering with EKF
  for k=1:size(Y,2)
    [m,P] = ekf_predict1(m,P,func_df,Q,func_f,[],{dt,b0,H0,Gm0,R0});
    [m,P] = ekf_update1(m,P,Y(:,k),func_dh,diag([vr va]),func_h,[],{xr,yr});
    MM_EKF(:,k) = m;
    PP_EKF(:,:,k) = P;
    VV_EKF(:,k) = diag(P);
    EE_EKF(:,k) = (X(:,k) - m).^2;
  end

  %
  % Calculate RMSE
  %
  ekf_rmse = sqrt(mean(sum((X(1:2,:)-MM_EKF(1:2,:)).^2)));
  ME_EKF = squeeze(PP_EKF(1,1,:)+PP_EKF(2,2,:));
  fprintf('Done!\n')

  fprintf('Running smoothers...');
  %
  % Smoother 1
  %
  [SM_ERTS,SP_ERTS] = erts_smooth1(MM_EKF,PP_EKF,func_df,Qc*dt,func_f,L,...
                          {dt,b0,H0,Gm0,R0});
  eks_rmse1 = sqrt(mean(sum((X(1:2,:)-SM_ERTS(1:2,:)).^2)));
  ME_ERTS = squeeze(SP_ERTS(1,1,:)+SP_ERTS(2,2,:));


  SV_ERTS = zeros(size(m,1),size(Y,2));
  SE_ERTS = zeros(size(m,1),size(Y,2));
  for k=1:size(Y,2)
    SV_ERTS(:,k) = diag(SP_ERTS(:,:,k));
    SE_ERTS(:,k) = (X(:,k) - SM_ERTS(:,k)).^2;
  end

  %
  % Smoother 2
  %
  [SM_ETF,SP_ETF] = etf_smooth1(MM_EKF,PP_EKF,Y,...
	func_df,Qc*dt,func_if,L,{dt,b0,H0,Gm0,R0},...
	func_dh,diag([vr va]),func_h,[],{xr,yr});
  
  eks_rmse2 = sqrt(mean(sum((X(1:2,:)-SM_ETF(1:2,:)).^2)));
  ME_ETF = squeeze(SP_ETF(1,1,:)+SP_ETF(2,2,:));
  fprintf('Done!\n');
  
  
%% Unscented Kalman filter

  fprintf('Running UKF...'); 

  % Initial values and space for (augmented) UKF 
  m = m0;
  P = P0;
  MM_UKF = zeros(size(m,1),size(Y,2));
  PP_UKF = zeros(size(m,1),size(m,1),size(Y,2));
  VV_UKF = zeros(size(m,1),size(Y,2));
  EE_UKF = zeros(size(m,1),size(Y,2));

  % Filtering with UKF
  for k=1:size(Y,2)
    % Non-augmented UKF
    %[m,P] = ukf_predict1(m,P,func_a,Q,{dt,b0,H0,Gm0,R0});

    % Augmented UKF with separate sigma points
    %[m,P] = ukf_predict2(m,P,func_f,Qc*dt,{dt,b0,H0,Gm0,R0,L});
    
    %[m,P] = ukf_update2(m,P,Y(:,k),func_h,diag([vr va]),{xr,yr});
    
    % Augmented UKF with same sigma points for predict and update steps
    [m,P,X_s,w] = ukf_predict3(m,P,func_f,Qc*dt,diag([vr va]),d_param);  
    [m,P] = ukf_update3(m,P,Y(:,k),func_h,diag([vr va]),X_s,w,h_param,h_param);    

    MM_UKF(:,k) = m;
    PP_UKF(:,:,k) = P;
    VV_UKF(:,k) = diag(P);
    EE_UKF(:,k) = (X(:,k) - m).^2;
  end

  %
  % Calculate RMSE of UKF
  %
  ukf_rmse = sqrt(mean(sum((X(1:2,:)-MM_UKF(1:2,:)).^2)));
  ME_UKF = squeeze(PP_UKF(1,1,:)+PP_UKF(2,2,:));
  fprintf('Done!\n');
  
  fprintf('Running smoothers...');
  % 
  % Smoother 1
  %
  [SM_URTS,SP_URTS] = urts_smooth1(MM_UKF,PP_UKF,func_f,Q,d_param);
  uks_rmse1 = sqrt(mean(sum((X(1:2,:)-SM_URTS(1:2,:)).^2)));
  ME_URTS = squeeze(SP_URTS(1,1,:)+SP_URTS(2,2,:));


  SV_URTS = zeros(size(m,1),size(Y,2));
  SE_URTS = zeros(size(m,1),size(Y,2));
  for k=1:size(Y,2)
    SV_URTS(:,k) = diag(SP_URTS(:,:,k));
    SE_URTS(:,k) = (X(:,k) - SM_URTS(:,k)).^2;
  end
  
  [SM_URTSb,SP_URTSb] = urts_smooth2(MM_UKF,PP_UKF,func_f,Qc*dt,d_param);
  uks_rmse1b = sqrt(mean(sum((X(1:2,:)-SM_URTSb(1:2,:)).^2)));
  ME_URTSb = squeeze(SP_URTSb(1,1,:)+SP_URTSb(2,2,:));

  %
  % Smoother 2
  %
  [SM_UTF,SP_UTF] = utf_smooth1(MM_UKF,PP_UKF,Y,...
	func_if,Qc*dt,d_param,...
	func_h,diag([vr va]),h_param);
  
  uks_rmse2 = sqrt(mean(sum((X(1:2,:)-SM_UTF(1:2,:)).^2)));
  ME_UTF = squeeze(SP_UTF(1,1,:)+SP_UTF(2,2,:));

  fprintf('Done!\n');

  
%% Gauss-Hermite Kalman filter
  
  fprintf('Running GHKF...'); 
  
  % Initial values and space for GHKF
  m = m0;
  P = P0;
  MM_GHKF = zeros(size(m,1),size(Y,2));
  PP_GHKF = zeros(size(m,1),size(m,1),size(Y,2));
  VV_GHKF = zeros(size(m,1),size(Y,2));
  EE_GHKF = zeros(size(m,1),size(Y,2));
  
  % Filtering with GHKF
  for k=1:size(Y,2)
    [m,P] = ghkf_predict(m,P,func_f,Q,{dt,b0,H0,Gm0,R0},3);
    [m,P] = ghkf_update(m,P,Y(:,k),func_h,diag([vr va]),{xr,yr},3);
    MM_GHKF(:,k) = m;
    PP_GHKF(:,:,k) = P;
    VV_GHKF(:,k) = diag(P);
    EE_GHKF(:,k) = (X(:,k) - m).^2;
  end

  %
  % Calculate RMSE
  %
  ghkf_rmse = sqrt(mean(sum((X(1:2,:)-MM_GHKF(1:2,:)).^2)));
  ME_GHKF = squeeze(PP_GHKF(1,1,:)+PP_GHKF(2,2,:));
  fprintf('Done!\n')

  fprintf('Running smoother...');
  %
  % Smoother
  %
  [SM_GHRTS,SP_GHRTS] = ghrts_smooth(MM_GHKF,PP_GHKF,func_f,Q,...
                          {dt,b0,H0,Gm0,R0},3);
  ghrts_rmse1 = sqrt(mean(sum((X(1:2,:)-SM_GHRTS(1:2,:)).^2)));
  ME_GHRTS = squeeze(SP_GHRTS(1,1,:)+SP_GHRTS(2,2,:));


  SV_GHRTS = zeros(size(m,1),size(Y,2));
  SE_GHRTS = zeros(size(m,1),size(Y,2));
  for k=1:size(Y,2)
    SV_GHRTS(:,k) = diag(SP_GHRTS(:,:,k));
    SE_GHRTS(:,k) = (X(:,k) - SM_GHRTS(:,k)).^2;
  end

  fprintf('Done!\n');

  
%% Cubature Kalman filter
  
  fprintf('Running CKF...'); 
  
  % Initial values and space for CKF
  m = m0;
  P = P0;
  MM_CKF = zeros(size(m,1),size(Y,2));
  PP_CKF = zeros(size(m,1),size(m,1),size(Y,2));
  VV_CKF = zeros(size(m,1),size(Y,2));
  EE_CKF = zeros(size(m,1),size(Y,2));
  
  % Filtering with CKF
  for k=1:size(Y,2)
    [m,P] = ckf_predict(m,P,func_f,Q,{dt,b0,H0,Gm0,R0});
    [m,P] = ckf_update(m,P,Y(:,k),func_h,diag([vr va]),{xr,yr});
    MM_CKF(:,k) = m;
    PP_CKF(:,:,k) = P;
    VV_CKF(:,k) = diag(P);
    EE_CKF(:,k) = (X(:,k) - m).^2;
  end

  %
  % Calculate RMSE
  %
  ckf_rmse = sqrt(mean(sum((X(1:2,:)-MM_CKF(1:2,:)).^2)));
  ME_CKF = squeeze(PP_CKF(1,1,:)+PP_CKF(2,2,:));
  fprintf('Done!\n')

  fprintf('Running smoother...');
  %
  % Smoother
  %
  [SM_CRTS,SP_CRTS] = crts_smooth(MM_CKF,PP_CKF,func_f,Q,...
                          {dt,b0,H0,Gm0,R0});
  crts_rmse1 = sqrt(mean(sum((X(1:2,:)-SM_CRTS(1:2,:)).^2)));
  ME_CRTS = squeeze(SP_CRTS(1,1,:)+SP_CRTS(2,2,:));


  SV_CRTS = zeros(size(m,1),size(Y,2));
  SE_CRTS = zeros(size(m,1),size(Y,2));
  for k=1:size(Y,2)
    SV_CRTS(:,k) = diag(SP_CRTS(:,:,k));
    SE_CRTS(:,k) = (X(:,k) - SM_CRTS(:,k)).^2;
  end

  fprintf('Done!\n');
  
%% Visualize methods
  
  if ~silent
    aa = 0.02*(-1:0.1:4);
    cx = R0 * cos(aa);
    cy = R0 * sin(aa);
    % Plot EKF estimate
    plot(xr,yr,'ko',cx,cy,'r-',X(1,:),X(2,:),'g-',...
         MM_EKF(1,:),MM_EKF(2,:),'k--');
    legend('Radar','Earth','True','Estimate');
    title('Filtering result with EKF');
    disp(' ');
    disp('Filtering result with EKF is now displayed.');
    disp(' ');
    disp('<press any key to see the estimation error of x_1>');    
    pause;
  
    % Error for x_1 with EKF
    semilogy(T,EE_EKF(1,:),'g-',T,VV_EKF(1,:),'b--',...
	     T,SE_ERTS(1,:),'r-',T,SV_ERTS(1,:),'k--');
    legend('EKF-RMSE','EKF-STDE',...
	   'ERTS-RMSE1','ERTS-STDE1');  
    title('RMSE of estimating x_1 with EKF and ERTS')
    clc;
    disp('RMSE of estimating x_1 with EKF and ERTS is now displayed.');
    disp(' ');
    disp('<press any key to see the estimation error of x_5>');
    pause  
    
    % Error for x_5 with EKF
    semilogy(T,EE_EKF(5,:),'g-',T,VV_EKF(5,:),'b--',...
	     T,SE_ERTS(5,:),'r-',T,SV_ERTS(5,:),'k--');
    legend('EKF-RMSE','EKF-STDE',...
	   'ERTS-RMSE1','ERTS-STDE1');
    title('RMSE of estimating x_5 with EKF and ERTS')
    clc;
    disp('RMSE of estimating x_5 with EKF and ERTS is now displayed.');
    disp(' ');
    disp('<press any key to see the estimation error of x_5>');
    pause  
    
    % Plot UKF estimate
    plot(xr,yr,'ko',cx,cy,'r-',X(1,:),X(2,:),'g-',...
         MM_UKF(1,:),MM_UKF(2,:),'k--');
    legend('Radar','Earth','True','Estimate');
    clc;
    disp('Filtering result with UKF is now displayed.');
    disp(' ');
    disp('<press any key to see the estimation error of x_1>');    
    pause;
  
    % Error for x_1 with UKF
    semilogy(T,EE_UKF(1,:),'g-',T,VV_UKF(1,:),'b--',...
             T,SE_URTS(1,:),'r-',T,SV_URTS(1,:),'k--');
    legend('UKF-RMSE','UKF-STDE','URTS-RMSE','URTS-STDE');
    title('RMSE of estimating x_1 with UKF and URTS')    
    clc;
    disp('RMSE of estimating x_1 with EKF and ERTS is now displayed.');
    disp(' ');
    disp('<press any key to see the estimation error of x_5>');
    pause;
  
    % Error for x_5 with UKF
    semilogy(T,EE_UKF(5,:),'g-',T,VV_UKF(5,:),'b--',...
             T,SE_URTS(5,:),'r-',T,SV_URTS(5,:),'k--');
    legend('UKF-RMSE','UKF-STDE','URTS-RMSE','URTS-STDE');
    title('RMSE of estimating x_5 with UKF and URTS')
    clc;
    disp('RMSE of estimating x_5 with EKF and ERTS is now displayed.');
    
    
    % Plot GHKF estimate
    plot(xr,yr,'ko',cx,cy,'r-',X(1,:),X(2,:),'g-',...
         MM_GHKF(1,:),MM_GHKF(2,:),'k--');
    legend('Radar','Earth','True','Estimate');
    clc;
    disp('Filtering result with GHKF is now displayed.');
    disp(' ');
    disp('<press any key to see the estimation error of x_1>');    
    pause;
  
    % Error for x_1 with GHKF
    semilogy(T,EE_UKF(1,:),'g-',T,VV_UKF(1,:),'b--',...
             T,SE_URTS(1,:),'r-',T,SV_URTS(1,:),'k--');
    legend('GHKF-RMSE','GHKF-STDE','GHRTS-RMSE','GHRTS-STDE');
    title('RMSE of estimating x_1 with GHKF and GHRTS')    
    clc;
    disp('RMSE of estimating x_1 with GHKF and GHRTS is now displayed.');
    disp(' ');
    disp('<press any key to see the estimation error of x_5>');
    pause;
  
    % Error for x_5 with GHKF
    semilogy(T,EE_GHKF(5,:),'g-',T,VV_GHKF(5,:),'b--',...
             T,SE_GHRTS(5,:),'r-',T,SV_GHRTS(5,:),'k--');
    legend('GHKF-RMSE','GHKF-STDE','GHRTS-RMSE','GHRTS-STDE');
    title('RMSE of estimating x_5 with GHKF and GHRTS')
    clc;
    disp('RMSE of estimating x_5 with EKF and GHRTS is now displayed.');
    
    
    % Plot CKF estimate
    plot(xr,yr,'ko',cx,cy,'r-',X(1,:),X(2,:),'g-',...
         MM_CKF(1,:),MM_CKF(2,:),'k--');
    legend('Radar','Earth','True','Estimate');
    clc;
    disp('Filtering result with CKF is now displayed.');
    disp(' ');
    disp('<press any key to see the estimation error of x_1>');    
    pause;
  
    % Error for x_1 with CKF
    semilogy(T,EE_UKF(1,:),'g-',T,VV_UKF(1,:),'b--',...
             T,SE_URTS(1,:),'r-',T,SV_URTS(1,:),'k--');
    legend('CKF-RMSE','CKF-STDE','CRTS-RMSE','CRTS-STDE');
    title('RMSE of estimating x_1 with CKF and CRTS')    
    clc;
    disp('RMSE of estimating x_1 with CKF and CRTS is now displayed.');
    disp(' ');
    disp('<press any key to see the estimation error of x_5>');
    pause;
  
    % Error for x_5 with CKF
    semilogy(T,EE_CKF(5,:),'g-',T,VV_CKF(5,:),'b--',...
             T,SE_CRTS(5,:),'r-',T,SV_CRTS(5,:),'k--');
    legend('CKF-RMSE','CKF-STDE','CRTS-RMSE','CRTS-STDE');
    title('RMSE of estimating x_5 with CKF and CRTS')
    clc;
    disp('RMSE of estimating x_5 with EKF and CRTS is now displayed.');

    
  end

  
%% Show RMSE errors for each method
  
  disp(' ');
  disp('RMS errors:');
  fprintf('EKF-RMSE   = %.6f [%.6f]\n',ekf_rmse,sqrt(mean(ME_EKF)));
  fprintf('ERTS-RMSE  = %.6f [%.6f]\n',eks_rmse1,sqrt(mean(ME_ERTS)));
  fprintf('ETF-RMSE   = %.6f [%.6f]\n',eks_rmse2,sqrt(mean(ME_ETF)));
  fprintf('UKF-RMSE   = %.6f [%.6f]\n',ukf_rmse,sqrt(mean(ME_UKF)));
  fprintf('URTS1-RMSE = %.6f [%.6f]\n',uks_rmse1,sqrt(mean(ME_URTS)));
  fprintf('URTS2-RMSE = %.6f [%.6f]\n',uks_rmse1b,sqrt(mean(ME_URTSb)));
  fprintf('UTF-RMSE   = %.6f [%.6f]\n',uks_rmse2,sqrt(mean(ME_UTF)));

  % Cubature methods
  fprintf('GHKF-RMSE  = %.6f [%.6f]\n',ghkf_rmse,sqrt(mean(ME_GHKF)));
  fprintf('GHRTS-RMSE = %.6f [%.6f]\n',ghrts_rmse1,sqrt(mean(ME_GHRTS)));
  fprintf('CKF-RMSE   = %.6f [%.6f]\n',ckf_rmse,sqrt(mean(ME_CKF)));
  fprintf('CRTS-RMSE  = %.6f [%.6f]\n',crts_rmse1,sqrt(mean(ME_CRTS)));
  
  
  