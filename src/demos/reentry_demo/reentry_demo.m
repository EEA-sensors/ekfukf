% Discrete-time reentry dynamics demonstration with EKF and UKF. 

% Copyright (C) 2006 Simo Särkkä
%               2007 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.
  function reentry_demo

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
  

  % Initial values and space for (augmented) UKF 
  m = m0;
  P = P0;
  MM_UKF = zeros(size(m,1),size(Y,2));
  PP_UKF = zeros(size(m,1),size(m,1),size(Y,2));
  VV_UKF = zeros(size(m,1),size(Y,2));
  EE_UKF = zeros(size(m,1),size(Y,2));
  
  fprintf('Running UKF...'); 

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
  
  %
  % Plot
  %
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
  end

  disp(' ');
  disp('RMS errors:');
  fprintf('EKF-RMSE = %.6f [%.6f]\n',ekf_rmse,sqrt(mean(ME_EKF)));
  fprintf('ERTS-RMSE = %.6f [%.6f]\n',eks_rmse1,sqrt(mean(ME_ERTS)));
  fprintf('ETF-RMSE = %.6f [%.6f]\n',eks_rmse2,sqrt(mean(ME_ETF)));
  fprintf('UKF-RMSE = %.6f [%.6f]\n',ukf_rmse,sqrt(mean(ME_UKF)));
  fprintf('URTS1-RMSE = %.6f [%.6f]\n',uks_rmse1,sqrt(mean(ME_URTS)));
  fprintf('URTS2-RMSE = %.6f [%.6f]\n',uks_rmse1b,sqrt(mean(ME_URTSb)));
  fprintf('UTF-RMSE = %.6f [%.6f]\n',uks_rmse2,sqrt(mean(ME_UTF)));

  