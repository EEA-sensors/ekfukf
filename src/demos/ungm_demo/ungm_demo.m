% Demonstration of univariate nonstationary growth model (UNGM) using EKF and UKF.
%
% Copyright (C) 2007 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

silent = 0;

% Derivatives
f_func = @ungm_f;
h_func = @ungm_h;
df_func = @ungm_df_dx;
dh_func = @ungm_dh_dx;
d2f_func = @ungm_d2f_dx2;
d2h_func = @ungm_d2h_dx2;

n = 500;
x_0 = .1;
P_0 = 1;

Y = zeros(1,n);

u_n = 1;
v_n = 1;

% Generate the true states with process noise
X = zeros(1,n);
X(1) = ungm_f(x_0,1) + gauss_rnd(0,u_n,1);
for i = 2:n
    X(i) = feval(f_func,X(i-1),i) + gauss_rnd(0,u_n,1);
end
    
% Generate the observations with measurement noise
for i = 1:n
    Y(i) = feval(h_func,X(i)) + gauss_rnd(0,v_n,1);
end

Y_r = feval(h_func,X);

% Parameters for dynamic model. Used by smoothers.
params = cell(size(Y,2));
for i = 1:size(Y,2)
   params{i} = i+1; 
end

% Filtering with non-augmented UKF
M = x_0;
P = P_0;

MM_UKF1 = zeros(size(M,1),size(Y,2));
PP_UKF1 = zeros(size(M,1),size(M,1),size(Y,2));

% Filtering loop for UKF1
for k = 1:size(Y,2)
   [M,P] = ukf_predict1(M,P,f_func,u_n,k);
   [M,P] = ukf_update1(M,P,Y(:,k),h_func,v_n);
   MM_UKF1(:,k)   = M;
   PP_UKF1(:,:,k) = P;    
end
[MMS_URTS1, PPS_URTS1] = urts_smooth1(MM_UKF1,PP_UKF1,f_func,u_n,params,[],[],[],[],0);

% MSE for UKF1
UKF1_MSE = sum((X-MM_UKF1).^2)/n;
fprintf('UKF1-MSE = %.4f\n',UKF1_MSE);
UKS1_MSE = sum((X-MMS_URTS1).^2)/n;
fprintf('UKS1-MSE = %.4f\n',UKS1_MSE);


% Filtering with augmented UKF
M = x_0;
P = P_0;

MM_UKF2 = zeros(size(M,1),size(Y,2));
PP_UKF2 = zeros(size(M,1),size(M,1),size(Y,2));

% Filtering loop for UKF2
for k = 1:size(Y,2)
    % Augmented UKF
   [M,P,X_s,w] = ukf_predict3(M,P,f_func,u_n,v_n,k);
   [M,P] = ukf_update3(M,P,Y(:,k),h_func,v_n,X_s,w,[]);
   MM_UKF2(:,k)   = M;
   PP_UKF2(:,:,k) = P;    
end

[MMS_URTS2, PPS_URTS2] = urts_smooth1(MM_UKF2,PP_UKF2,f_func,u_n,params,[],[],[],[],0);

% MSE for UKF2
UKF2_MSE = sum((X-MM_UKF2).^2)/n;
fprintf('UKF2-MSE = %.4f\n',UKF2_MSE);
UKS2_MSE = sum((X-MMS_URTS2).^2)/n;
fprintf('UKS2-MSE = %.4f\n',UKS2_MSE);

% Filtering with EKF
M = x_0;
P = P_0;

MM_EKF = zeros(size(M,1),size(Y,2));
PP_EKF = zeros(size(M,1),size(M,1),size(Y,2));

% Filtering loop for EKF
for k = 1:size(Y,2)
   [M,P] = ekf_predict1(M,P,df_func,u_n,f_func,[],k);
   [M,P] = ekf_update1(M,P,Y(:,k),dh_func,v_n,h_func,[],[]);
   MM_EKF(:,k)   = M;
   PP_EKF(:,:,k) = P;    
end

params = cell(size(Y,2));
for i = 1:size(Y,2)
   params{i} = i+1; 
end

[MMS_ERTS,PPS_ERTS] = erts_smooth1(MM_EKF,PP_EKF,df_func,u_n,f_func,[],params,0);


EKF_MSE = sum((X-MM_EKF).^2)/n;
fprintf('EKF-MSE = %.4f\n',EKF_MSE);
ERTS_MSE = sum((X-MMS_ERTS).^2)/n;
fprintf('ERTS-MSE = %.4f\n',ERTS_MSE);

% Filtering with bootstrap filter
M = x_0;
P = P_0;
n_particles = 1000;
SX = gauss_rnd(M,P,n_particles);

MM_BS = zeros(size(M,1),size(Y,2));
PP_BS = zeros(size(M,1),size(M,1),size(Y,2));

% Filtering loop for bootstrap filter
for k = 1:size(Y,2)
   SX = ungm_f(SX,k) + gauss_rnd(0,1,size(SX,2));
   W  = gauss_pdf(Y(:,k),ungm_h(SX,k),1);
   ind = resampstr(W);
   SX = SX(:,ind);
   M = mean(SX);
   P = var(SX);
   MM_BS(:,k)   = M;
   PP_BS(:,:,k) = P;    
end

BS_MSE = sum((X-MM_BS).^2)/n;
fprintf('BS-MSE = %.4f\n',BS_MSE);

if ~silent
  subplot(3,1,1);
  plot(1:100,X(1:100),'-kx',1:100,MM_UKF2(1:100),'--bo')
  title('UKF2 filtering result');
  xlim([0 100]);
  ylim([-20 20]);
  legend('Real signal', 'UKF2 filtered estimate');

  subplot(3,1,2);
  plot(1:100,X(1:100),'-kx',1:100,MM_EKF(1:100),'--bo')
  title('EKF filtering result');
  xlim([0 100]);
  ylim([-20 20]);
  legend('Real signal', 'EKF filtered estimate');

  subplot(3,1,3);
  plot(1:100,X(1:100),'-kx',1:100,MM_BS(1:100),'--bo')
  title('Bootstrap filtering result')
  xlim([0 100]);
  ylim([-20 20]);
  legend('Real signal','Filtered estimates')

  %print -dpsc ungm_states.ps

  pause;
  
  subplot(2,1,1);
  plot(1:100,X(1:100),'-kx',1:100,MM_UKF1(1:100),'--bo')
  title('UKF1 filtering result');
  xlim([0 100]);
  ylim([-20 20]);
  legend('Real signal', 'UKF1 filtered estimate');
  
  subplot(2,1,2);
  plot(1:100,X(1:100),'-kx',1:100,MM_UKF2(1:100),'--bo')
  title('UKF2 filtering result');
  xlim([0 100]);
  ylim([-20 20]);
  legend('Real signal', 'UKF2 filtered estimate');
  
  %print -dpsc ungm_ukf_comp.ps
  
  pause

  subplot(3,1,1);
  E = X-MM_UKF2;
  PP_UKF2 = squeeze(PP_UKF2);
  plot(1:100,E(1:100),'-kx',1:100,3*sqrt(PP_UKF2(1:100))','--r',1:100,-3*sqrt(PP_UKF2(1:100))','--r');
  title('UKF2 error');
  legend('Estimation error of UKF2','3\sigma interval')

  subplot(3,1,2);
  E = X-MM_EKF;
  PP_EKF = squeeze(PP_EKF);
  plot(1:100,E(1:100),'-kx',1:100,3*sqrt(PP_EKF(1:100))','--r',1:100,-3*sqrt(PP_EKF(1:100))','--r');
  title('EKF error');
  legend('Estimation error of EKF','3\sigma interval');

  subplot(3,1,3);
  E = X-MM_BS;
  PP_BS = squeeze(PP_BS);
  plot(1:100,E(1:100),'-kx',1:100,3*sqrt(PP_BS(1:100))','--r',1:100,-3*sqrt(PP_BS(1:100))','--r');
  title('Bootstrap filtering error');
  legend('Estimation error with BS','3\sigma interval')

  %print -dpsc ungm_c_errors.ps
  
  pause

  subplot(2,1,1);
  E = X-MM_UKF1;
  PP_UKF1 = squeeze(PP_UKF1);
  plot(1:100,E(1:100),'-kx',1:100,3*sqrt(PP_UKF1(1:100))','--r',1:100,-3*sqrt(PP_UKF1(1:100))','--r');
  title('UKF1 error');
  legend('Estimation error of UKF1','3\sigma interval')

  subplot(2,1,2);
  E = X-MM_UKF2;
  PP_UKF2 = squeeze(PP_UKF2);
  plot(1:100,E(1:100),'-kx',1:100,3*sqrt(PP_UKF2(1:100))','--r',1:100,-3*sqrt(PP_UKF2(1:100))','--r');
  title('UKF2 error');
  legend('Estimation error of UKF2','3\sigma interval')
  
  %print -dpsc ungm_ukf_comp_error.ps
end