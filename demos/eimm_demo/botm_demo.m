%EIMM_DEMO2 Bearings Only Tracking of a Manouvering Target demonstration
% 
% Simple demonstration for non-linear IMM using the following models:
%  1. Standard Wiener process velocity model
%  2. Coordinated turn model
%
% The measurement model is non-linear bearings only model.
% See BOT-demo or the documentation for more details.
% 
% Copyright (C) 2007-2008 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

print_figures = 1;
save_figures = 1;

% Dimensionality of the state space
dims = 5;

% The number of models in use
nmodels = 2;

% Step size
dt = 0.1;

% Space for function handles and parameters
a_func = {};
ia_func = {};
a_param = {};
h_func = {};
h_param = {};
dh_dx_func = {};

a_func{1} = [];
a_func{2} = @f_turn;
ia_func{1} = [];
ia_func{2} = @f_turn_inv;
a_param{1} = [];
a_param{2} = {dt};


% Space for model parameters
ind = cell(1,nmodels);
F   = cell(1,nmodels);
L   = cell(1,nmodels);
Qc  = cell(1,nmodels);
A   = cell(1,nmodels);
Q   = cell(1,nmodels);
H   = cell(1,nmodels);
R   = cell(1,nmodels);

% Index vector of model 1
ind{1} = [1 2 3 4]';

% Transition matrix for the continous-time velocity model.
F{1} = [0 0 1 0;
        0 0 0 1;
        0 0 0 0;
        0 0 0 0];

% Noise effect matrix for the continous-time system.
L{1} = [0 0;
        0 0;
        1 0;
        0 1];

% Process noise variance
q1 = 0.01;
Qc{1} = diag([q1 q1]);

% Discretization of the continous-time system.
[A{1},Q{1}] = lti_disc(F{1},L{1},Qc{1},dt);

% Process noise variance for model 1 using EKF
EKF_q1 = .05;
EKF_Qc1 = diag([EKF_q1 EKF_q1]);

% Discretization of the continous-time system.
[EKF_A1,EKF_Q1] = lti_disc(F{1},L{1},EKF_Qc1,dt);


%%% Specification of the turning model

% System components. 5th parameter is the turning rate 
ind{2} = [1 2 3 4 5]';

% Derivative of the dynamic function 
A{2} = @f_turn_dx;
% Process noise for the turning rate
Qc{2} = 0.15;

% Noise effect matrix
L{2} = [0 0 0 0 1]';

% Process noise covariance
Q{2} = L{2}*Qc{2}*L{2}'*dt;

hdims = 2;

%mu_ip = [0.90 0.05 0.05];
mu_ip = [0.95 0.05];
mu_0j = mu_ip;
%p_ij = [0.65 0.35;
%        0.10 0.90];

p_ij = [0.90 0.10;
        0.10 0.90];

% Number of data points
n = 200;

% Space for real states and modes
X_r = zeros(dims,n);
mstate = zeros(1,n);

%%%%%%% Creation of trajectory %%%%%%%

% Start with constant velocity 1 toward right
mstate(1:40) = 1;
X_r(:,1) = [0 0 1 0 0]';

% At 4s make a turn left with rate 1 
mstate(41:90) = 2;
X_r(5,40) = 1;

% At 9s move straight for 2 seconds
mstate(91:110) = 1;

% At 11s commence another turn right with rate -1
mstate(111:160) = 2;
X_r(5,110) = -1;

% At 16s move straight for 4 seconds
mstate(161:200) = 1;

% Generate object state values
for i = 2:n
   st = mstate(i);
   if isstr(a_func{st}) | strcmp(class(a_func{st}),'function_handle')
       X_r(ind{st},i) = feval(a_func{st},X_r(ind{st},i-1),a_param{st});
   else 
       X_r(ind{st},i) = A{st}*X_r(ind{st},i-1);
   end
end

% Positions of sensors
% $$$ S1 = [-0.5; 3];
% $$$ S2 = [-0.5;-3];
% $$$ S3 = [   7;-3];
% $$$ S4 = [   7; 3];

S1 = [-0.5; 3.5];
S2 = [-0.5;-3.5];
S3 = [   7;-3.5];
S4 = [   7; 3.5];

s = [S1 S2 S3 S4];

% Handles to measurement models
H{1} = @bot_dh_dx;
H{2} = @bot_dh_dx;

h_func{1} = @bot_h;
h_func{2} = @bot_h;

h_param{1} = s;
h_param{2} = s;

% Standard deviation
sd = 0.1*ones(1,size(s,2));
R = {};
R{1} = diag(sd.^2);
R{2} = R{1};
% Generate measurements
Y = bot_h(X_r,s);
% Add noise
for i = 1:size(Y,2)
    st = mstate(i);
    for j = 1:size(Y,1)
        Y(j,i) = Y(j,i) + sqrt(R{st}(j,j)) * randn; 
    end
end

% Print the trajectory
if print_figures
    h = plot(X_r(1,:),X_r(2,:),'-g',...
             s(1,:),s(2,:),'k^',... 
             X_r(1,1),X_r(2,1),'ro','MarkerSize',12);    
    legend('Real trajectory',...
           'Positions of sensors',...
           'Starting position', 'Location', 'North');
    set(h,'markersize',5);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',8);
    xlim([-1 7.5]) 
    ylim([-4 4]) 
    if save_figures
        print('-depsc','eimm2_trajectory.eps');
    end
    pause
end


m = [0 0 0 -1 0]';
P = diag([10.1 10.1 1.1 1.1 1]);

%% Space for the estimates.

% EKF with model 1
EKF_MM = zeros(size(A{1},1), size(Y,2));
EKF_PP = zeros(size(A{1},1), size(A{1},1), size(Y,2));

% UKF with model 1
UKF_MM = zeros(size(A{1},1), size(Y,2));
UKF_PP = zeros(size(A{1},1), size(A{1},1), size(Y,2));

% EKF based IMM
EIMM_MM = zeros(size(m,1), size(Y,2));
EIMM_PP = zeros(size(m,1), size(m,1), size(Y,2));
EIMM_MM_i = cell(nmodels,n);
EIMM_PP_i = cell(nmodels,n);
EIMM_MU = zeros(nmodels,size(Y,2));

% UKF based IMM
UIMM_MM = zeros(size(m,1), size(Y,2));
UIMM_PP = zeros(size(m,1), size(m,1), size(Y,2));
UIMM_MM_i = cell(nmodels,n);
UIMM_PP_i = cell(nmodels,n);
UIMM_MU = zeros(nmodels,size(Y,2));

%%% Initial estimates %%%

% EKF with model 1
EKF_M = [0 0 0 -1]';
EKF_P = diag([1.1 1.1 0.1 0.1]);

% UKF with model 1
UKF_M = [0 0 0 -1]';
UKF_P = diag([1.1 1.1 0.1 0.1]);

% EKF based IMM
x_ip1{1} = [0 0 1 0]';
x_ip1{2} = [0 0 1 0 0]';
mu_ip1 = mu_ip;

P_ip1{1} = diag([0.1 0.1 0.1 0.1]);
P_ip1{2} = diag([0.1 0.1 0.1 0.1 1]);

% UKF based IMM
x_ip2{1} = [0 0 1 0]';
x_ip2{2} = [0 0 1 0 0]';
mu_ip2 = mu_ip;

P_ip2{1} = diag([0.1 0.1 0.1 0.1]);
P_ip2{2} = diag([0.1 0.1 0.1 0.1 1]);


% Filtering steps.
for i = 1:size(Y,2)
    % EKF with model 1
    [EKF_M,EKF_P] = kf_predict(EKF_M,EKF_P,EKF_A1,EKF_Q1);
    [EKF_M,EKF_P] = ekf_update1(EKF_M,EKF_P,Y(:,i),H{1},R{1},h_func{1},[],h_param{1});
    
    EKF_MM(:,i)   = EKF_M;
    EKF_PP(:,:,i) = EKF_P;

    % UKF with model 2
    [UKF_M,UKF_P] = kf_predict(UKF_M,UKF_P,EKF_A1,EKF_Q1);
    [UKF_M,UKF_P] = ukf_update1(UKF_M,UKF_P,Y(:,i),h_func{1},R{1},h_param{1});
    
    UKF_MM(:,i)   = UKF_M;
    UKF_PP(:,:,i) = UKF_P;
    
    % EKF based IMM
    [x_p1,P_p1,c_j1] = eimm_predict(x_ip1,P_ip1,mu_ip1,p_ij,ind,dims,A,a_func,a_param,Q);
    [x_ip1,P_ip1,mu_ip1,m1,P1] = eimm_update(x_p1,P_p1,c_j1,ind,dims,Y(:,i),H,h_func,R,h_param);
    EIMM_MM(:,i)   = m1;
    EIMM_PP(:,:,i) = P1;
    EIMM_MU(:,i)   = mu_ip1';
    EIMM_MM_i(:,i) = x_ip1';
    EIMM_PP_i(:,i) = P_ip1';

    % UKF based IMM    
    [x_p2,P_p2,c_j2] = uimm_predict(x_ip2,P_ip2,mu_ip2,p_ij,ind,dims,A,a_func,a_param,Q);
    [x_ip2,P_ip2,mu_ip2,m2,P2] = uimm_update(x_p2,P_p2,c_j2,ind,dims,Y(:,i),H,h_func,R,h_param);
    
    UIMM_MM(:,i)   = m2;
    UIMM_PP(:,:,i) = P2;
    UIMM_MU(:,i)   = mu_ip2';
    UIMM_MM_i(:,i) = x_ip2';
    UIMM_PP_i(:,i) = P_ip2';
    
    % Plot the estimates so far
    if print_figures
        plot(EKF_MM(1,1:i),EKF_MM(2,1:i),'y-',...
             EIMM_MM(1,1:i),EIMM_MM(2,1:i),'r-',...
             UIMM_MM(1,1:i),UIMM_MM(2,1:i),'b-',...
             X_r(1,1:i),X_r(2,1:i),'g-');
        
        % Measurement directions
        hold on
        for k = 1:size(s,2)
            len = sqrt(sum((X_r(1:2,i)-s(:,k)).^2,1));
            dx = len*cos(Y(k,i));
            dy = len*sin(Y(k,i));
            
            plot([s(1,k);s(1,k)+dx], [s(2,k);s(2,k)+dy], 'k--')
        end
        plot(s(1,:),s(2,:),'k^')
        xlim([-1 7.5]) 
        ylim([-4 4]) 
        hold off
        drawnow
    end
end

% Smooth with EKF based IMM smoother
[SMI_1,SPI_1,SMI_i_1,SPI_i_1,MU_S_1] = eimm_smooth(EIMM_MM,EIMM_PP,EIMM_MM_i,EIMM_PP_i,EIMM_MU,p_ij,mu_0j,ind,dims,A,ia_func,a_param,Q,R,H,h_func,h_param,Y);

% Smooth with UKF based IMM smoother
[SMI_2,SPI_2,SMI_i_2,SPI_i_2,MU_S_2] = uimm_smooth(UIMM_MM,UIMM_PP,UIMM_MM_i,UIMM_PP_i,UIMM_MU,p_ij,mu_0j,ind,dims,A,ia_func,a_param,Q,R,H,h_func,h_param,Y);

% Smooth the EKF estimates with RTS smoother using model 1
[SM1,SP1] = rts_smooth(EKF_MM,EKF_PP,EKF_A1,EKF_Q1);

% Smooth the UKF estimates with RTS smoother using model 1
[SM2,SP2] = rts_smooth(UKF_MM,UKF_PP,EKF_A1,EKF_Q1);

% Calculate the MSEs
MSE_EKF1_1 = mean((X_r(1,:)-EKF_MM(1,:)).^2);
MSE_EKF1_2 = mean((X_r(2,:)-EKF_MM(2,:)).^2);
MSE_EKF1 = 1/2*(MSE_EKF1_1 + MSE_EKF1_2);

MSE_EKS1_1 = mean((X_r(1,:)-SM1(1,:)).^2);
MSE_EKS1_2 = mean((X_r(2,:)-SM1(2,:)).^2);
MSE_EKS1 = 1/2*(MSE_EKS1_1 + MSE_EKS1_2);

% Calculate the MSEs
MSE_UKF1_1 = mean((X_r(1,:)-UKF_MM(1,:)).^2);
MSE_UKF1_2 = mean((X_r(2,:)-UKF_MM(2,:)).^2);
MSE_UKF1 = 1/2*(MSE_UKF1_1 + MSE_UKF1_2);

MSE_UKS1_1 = mean((X_r(1,:)-SM2(1,:)).^2);
MSE_UKS1_2 = mean((X_r(2,:)-SM2(2,:)).^2);
MSE_UKS1 = 1/2*(MSE_UKS1_1 + MSE_UKS1_2);

MSE_EIMM1 = mean((X_r(1,:)-EIMM_MM(1,:)).^2);
MSE_EIMM2 = mean((X_r(2,:)-EIMM_MM(2,:)).^2);
MSE_EIMM = 1/2*(MSE_EIMM1 + MSE_EIMM2);

MSE_EIMMS1 = mean((X_r(1,:)-SMI_1(1,:)).^2);
MSE_EIMMS2 = mean((X_r(2,:)-SMI_1(2,:)).^2);
MSE_EIMMS = 1/2*(MSE_EIMMS1 + MSE_EIMMS2);

MSE_UIMM1 = mean((X_r(1,:)-UIMM_MM(1,:)).^2);
MSE_UIMM2 = mean((X_r(2,:)-UIMM_MM(2,:)).^2);
MSE_UIMM = 1/2*(MSE_UIMM1 + MSE_UIMM2);

MSE_UIMMS1 = mean((X_r(1,:)-SMI_2(1,:)).^2);
MSE_UIMMS2 = mean((X_r(2,:)-SMI_2(2,:)).^2);
MSE_UIMMS = 1/2*(MSE_UIMMS1 + MSE_UIMMS2);


fprintf('Mean square errors of position estimates:\n');
fprintf('EKF1-RMSE = %.4f\n',MSE_EKF1);
fprintf('EKS1-RMSE = %.4f\n',MSE_EKS1);
fprintf('UKF1-RMSE = %.4f\n',MSE_UKF1);
fprintf('UKS1-RMSE = %.4f\n',MSE_UKS1);
fprintf('EIMM-RMSE = %.4f\n',MSE_EIMM);
fprintf('EIMMS-RMSE = %.4f\n',MSE_EIMMS);
fprintf('UIMM-RMSE = %.4f\n',MSE_UIMM);
fprintf('UIMMS-RMSE = %.4f\n',MSE_UIMMS);


% Plot the final filtering and smoothing results
if print_figures
    h = plot(X_r(1,:),X_r(2,:),'g-',...         
             EKF_MM(1,:),EKF_MM(2,:),'-k',...
             SM1(1,:),SM1(2,:),'--k',...             
             EIMM_MM(1,:),EIMM_MM(2,:),'-r',...             
             SMI_1(1,:),SMI_1(2,:),'--r',...
             UIMM_MM(1,:),UIMM_MM(2,:),'-b',...
             SMI_2(1,:),SMI_2(2,:),'--b');
    legend('True trajectory',...
           'EKF',...
           'ERTS',...
           'IMM-EKF',...
           'IMM-EKS',...
           'IMM-UKF',...           
           'IMM-UKS');
    %title('Estimates produced by IMM-filter.')
    set(h,'markersize',2);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',8);
    xlim([-1 7.5]) 
    ylim([-4 4]) 
    if save_figures
        print('-depsc','eimm2_1.eps');
    end
    pause

    
    % Determine the real model probabilities
    p_models = zeros(nmodels,n);
    I1 = find(mstate == 1);
    p_models(1,I1) = 1;
    I2 = find(mstate == 2);
    p_models(2,I2) = 1;
 
    % Plot model 2 probability for filters
    h = plot(1:n, p_models(2,:),'g--',...
             1:n,EIMM_MU(2,:)','-r',...
             1:n,MU_S_1(2,:)','--r',...
             1:n,UIMM_MU(2,:)','-b',...
             1:n,MU_S_2(2,:)','--b');
    %1:n,MU_S_2(2,:)','b-');
    legend('True',...
           'IMM-EKF',...
           'IMM-EKS',...
           'IMM-UKF',...
           'IMM-UKS');
    %title('Probability of model 2');
    ylim([-0.1,1.1]);
    set(h,'markersize',2);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',8);
    if save_figures
        print('-depsc','eimm2_2.eps');
    end
    pause
    
    % Collect the turn rate estimates
    EMM_t = zeros(5,n);
    EMM_s = zeros(5,n);
    UMM_t = zeros(5,n);
    UMM_s = zeros(5,n); 
    for i = 1:n
        EMM_t(:,i) = EIMM_MM_i{2,i};
        EMM_s(:,i) = SMI_i_1{2,i};
        UMM_t(:,i) = UIMM_MM_i{2,i};
        UMM_s(:,i) = SMI_i_2{2,i};
    end
    
    % Plot the filtered turn rates
    h = plot(1:n, X_r(5,:),'-g',...
             1:n,EMM_t(5,:),'-r',...     
             1:n,EMM_s(5,:),'--r',...     
             1:n,UMM_t(5,:),'-b',...             
             1:n,UMM_s(5,:),'--b');
    %title('Turn rate estimates')
    legend('True',...
           'IMM-EKF',...
           'IMM-EKS',...
           'IMM-UKF',...
           'IMM-UKS');
    set(h,'markersize',2);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',8);
    if save_figures
        print('-depsc','eimm2_3.eps');
    end
  
end

