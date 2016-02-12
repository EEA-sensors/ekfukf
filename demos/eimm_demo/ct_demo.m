%EIMM_DEMO1 Coordinated turn model demonstration 
% 
% Simple demonstration for non-linear IMM using the following models:
%  1. Standard Wiener process velocity model
%  2. Coordinated turn model
%
% The measurement model is linear, which gives noisy measurements 
% of target's position.
% 
% Copyright (C) 2007-2008 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

print_figures = 1;
save_figures = 0;

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

% Process noise variance
KF_q1 = .05;
KF_Qc1 = diag([KF_q1 KF_q1]);

% Discretization of the continous-time system.
[KF_A1,KF_Q1] = lti_disc(F{1},L{1},KF_Qc1,dt);


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

% Check the derivatives
%der_check(a_func{2}, A{2}, 1, [1 1 1 1 0]',{dt});


% Measurement models.
H{1} = [1 0 0 0;
        0 1 0 0];

H{2} = [1 0 0 0 0;
        0 1 0 0 0];

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

% Noise variances of the measurement models

% Model 1 
r1 = .05;
r2 = .05;
R{1} = diag([r1 r2]);

% Model 2
r1 = .05;
r2 = .05;
R{2} = diag([r1 r2]);

% Generate the measurements.
Y = zeros(hdims,n);
for i = 1:n
    Y(:,i) = H{mstate(i)}*X_r(ind{mstate(i)},i) + gauss_rnd(zeros(size(Y,1),1), R{mstate(i)});
end

ac = floor(n/2)+1:floor(n/2)+2;
clf; clc;
%fprintf('Filtering with KF...');

if print_figures
    h = plot(X_r(1,:),X_r(2,:),'-g',...
             Y(1,:),Y(2,:),'.',... 
             X_r(1,1),X_r(2,1),'ro','MarkerSize',12);    
    legend('Real trajectory',...
           'Measurements',...
           'Starting position');
    set(h,'markersize',5);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',8);
    xlim([-1 7.5]) 
    ylim([-3.5 3.5]) 
    %title('Trajectory of the object');
    if save_figures
        print('-dpsc','eimm1_trajectory.ps');
    end
    pause
end


m = [0 0 0 -1 0]';
P = diag([10.1 10.1 1.1 1.1 1]);

%% Space for the estimates.

% KF with model 1
KF_MM = zeros(size(A{1},1), size(Y,2));
KF_PP = zeros(size(A{1},1), size(A{1},1), size(Y,2));

% IMM
IMM1_MM = zeros(size(m,1), size(Y,2));
IMM1_PP = zeros(size(m,1), size(m,1), size(Y,2));
IMM1_MM_i = cell(nmodels,n);
IMM1_PP_i = cell(nmodels,n);
IMM1_MU = zeros(nmodels,size(Y,2));

%%% Initial estimates %%%

% KF with model 1
KF_M = [0 0 0 -1]';
KF_P = diag([1.1 1.1 0.1 0.1]);

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
    % KF with model 1
    [KF_M,KF_P] = kf_predict(KF_M,KF_P,KF_A1,KF_Q1);
    [KF_M,KF_P] = kf_update(KF_M,KF_P,Y(:,i),H{1},R{1});
    KF_MM(:,i)   = KF_M;
    KF_PP(:,:,i) = KF_P;
    
    % EKF based IMM
    [x_p1,P_p1,c_j1] = eimm_predict(x_ip1,P_ip1,mu_ip1,p_ij,ind,dims,A,a_func,a_param,Q);
    [x_ip1,P_ip1,mu_ip1,m1,P1] = imm_update(x_p1,P_p1,c_j1,ind,dims,Y(:,i),H,R);

    EIMM_MM(:,i)   = m1;
    EIMM_PP(:,:,i) = P1;
    EIMM_MU(:,i)   = mu_ip1';
    EIMM_MM_i(:,i) = x_ip1';
    EIMM_PP_i(:,i) = P_ip1';
    
    % UKF based IMM
    [x_p2,P_p2,c_j2] = uimm_predict(x_ip2,P_ip2,mu_ip2,p_ij,ind,dims,A,a_func,a_param,Q);    
    [x_ip2,P_ip2,mu_ip2,m2,P2] = imm_update(x_p2,P_p2,c_j2,ind,dims,Y(:,i),H,R);

    UIMM_MM(:,i)   = m2;
    UIMM_PP(:,:,i) = P2;
    UIMM_MU(:,i)   = mu_ip2';
    UIMM_MM_i(:,i) = x_ip2';
    UIMM_PP_i(:,i) = P_ip2';

    if print_figures
        % Plot the estimates obtained so far
        plot(Y(1,1:i),Y(2,1:i),'k.',...
             KF_MM(1,1:i),KF_MM(2,1:i),'y-',...
             EIMM_MM(1,1:i),EIMM_MM(2,1:i),'r-',...
             UIMM_MM(1,1:i),UIMM_MM(2,1:i),'b-',...         
             X_r(1,1:i),X_r(2,1:i),'g-');
        xlim([-1 7.5]) 
        ylim([-3.5 3.5]) 
        
        drawnow
    end
end

% Smooth with EKF based IMM
[SMI_1,SPI_1,SMI_i_1,SPI_i_1,MU_S_1] = eimm_smooth(EIMM_MM,EIMM_PP,EIMM_MM_i,EIMM_PP_i,EIMM_MU,p_ij,mu_0j,ind,dims,A,ia_func,a_param,Q,R,H,[],[],Y);

% Smooth with UKF based IMM
[SMI_2,SPI_2,SMI_i_2,SPI_i_2,MU_S_2] = uimm_smooth(UIMM_MM,UIMM_PP,UIMM_MM_i,UIMM_PP_i,UIMM_MU,p_ij,mu_0j,ind,dims,A,ia_func,a_param,Q,R,H,[],[],Y);

% Smooth with KF
[SM1,SP1] = erts_smooth1(KF_MM,KF_PP,KF_A1,KF_Q1);

% Calculate the MSEs
MSE_KF1_1 = mean((X_r(1,:)-KF_MM(1,:)).^2);
MSE_KF1_2 = mean((X_r(2,:)-KF_MM(2,:)).^2);
MSE_KF1 = 1/2*(MSE_KF1_1 + MSE_KF1_2);

MSE_KS1_1 = mean((X_r(1,:)-SM1(1,:)).^2);
MSE_KS1_2 = mean((X_r(2,:)-SM1(2,:)).^2);
MSE_KS1 = 1/2*(MSE_KS1_1 + MSE_KS1_2);

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
fprintf('KF1-MSE = %.4f\n',MSE_KF1);
fprintf('KS1-MSE = %.4f\n',MSE_KS1);
fprintf('EIMM-MSE = %.4f\n',MSE_EIMM);
fprintf('EIMMS-MSE = %.4f\n',MSE_EIMMS);
fprintf('UIMM-MSE = %.4f\n',MSE_UIMM);
fprintf('UIMMS-MSE = %.4f\n',MSE_UIMMS);


% Plot the final filtering and smoothing results
if print_figures
    h = plot(X_r(1,:),X_r(2,:),'g-',...         
             KF_MM(1,:),KF_MM(2,:),'-k',...
             SM1(1,:),SM1(2,:),'--k',...
             EIMM_MM(1,:),EIMM_MM(2,:),'-r',...             
             SMI_1(1,:),SMI_1(2,:),'--r',...
             UIMM_MM(1,:),UIMM_MM(2,:),'-b',...
             SMI_2(1,:),SMI_2(2,:),'--b');
    legend('True trajectory',...
           'KF',...
           'RTS',...
           'IMM-EKF',...
           'IMM-EKS',...
           'IMM-UKF',...
           'IMM-UKS');
    %title('Estimates produced by Kalman filter using the model 1.');
    set(h,'markersize',2);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',8);
    xlim([-1 7.5]) 
    ylim([-3.5 3.5]) 
    if save_figures
        print('-dpsc','eimm1_1.ps');
    end
    pause
    
    % Determine the real model probabilities
    p_models = zeros(nmodels,n);
    I1 = find(mstate == 1);
    p_models(1,I1) = 1;
    I2 = find(mstate == 2);
    p_models(2,I2) = 1;

    % Plot model 2 probability for each step
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
        print('-dpsc','eimm1_2.ps');
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
        print('-dpsc','eimm1_3.ps');
    end

end

