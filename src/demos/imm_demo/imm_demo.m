%IMM_DEMO  Tracking a Target with Simple Manouvers demonstration
% 
% Simple demonstration for linear IMM using the following models:
%  1. Standard Wiener process velocity model
%  2. Standard Wiener process acceleration model
%
% The measurement model is linear, which gives noisy measurements 
% of target's position.
% 
% Copyright (C) 2007-2008 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.


% Simple demonstration for IMM using velocity and acceleration models

% Save figures
print_figures = 1;
save_figures = 1;


% Dimensionality of the state space
dims = 6;

nmodels = 2;

% Space for models
ind = cell(1,2);
A = cell(1,2);
Q = cell(1,2);
H = cell(1,2);
R = cell(1,2);

% Stepsize
dt = 0.1;

ind{2} = [1 2 3 4 5 6]';

% Transition matrix for the continous-time acceleration model.
F2 = [0 0 1 0 0 0;
      0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1;
      0 0 0 0 0 0;
      0 0 0 0 0 0];

% Noise effect matrix for the continous-time system.
L2 = [0 0;
      0 0;
      0 0;
      0 0;
      1 0;
      0 1];

% Process noise variance
q2 = 1.00;
Qc2 = diag([q2 q2]);

% Discretization of the continous-time system.
[A{2},Q{2}] = lti_disc(F2,L2,Qc2,dt);

ind{1} = [1 2 3 4]';

% Transition matrix for the continous-time velocity model.
F1 = [0 0 1 0;
      0 0 0 1;
      0 0 0 0;
      0 0 0 0];

% Noise effect matrix for the continous-time system.
L1 = [0 0;
      0 0;
      1 0;
      0 1];

% Process noise variance
q1 = .01;
Qc1 = diag([q1 q1]);

% Discretization of the continous-time system.
[A{1},Q{1}] = lti_disc(F1,L1,Qc1,dt);


H{1} = [1 0 0 0;
        0 1 0 0];

% Measurement models.
H{2} = [1 0 0 0 0 0;
        0 1 0 0 0 0];

hdims = 2;

% Variance in the measurements.
r1 = .1;
r2 = .1;
R{1} = diag([r1 r2]);

r1 = .1;
r2 = .1;
R{2} = diag([r1 r2]);


% Generate the data.
n = 200;
Y = zeros(hdims,n);
X_r = zeros(dims,n);
X_r(:,1) = [0 0 0 -1 0 0]';
mstate = zeros(1,n);

mu_ip = [0.9 0.1];
mu_0j = mu_ip;
%p_ij = [0.65 0.35;
%        0.10 0.90];

p_ij = [0.98 0.02;
        0.02 0.98];

% $$$ % Generate random mode transitions
% $$$ mstate(1) = 2;
% $$$ for i = 2:n
% $$$     r = rand;
% $$$     for j = size(p_ij,1)-1;
% $$$         if r < p_ij(mstate(i-1),j);
% $$$             mstate(i) = j;
% $$$         end
% $$$     end
% $$$     if mstate(i) == 0
% $$$         mstate(i) = size(p_ij,1);
% $$$     end
% $$$ end

% Forced mode transitions 
mstate(1:50) = 1;
mstate(51:70) = 2;
mstate(71:120) = 1;
mstate(121:150) = 2;
mstate(151:200) = 1;

% Acceleration model
for i = 2:n
    st = mstate(i);
    X_r(ind{st},i) = A{st}*X_r(ind{st},i-1) + gauss_rnd(zeros(size(A{st},1),1), Q{st});
end

% Generate the measurements.
for i = 1:n
    Y(:,i) = H{mstate(i)}*X_r(ind{mstate(i)},i) + gauss_rnd(zeros(size(Y,1),1), R{mstate(i)});
end

ac = floor(n/2)+1:floor(n/2)+2;
clf; clc;
%fprintf('Filtering with KF...');

% $$$ plot(X_r(1,:),X_r(2,:),Y(1,:),Y(2,:),'.',X_r(1,1),...
% $$$      X_r(2,1),'ro','MarkerSize',12);
% $$$ legend('Real trajectory', 'Measurements');
% $$$ title('Position');

m = [0 0 0 -1 0 0]';
P = diag([0.1 0.1 0.1 0.1 0.5 0.5]);

%%%% Space for the estimates %%%% 

% KF using model 1 
MM1 = zeros(size(A{1},1), size(Y,2));
PP1 = zeros(size(A{1},1), size(A{1},1), size(Y,2));

% KF using model 2
MM2 = zeros(size(A{2},1), size(Y,2));
PP2 = zeros(size(A{2},1), size(A{2},1), size(Y,2));

% Overall estimates of IMM filter
MM = zeros(size(m,1), size(Y,2));
PP = zeros(size(m,1), size(m,1), size(Y,2));
% Model-conditioned estimates of IMM
MM_i = cell(2,n);
PP_i = cell(2,n);

% Model probabilities 
MU = zeros(2,size(Y,2));

%%%% Prior estimates %%%%

% KF1
M1 = [0 0 0 -1]';
P1 = diag([0.1 0.1 0.1 0.1]);

% KF2
M2 = [0 0 0 -1 0 0]';
P2 = diag([0.1 0.1 0.1 0.1 0.5 0.5]);


% IMM
x_ip{1} = [0 0 0 -1]';
P_ip{1} = diag([0.1 0.1 0.1 0.1]);
x_ip{2} = [0 0 0 -1 0 0]';
P_ip{2} = diag([0.1 0.1 0.1 0.1 0.5 0.5]);


% Filtering steps.
for i = 1:size(Y,2)
    % KF with model 1
    [M1,P1] = kf_predict(M1,P1,A{1},Q{1});
    [M1,P1] = kf_update(M1,P1,Y(:,i),H{1},R{1});
    MM1(:,i) = M1;
    PP1(:,:,i) = P1;

    % KF with model 2
    [M2,P2] = kf_predict(M2,P2,A{2},Q{2});
    [M2,P2] = kf_update(M2,P2,Y(:,i),H{2},R{2});
    MM2(:,i) = M2;
    PP2(:,:,i) = P2;
    
    % IMM
    [x_p,P_p,c_j] = imm_predict(x_ip,P_ip,mu_ip,p_ij,ind,dims,A,Q);
    [x_ip,P_ip,mu_ip,m,P] = imm_update(x_p,P_p,c_j,ind,dims,Y(:,i),H,R);
    MM(:,i)   = m;
    PP(:,:,i) = P;
    MU(:,i)   = mu_ip';
    MM_i(:,i) = x_ip';
    PP_i(:,i) = P_ip';
    
    % Plot the estimates
    if print_figures
        plot(X_r(1,1:i),X_r(2,1:i),'g-',...         
             MM(1,1:i),MM(2,1:i),'r-',...
             MM1(1,1:i),MM1(2,1:i),'y-',...
             MM2(1,1:i),MM2(2,1:i),'k-',...
             Y(1,1:i),Y(2,1:i),'ko');
        drawnow
    end
end

% Smooth the filtered estimates
[SM1,SP1,SS1] = rts_smooth(MM1,PP1,A{1},Q{1});
[SM2,SP2,SS2] = rts_smooth(MM2,PP2,A{2},Q{2});
[SM3,SP3,SM3_i,SP3_i,MU_S] = imm_smooth(MM,PP,MM_i,PP_i,MU,p_ij,mu_0j,ind,dims,A,Q,R,H,Y);

% Calculate the MSEs
MSE_KF1_1 = mean((X_r(1,:)-MM1(1,:)).^2);
MSE_KF1_2 = mean((X_r(2,:)-MM1(2,:)).^2);
MSE_KF1 = 1/2*(MSE_KF1_1 + MSE_KF1_2);

MSE_KS1_1 = mean((X_r(1,:)-SM1(1,:)).^2);
MSE_KS1_2 = mean((X_r(2,:)-SM1(2,:)).^2);
MSE_KS1 = 1/2*(MSE_KS1_1 + MSE_KS1_2);

MSE_KF2_1 = mean((X_r(1,:)-MM2(1,:)).^2);
MSE_KF2_2 = mean((X_r(2,:)-MM2(2,:)).^2);
MSE_KF2 = 1/2*(MSE_KF2_1 + MSE_KF2_2);

MSE_KS2_1 = mean((X_r(1,:)-SM2(1,:)).^2);
MSE_KS2_2 = mean((X_r(2,:)-SM2(2,:)).^2);
MSE_KS2 = 1/2*(MSE_KS2_1 + MSE_KS2_2);

MSE_IMM1 = mean((X_r(1,:)-MM(1,:)).^2);
MSE_IMM2 = mean((X_r(2,:)-MM(2,:)).^2);
MSE_IMM = 1/2*(MSE_IMM1 + MSE_IMM2);

MSE_IMMS1 = mean((X_r(1,:)-SM3(1,:)).^2);
MSE_IMMS2 = mean((X_r(2,:)-SM3(2,:)).^2);
MSE_IMMS = 1/2*(MSE_IMMS1 + MSE_IMMS2);

fprintf('Mean square errors of position estimates:\n');
fprintf('KF1-RMSE = %.4f\n',MSE_KF1);
fprintf('KS1-RMSE = %.4f\n',MSE_KS1);
fprintf('KF2-RMSE = %.4f\n',MSE_KF2);
fprintf('KS2-RMSE = %.4f\n',MSE_KS2);
fprintf('IMM-RMSE = %.4f\n',MSE_IMM);
fprintf('IMMS-RMSE = %.4f\n',MSE_IMMS);

% Plot the final filtering and smoothing results
if print_figures
    h = plot(Y(1,:),Y(2,:),'ko',...
             X_r(1,:),X_r(2,:),'g-',...         
             MM1(1,:),MM1(2,:),'r-',...
             SM1(1,:),SM1(2,:),'b-');
    legend('Measurement',...
           'True trajectory',...
           'Filtered',...
           'Smoothed');
    title('Estimates produced by Kalman filter using the model 1.');
    set(h,'markersize',2);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',6);
    if save_figures
        print('-depsc','imm1.eps');
    end
    
    h = plot(Y(1,:),Y(2,:),'ko',...
             X_r(1,:),X_r(2,:),'g-',...         
             MM2(1,:),MM2(2,:),'r-',...
             SM2(1,:),SM2(2,:),'b-');
    legend('Measurement',...
           'True trajectory',...
           'Filtered',...
           'Smoothed');
    title('Estimates produced by Kalman filter using the model 2.');
    set(h,'markersize',2);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',6);
    if save_figures
        print('-depsc','imm2.eps');
    end
    
    h = plot(Y(1,:),Y(2,:),'ko',...
             X_r(1,:),X_r(2,:),'g-',...         
             MM(1,:),MM(2,:),'r-',...
             SM2(1,:),SM2(2,:),'b-');
    legend('Measurement',...
           'True trajectory',...
           'Filtered',...
           'Smoothed');
    title('Estimates produced by IMM-filter.')
    set(h,'markersize',2);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',6);
    if save_figures
        print('-depsc','imm3.eps');
    end
    
    h = plot(1:n,2-mstate,'g--',...
             1:n,MU(1,:)','r-',...  
             1:n,MU_S(1,:)','b-');         
    legend('True',...
           'Filtered',...
           'Smoothed');
    title('Probability of model 1');
    ylim([-0.1,1.1]);
    set(h,'markersize',2);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',6);
    if save_figures
        print('-depsc','imm4.eps');
    end
end