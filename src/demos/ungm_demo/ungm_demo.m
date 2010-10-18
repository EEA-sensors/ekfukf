%% Demonstration of univariate nonstationary growth model (UNGM)
%
%  Description:
%    In this example various different non-linear filters and smoothers are
%    applied to the univariate nonstationary growth model (UNGM). The
%    filters used in this demonstration are:
%      * Extended Kalman filter
%      * Unscented Kalman filter (augmented and non-augmented forms)
%      * Bootstrap filter
%      * Gauss-Hermite Kalman filter (degree 10)
%      * Cubature Kalman filter
%    Additionally, the corresponding Rauch-Tung-Striebel smoother results
%    are also presented.
%
%  References:
%    Refer to the Toolbox documentation for details on the model.
%
%  See also:
%   adf_predict, adf_update, adrts_smooth, lin_transform, ut_transform
%   cd_transform, gh_transform, ckf_transform
%
%  Author:
%    Copyright (C) 2007, 2010 Jouni Hartikainen,
%                  2010, Arno Solin
%
%  Licence:
%    This software is distributed under the GNU General Public
%    Licence (version 2 or later); please refer to the file
%    Licence.txt, included with the software, for details.

%% Set up the model parameters

silent = 0;

clf; clc;
disp('This is a demonstration for Assumed Density Filters and Smoothers');
disp('by using the univariate nonstationary growth model (UNGM).');
disp(' ');

% Handles to dynamic and measurement model functions,
% and their derivatives
f_func = @ungm_f;
h_func = @ungm_h;
df_func = @ungm_df_dx;
dh_func = @ungm_dh_dx;
d2f_func = @ungm_d2f_dx2;
d2h_func = @ungm_d2h_dx2;

% Number of samples
n = 500;

% Initial state and covariance
M_0 = .1;
P_0 = 1;

% Space for measurements
Y = zeros(1,n);

% Strengths of perturbations
Q = 1;
R = 1;

fprintf('Generating real states and measurements...');
% Generate the true states with process noise
X = zeros(1,n);
X(1) = ungm_f(M_0,1) + gauss_rnd(0,P_0,1);
%X(1) = ungm_f(x_0,1) + gauss_rnd(0,u_n,1);
for i = 2:n
    X(i) = feval(f_func,X(i-1),i) + gauss_rnd(0,Q,1);
end

% Generate the observations with measurement noise
for i = 1:n
    Y(i) = feval(h_func,X(i)) + gauss_rnd(0,R,1);
end

Y_r = feval(h_func,X);

fprintf('Done!\n');

% Parameters for dynamic model. Used by smoothers.
params = cell(size(Y,2));
for i = 1:size(Y,2)
    params{i} = i+1;
end


%% Run the various filters and smoothers
tr_method  = cell(1,6);
tr_param_f = cell(1,6);
tr_param_h = cell(1,6);
tr_name    = cell(1,6);

% EKF
tr_method{1}  = @lin_transform;
tr_param_f{1} = {df_func};
tr_param_h{1} = {dh_func};
tr_name{1}    = {'EKF' 'ERTS'};

% UKF
tr_method{2}  = @ut_transform;
tr_param_f{2} = {};
tr_param_h{2} = {};
tr_name{2}    = {'UKF' 'URTS'};

% Gauss-Hermite Kalman filter
tr_method{3}  = @gh_transform;
tr_param_f{3} = {10};
tr_param_h{3} = {10};
tr_name{3}    = {'GHKF' 'GHRTS'};

% Cubature Kalman filter
tr_method{4}  = @ckf_transform;
tr_param_f{4} = {};
tr_param_h{4} = {};
tr_name{4}    = {'CKF' 'CRTS'};

% Central-difference Kalman filter
tr_method{5}  = @cd_transform;
tr_param_f{5} = {sqrt(3)};
tr_param_h{5} = {sqrt(3)};
tr_name{5}    = {'CDKF' 'CDRTS'};

% Monte Carlo Kalman filter
% (can be numerically unstable with low amount of samples)
tr_method{6}  = @mc_transform;
tr_param_f{6} = {5000};
tr_param_h{6} = {5000};
tr_name{6}    = {'MCKF' 'MCRTS'};

%%% Filtering and smoothing

% Space for means and covariances
ntr = length(tr_method);
MM_ADF = cell(1,ntr);
PP_ADF = cell(1,ntr);
MMS_ADRTS = cell(1,ntr);
PPS_ADRTS = cell(1,ntr);

FMSE = zeros(1,ntr);
SMSE = zeros(1,ntr);


for i = 1:ntr
    fprintf(['Filtering and smoothing with ' tr_name{i}{1} '/' tr_name{i}{2} '...']);
    
    % Initialize the mean and covariance
    M = M_0;
    P = P_0;
    
    MM = zeros(size(M,1),size(Y,2));
    PP = zeros(size(M,1),size(M,1),size(Y,2));
    
    % Filtering loop
    for k = 1:size(Y,2)
        if k > 1
            [M,P] = adf_predict(M,P,f_func,Q,k,tr_method{i},tr_param_f{i});
        end
        [M,P] = adf_update(M,P,Y(:,k),h_func,R,[],tr_method{i},tr_param_h{i});
        MM(:,k)   = M;
        PP(:,:,k) = P;
    end
    
    % Smooth
    [MMS, PPS] = adrts_smooth(MM,PP,f_func,Q,params,tr_method{i},tr_param_f{i},0);
    
    % Save estimates
    MM_ADF{i} = MM;
    PP_ADF{i} = PP;
    MMS_ADRTS{i} = MMS;
    PPS_ADRTS{i} = PPS;
    
    % Calculate MSE values
    FMSE(i) = sum((X-MM).^2)/n;
    SMSE(i) = sum((X-MMS).^2)/n;
    
    fprintf('Done!\n');
end

fprintf('Filtering with bootstrap filter...');

% Filtering with bootstrap filter
M = M_0;
P = P_0;
n_particles = 1000;
SX = gauss_rnd(M,P,n_particles);
SXX = zeros(n_particles,size(Y,2));

MM_BSF = zeros(size(M,1),size(Y,2));
PP_BSF = zeros(size(M,1),size(M,1),size(Y,2));

% Filtering loop for bootstrap filter
for k = 1:size(Y,2)
    SX = feval(f_func,SX,k) + gauss_rnd(0,1,size(SX,2));
    W  = gauss_pdf(Y(:,k),feval(h_func,SX,k),1);
    ind = resampstr(W);
    SX = SX(:,ind);
    M = mean(SX);
    P = var(SX);
    MM_BSF(:,k)   = M;
    PP_BSF(:,:,k) = P;
    SXX(:,k) = SX';
end

BSF_MSE = sum((X-MM_BSF).^2)/n;

fprintf('Done!\n');


%% Visualize results
if ~silent
    subplot(3,1,1);
    plot(1:100,X(1:100),'-kx',1:100,MM_ADF{1}(1:100),'--bo')
    title([tr_name{1}{1} ' filtering result']);
    xlim([0 100]);
    ylim([-20 20]);
    legend('Real signal', [tr_name{1}{1} ' estimate']);
    
    subplot(3,1,2);
    plot(1:100,X(1:100),'-kx',1:100,MM_ADF{3}(1:100),'--bo')
    title([tr_name{3}{1} ' filtering result']);
    xlim([0 100]);
    ylim([-20 20]);
    legend('Real signal', [tr_name{3}{1} ' estimate']);
    
    subplot(3,1,3);
    plot(1:100,X(1:100),'-kx',1:100,MM_BSF(1:100),'--bo')
    title('Bootstrap filtering result')
    xlim([0 100]);
    ylim([-20 20]);
    legend('Real signal','BSF estimates')
    
    disp(' ');
    disp('First 100 values of the filtering results with EKF, GHKF and BS are now displayed.')
    disp(' ');
    disp('<press any key to see the estimation error of EKF, GHKF and BS>');
    pause;
    
    subplot(3,1,1);
    E = X-MM_ADF{1};
    PP_EKF = squeeze(PP_ADF{1});
    plot(1:100,E(1:100),'-kx',1:100,3*sqrt(PP_EKF(1:100))','--r',1:100,-3*sqrt(PP_EKF(1:100))','--r');
    title([tr_name{1}{1} ' error']);
    legend([tr_name{1}{1} ' error'],'3\sigma interval');
    
    subplot(3,1,2);
    E = X-MM_ADF{3};
    PP_GHKF = squeeze(PP_ADF{3});
    plot(1:100,E(1:100),'-kx',1:100,3*sqrt(PP_GHKF(1:100))','--r',1:100,-3*sqrt(PP_GHKF(1:100))','--r');
    title([tr_name{3}{1} ' error']);
    legend([tr_name{3}{1} ' error'],'3\sigma interval');
    
    subplot(3,1,3);
    E = X-MM_BSF;
    PP_BSF = squeeze(PP_BSF);
    plot(1:100,E(1:100),'-kx',1:100,3*sqrt(PP_BSF(1:100))','--r',1:100,-3*sqrt(PP_BSF(1:100))','--r');
    title('BSF error');
    legend('BSF error','3\sigma interval')

    clc;
    disp('The estimation errors of EKF, GHKF and BS are now displayed.');
    disp(' ');
    disp('<press any key to see the filtering results of UKF1 and UKF2>');
end

%%% Print errors
disp('MS errors:');
for i = 1:ntr
    fprintf([tr_name{i}{1} sprintf(' = %.4f\n',FMSE(i))]);
    fprintf([tr_name{i}{2} sprintf(' = %.4f\n',SMSE(i))]);
end
fprintf('BSF = %.4f\n',BSF_MSE)