% Demonstration for ADF methods using a random sine signal model. 
%
% A Very simple demonstration for ADF methods, which are
% used to track a random single-component sinusoid signal,
% which is modelled as x_k = a_k*sin(\theta_k), dtheta/dt = omega_k.
% The signal is also filtered with unscented Kalman filter (UKF) for
% comparison.
%
% Copyright (C) 2007, 2010 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

clc;
disp('Filtering the signal with EKF...');

% Generate the real signal.clc;
disp('Filtering the signal with EKF...');

save_plots = 0;

% Measurement model and it's derivative
h_func = @sine_h;
dh_dx_func = @sine_dh_dx;
d2h_dx2_func = @sine_d2h_dx2;

% Initial values for the signal.
f = 0;
w = 10;
a = 1;
  
% Number of samples and stepsize.
d = 5;
n = 500;
dt = d/n;
xx = 1:n;

% Check the derivative of the measurement function.
der_check(h_func, dh_dx_func, 1, [f w a]');

% Dynamic state transition matrix in continous-time domain.
F = [0 1 0;
     0 0 0;
     0 0 0];
  
% Noise effect matrix in continous-time domain.
L = [0 0;
     1 0;
     0 1];
  
% Spectral power density of the white noise.
q1 = 2*0.2;
q2 = 2*0.1;
Qc = diag([q1 q2]);
  
% Discretize the plant equation.
[A,Q] = lti_disc(F,L,Qc,dt);
  
% State prior
M0 = [f w a]';
P0 = diag([.3 .3 .3]);

% Generate the real signal.
X = zeros(3, n);
x = gauss_rnd(M0,P0,1);

for i = 1:n
   x = A*x + gauss_rnd([0 0 0]', Q);
   X(:,i) = x;
end  

X = zeros(3, n);
x = gauss_rnd(M0,P0,1);

for i = 1:n
   x = A*x + gauss_rnd([0 0 0]', Q);
   X(:,i) = x;
end  
  
% Generate the observations with Gaussian noise.
sd = 1;
R = sd^2;

Y_real = feval(h_func,X);     
Y = Y_real + gauss_rnd(0,R,n);
  
plot(xx,Y,'.'); 
hold on;
plot(xx,Y_real,'k-','LineWidth',2)
hold off;

%% Run the various filters and smoothers

% First setup the needed parameters
tr_method  = cell(1,6);
tr_param_f = cell(1,6);
tr_param_h = cell(1,6);
tr_name    = cell(1,6);

% EKF1
tr_method{1}  = @lin_transform;
tr_param_f{1} = {};
tr_param_h{1} = {dh_dx_func};
tr_name{1}    = {'EKF1' 'ERTS1'};

% EKF2
tr_method{2}  = @quad_transform;
tr_param_f{2} = {};
tr_param_h{2} = {dh_dx_func d2h_dx2_func};
tr_name{2}    = {'EKF2' 'ERTS2'};

% UKF
tr_method{3}  = @ut_transform;
tr_param_f{3} = {};
tr_param_h{3} = {};
tr_name{3}    = {'UKF' 'URTS'};

% Gauss-Hermite Kalman filter
tr_method{4}  = @gh_transform;
tr_param_f{4} = {3};
tr_param_h{4} = {3};
tr_name{4}    = {'GHKF' 'GHRTS'};

% Cubature Kalman filter
tr_method{5}  = @ckf_transform;
tr_param_f{5} = {};
tr_param_h{5} = {};
tr_name{5}    = {'CKF' 'CRTS'};

% Central-difference Kalman filter
tr_method{6}  = @cd_transform;
tr_param_f{6} = {sqrt(3)};
tr_param_h{6} = {sqrt(3)};
tr_name{6}    = {'CDKF' 'CDRTS'};


%%% Filtering and smoothing

% Space for means and covariances
ntr = length(tr_method);
MM_ADF = cell(1,ntr);
PP_ADF = cell(1,ntr);
MMS_ADRTS = cell(1,ntr);
PPS_ADRTS = cell(1,ntr);

FMSE = zeros(ntr,length(M0));
SMSE = zeros(ntr,length(M0));


for i = 1:ntr
    fprintf(['Filtering and smoothing with ' tr_name{i}{1} '/' tr_name{i}{2} '...']);
    
    % Initialize the mean and covariance
    M = M0;
    P = P0;
    
    MM = zeros(size(M,1),size(Y,2));
    PP = zeros(size(M,1),size(M,1),size(Y,2));
    
    % Filtering loop
    for k = 1:size(Y,2)
        % The dynamic model is linear, so ordinary Kalman filter prediction
        % is sufficient
        [M,P] = kf_predict(M,P,A,Q);
        [M,P] = adf_update(M,P,Y(:,k),h_func,R,s,tr_method{i},tr_param_h{i});
        MM(:,k)   = M;
        PP(:,:,k) = P;
    end
    
    % Smooth
    [MMS, PPS] = rts_smooth(MM,PP,A,Q);
    
    % Save estimates
    MM_ADF{i} = MM;
    PP_ADF{i} = PP;
    MMS_ADRTS{i} = MMS;
    PPS_ADRTS{i} = PPS;
    
    % Calculate RMSE values
    FMSE(i,:) = sqrt(mean((X-MM).^2,2));
    SMSE(i,:) = sqrt(mean((X-MMS).^2,2));
    
    fprintf('Done!\n');
end

for i = 1:ntr 
    clf; clc;
    disp(['The filtering and smoothing results with ' tr_name{i}{1} ' and ' tr_name{i}{2} ' are now displayed'])

    % Project the estimates to measurement space 
    Y_fm = feval(h_func, MM_ADF{i}); 
    Y_sm = feval(h_func, MMS_ADRTS{i});
    h1 = plot(xx,Y,'k.'); hold on;
    h2 = plot(xx,Y_real,'k-','LineWidth',3); 
    h3 = plot(xx,Y_fm,'r-','LineWidth',2); 
    h4 = plot(xx,Y_sm,'b-','LineWidth',2); hold off;
    
    legend([h1;h2;h3;h4],'Measurements','Real signal', 'Filtered estimate', 'Smoothed estimate');
    xlim([0 ceil(max(xx))]);
    title(['Estimating a random Sine signal with ' tr_name{i}{1}]);

    disp(' ');
    disp('<push any button to continue>');
    pause
end


%%% Print errors
clc;
disp('RMS errors:');
for i = 1:ntr
    fprintf([tr_name{i}{1} sprintf(' = %.4f\n',FMSE(i))]);
    fprintf([tr_name{i}{2} sprintf(' = %.4f\n',SMSE(i))]);
end
  

