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
%      * Central-Difference Kalman filter
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
%  Authors:
%    Copyright (C) 2002, 2003 Simo Särkkä
%                  2007, 2010 Jouni Hartikainen
%                  2010       Arno Solin
%
%  Licence:
%    This software is distributed under the GNU General Public
%    Licence (version 2 or later); please refer to the file
%    Licence.txt, included with the software, for details.

%% Set parameters

silent = 0;


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
s = [S1 S2];
sd = 0.05;
dt = 0.01;
R  = sd.^2*eye(size(s,2));

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


clc;  clf;
disp(['In this demonstration we track a moving object with two sensors, ',...
    'which gives only bearings of the object with respect to sensors position. ',...
    'The state of the system is estimated with 1st and 2nd order EKF, and UKF.'])
disp(' ');

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
    FMSE(i) = sqrt(mean((X(1,:)-MM(1,:)).^2+(X(2,:)-MM(2,:)).^2));
    SMSE(i) = sqrt(mean((X(1,:)-MMS(1,:)).^2+(X(2,:)-MMS(2,:)).^2));
    
    fprintf('Done!\n');
end

if ~silent
    for i = 1:ntr
        plot(X(1,:),X(2,:),'k-',...
            MM_ADF{i}(1,:),MM_ADF{i}(2,:),'b--',...
            MMS_ADRTS{i}(1,:),MMS_ADRTS{i}(2,:),'r-.',...
            S1(1),S1(2),'k^',S2(1),S2(2),'k^');
        axis([-1.5 1.5 -2.5 1.5]);
        legend('Real trajectory',...
            [tr_name{i}{1} ' estimate'],...
            [tr_name{i}{2} ' estimate'],...
            'Positions of sensors',...
            'Location', 'NorthWest');
        title(['Filtering and smoothing result with ' tr_name{i}{1} '/' tr_name{i}{2}]);
        clc;
        disp(['Results with ' tr_name{i}{1} ' and ' tr_name{i}{2} ' are now displayed.'])
        disp(' ');
        disp('<push any key to continue>')
        pause;
    end
end


%%% Print errors
disp('RMS errors:');
for i = 1:ntr
    fprintf([tr_name{i}{1} sprintf(' = %.4f\n',FMSE(i))]);
    fprintf([tr_name{i}{2} sprintf(' = %.4f\n',SMSE(i))]);
end

  