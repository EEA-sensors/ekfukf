% A Very simple demonstration for extended Kalman filter (EKF), which is
% used to track a random single-component sinusoid signal,
% which is modelled as x_k = a_k*sin(\theta_k), dtheta/dt = omega_k.
% The signal is also filtered with unscented Kalman filter (UKF) for
% comparison.
%
% Copyright (C) 2007 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

%set(0,'DefaultAxesColorOrder',[0 0 0])
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultLineLineWidth',1)

clf
a4
%figscale=[1 0.28];
figscale=[2 1];
set(gcf,'pos',[27.4 3.8 10.795.*figscale(1) 15.26.*figscale(2)]);
set(gcf,'paperposition',[5.1 6.9 10.795.*figscale(1) 15.26.*figscale(2)]);
set(gcf,'DefaultAxesFontSize',10)
set(gcf,'DefaultTextFontSize',10)
set(gca,'Box','on');

clc;
disp('Filtering the signal with EKF...');

% Measurement model and it's derivative
f_func = @ekf_sine_f;
h_func = @ekf_sine_h;
dh_dx_func = @ekf_sine_dh_dx;

% Initial values for the signal.
f = 0;
w = 10;
a = 1;
  
% Number of samples and stepsize.
d = 5;
n = 500;
dt = d/n;
ti = f+dt:dt:d


% Check the derivative of the measurement function.
der_check(h_func, dh_dx_func, 1, [f w a]', []);

% Dynamic state transition matrix in continous-time domain.
F = [0 1 0;
     0 0 0;
     0 0 0];
  
% Noise effect matrix in continous-time domain.
L = [0 0;
     1 0;
     0 1];
  
% Spectral power density of the white noise.
q1 = 0.2;
q2 = 0.1;
Qc = diag([q1 q2]);
  
% Discretize the plant equation.
[A,Q] = lti_disc(F,L,Qc,dt);
  
% Generate the real signal.
X = zeros(3, n);
X(:,1) = [f w a]';
for i = 2:n
   X(:,i) = A*X(:,i-1) + gauss_rnd([0 0 0]', Q);
end  
  
% Generate the observations with Gaussian noise. Variance
% of the noise is set greater in samples (1/2n,3/4n) than in other
% samples.
sd = 1;
R = sd^2;

Y = zeros(1,n);
Y_real = feval(h_func,X);     
Y = Y_real + gauss_rnd(0,R,n);
  
plot(X(1,:),Y,'.', X(1,:), Y_real)
  
% Initial guesses for the state mean and covariance.
M = [f w a]';
P = diag([3 3 3]);    
  
% Reserve space for estimates.
MM = zeros(size(M,1),size(Y,2));
PP = zeros(size(M,1),size(M,1),size(Y,2));

% Estimate with EKF
for k=1:size(Y,2)
   [M,P] = ekf_predict1(M,P,A,Q);
   [M,P] = ekf_update1(M,P,Y(:,k),dh_dx_func,R*eye(1),h_func,[],[]);
   MM(:,k)   = M;
   PP(:,:,k) = P;
end


clf; clc;
disp('The filtering results using the EKF is now displayed')

% Project the estimates to measurement space 
Y_m = feval(h_func, MM);
plot(X(1,:),Y,'.', X(1,:), Y_real,'--',MM(1,:),Y_m);
legend('Measurements','Real signal', 'Filtered estimate');
xlim([0 d]);
title('Estimating a random Sine signal with extended Kalman filter.');
print -dps demo2_f1.ps 

clc;
disp('The filtering result using the EKF is now displayed.');
fprintf('Smoothing the estimates using the RTS smoother...');

[SM1,SP1] = erts_smooth1(MM,PP,A,Q);

[SM2,SP2] = efbf_smooth1(MM,PP,Y,A,Q,[],[],[],...
                         dh_dx_func,R*eye(1),h_func,[],[]);  

fprintf('ready.\n');
disp(' ');
disp('Push any button to display the smoothing results.');
pause

Y_s1 = feval(h_func, SM1);
plot(X(1,:),Y,'.', X(1,:), Y_real,'--',SM1(1,:),Y_s1);
legend('Measurements','Real signal','Smoothed estimate');
xlim([0 d]);
title(['Smoothing a random Sine signal with extended',...
       'Kalman (RTS) smoother.']);
print -dps demo2_f2.ps

clc;
disp('The smoothing results using the ERTS smoother is now displayed.');
disp(' ');
disp('Push any button to see the smoothing results of a EFBF smoother.');

pause
  
Y_s2 = feval(h_func, SM2);
plot(X(1,:),Y,'.', X(1,:), Y_real,'--',SM1(1,:),Y_s2);
legend('Measurements','Real signal','Smoothed estimate');
title(['Smoothing a random Sine signal with extended',...
       'Kalman (Forward-Backward) smoother.']);
xlim([0 d]);
print -dps demo2_f3.ps

clc;
disp('The smoothing results using the EFBF smoother is now displayed.');

% Errors.  
EMM_Y = sum((Y_m-Y_real).^2)/n;
ESM1_Y = sum((Y_s1-Y_real).^2)/n;
ESM2_Y = sum((Y_s2-Y_real).^2)/n;
  
disp(' ');
fprintf('Filtering now with UKF...');

% In the rest the signal is filtered with UKF for comparison.

% Initial guesses for the state mean and covariance.
M = [f w a]';
P = diag([3 3 3]);    
  
% Reserve space for estimates.
MM2 = zeros(size(M,1),size(Y,2));
PP2 = zeros(size(M,1),size(M,1),size(Y,2));

% Estimate with UKF
for k=1:size(Y,2)
   [M,P,X_s,w] = ukf_predict3(M,P,f_func,Q,R*eye(1),dt);
   [M,P] = ukf_update3(M,P,Y(:,k),h_func,R*eye(1),X_s,w,[]);
   MM2(:,k)   = M;
   PP2(:,:,k) = P;
end
[SM3, SP3] = urts_smooth1(MM2,PP2,f_func,Q,dt);

fprintf('ready.\n')
disp(' ');
disp('Push any button to see the filtering results.');

pause

Y_m2 = feval(h_func, MM2);
plot(ti,Y,'.', ti, Y_real,'--',ti,Y_m2);
legend('Measurements','Real signal', 'Filtered estimate');
xlim([0 d]);
title('Estimating a random Sine signal with unscented Kalman filter.');

clc;
disp('The filtering results of a UKF is now displayed.')
disp(' ');
disp('Push any button to display the smoothing results.');

pause
  
Y_m3 = feval(h_func, SM3);
plot(X(1,:),Y,'.', X(1,:), Y_real,'--',SM3(1,:),Y_m3);
legend('Measurements','Real signal', 'Filtered estimate');
xlim([0 d]);
title('Estimating a random Sine signal with unscented Kalman smoother (RTS).');

clc;
disp('The smoothing results of a ERTS smoother is now displayed');

UKF_EMM_Y = sum((Y_m2-Y_real).^2)/n;
URTS_EMM_Y = sum((Y_m3-Y_real).^2)/n;

disp(' ');
disp('Mean square errors of all estimates:');
fprintf('EKF-MSE = %.4f\n',sqrt(EMM_Y));
fprintf('ERTS-MSE = %.4f\n',sqrt(ESM1_Y));
fprintf('EFBF-MSE = %.4f\n',sqrt(ESM2_Y));
fprintf('UKF-RMSE = %.4f\n',sqrt(UKF_EMM_Y));
fprintf('URTS-RMSE = %.4f\n',sqrt(URTS_EMM_Y));

