% Demonstration for EKF using a random sine signal model. 
%
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

clc;
disp('Filtering the signal with EKF...');

save_plots = 1;

% Measurement model and it's derivative
f_func = @ekf_sine_f;
h_func = @ekf_sine_h;
dh_dx_func = @ekf_sine_dh_dx;
d2h_dx2_func = @ekf_sine_d2h_dx2;

% Initial values for the signal.
f = 0;
w = 10;
a = 1;
  
% Number of samples and stepsize.
d = 5;
n = 500;
dt = d/n;
x = 1:n;

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
  
% Generate the observations with Gaussian noise.
sd = 1;
R = sd^2;

Y = zeros(1,n);
Y_real = feval(h_func,X);     
Y = Y_real + gauss_rnd(0,R,n);
  
plot(x,Y,'.',x,Y_real)
  
% Initial guesses for the state mean and covariance.
M = [f w a]';
P = diag([3 3 3]);    
  
% Reserve space for estimates.
MM = zeros(size(M,1),size(Y,2));
PP = zeros(size(M,1),size(M,1),size(Y,2));

% Estimate with EKF
for k=1:size(Y,2)
   [M,P] = ekf_predict1(M,P,A,Q);
   [M,P] = ekf_update1(M,P,Y(:,k),dh_dx_func,R*eye(1),h_func);
   MM(:,k)   = M;
   PP(:,:,k) = P;
end

% Initial guesses for the state mean and covariance.
M = [f w a]';
P = diag([3 3 3]);    
  
% Reserve space for estimates.
MM2 = zeros(size(M,1),size(Y,2));
PP2 = zeros(size(M,1),size(M,1),size(Y,2));

% Estimate with EKF
for k=1:size(Y,2)
   [M,P] = ekf_predict1(M,P,A,Q);
   [M,P] = ekf_update2(M,P,Y(:,k),dh_dx_func,d2h_dx2_func,R*eye(1),h_func);
   MM2(:,k)   = M;
   PP2(:,:,k) = P;
end

clf; clc;
disp('The filtering results using the 1st order EKF is now displayed')

% Project the estimates to measurement space 
Y_m = feval(h_func, MM);
Y_m2 = feval(h_func, MM2);
plot(x,Y,'.', x,Y_real,'--',x,Y_m);
legend('Measurements','Real signal', 'Filtered estimate');
xlim([0 ceil(max(x))]);
title('Estimating a random Sine signal with extended Kalman filter.');

if save_plots
    print -dpsc demo2_f1.ps 
end

clc;
disp('The filtering result using the 1st order EKF is now displayed.');
fprintf('Smoothing the estimates using the RTS smoother...');

[SM1,SP1] = erts_smooth1(MM,PP,A,Q);

[SM2,SP2] = etf_smooth1(MM,PP,Y,A,Q,[],[],[],...
                         dh_dx_func,R*eye(1),h_func);  

[SM1_2,SP1_2] = erts_smooth1(MM2,PP2,A,Q);

[SM2_2,SP2_2] = etf_smooth1(MM2,PP2,Y,A,Q,[],[],[],...
                         dh_dx_func,R*eye(1),h_func);  

fprintf('ready.\n');
disp(' ');
disp('Push any button to display the smoothing results.');
pause

Y_s1 = feval(h_func, SM1);
plot(x,Y,'.', x, Y_real,'--',x,Y_s1);
legend('Measurements','Real signal','Smoothed estimate');
xlim([0 ceil(max(x))]);
title(['Smoothing a random Sine signal with extended ',...
       'Kalman (RTS) smoother.']);

if save_plots
    print -dpsc demo2_f2.ps
end

clc;
disp('The smoothing results using the ERTS smoother is now displayed.');
disp(' ');
disp('Push any button to see the smoothing results of a ETF smoother.');
pause
  
Y_s2 = feval(h_func, SM2);
plot(x,Y,'.', x, Y_real,'--',x,Y_s2);
legend('Measurements','Real signal','Smoothed estimate');
title(['Smoothing a random Sine signal with extended ',...
       'Kalman (Two Filter) smoother.']);
xlim([0 ceil(max(x))]);

if save_plots
    print -dpsc demo2_f3.ps
end

clc;
disp('The smoothing results using the ETF smoother is now displayed.');


Y_s1_2 = feval(h_func, SM1_2);
Y_s2_2 = feval(h_func, SM2_2);
% Errors.  
EMM_Y = sum((Y_m-Y_real).^2)/n;
EMM_Y2 = sum((Y_m2-Y_real).^2)/n;
ESM1_Y = sum((Y_s1-Y_real).^2)/n;
ESM2_Y = sum((Y_s2-Y_real).^2)/n;
ESM1_2_Y = sum((Y_s1_2-Y_real).^2)/n;
ESM2_2_Y = sum((Y_s2_2-Y_real).^2)/n;
  
disp(' ');
fprintf('Filtering now with UKF...');

% In the rest the signal is filtered with UKF for comparison.

% Initial guesses for the state mean and covariance.
M = [f w a]';
P = diag([3 3 3]);    
  
% Reserve space for estimates.
U_MM = zeros(size(M,1),size(Y,2));
U_PP = zeros(size(M,1),size(M,1),size(Y,2));

% Estimate with UKF
for k=1:size(Y,2)
   [M,P,X_s,w] = ukf_predict3(M,P,f_func,Q,R*eye(1),dt);
   [M,P] = ukf_update3(M,P,Y(:,k),h_func,R*eye(1),X_s,w,[]);
   U_MM(:,k)   = M;
   U_PP(:,:,k) = P;
end
[U_SM, U_SP] = urts_smooth1(U_MM,U_PP,f_func,Q,dt);

fprintf('ready.\n')
disp(' ');
disp('Push any button to see the filtering results.');
pause

Y_m_u = feval(h_func, U_MM);
plot(x,Y,'.', x, Y_real,'--',x,Y_m_u);
legend('Measurements','Real signal', 'Filtered estimate');
xlim([0 ceil(max(x))]);
title('Estimating a random Sine signal with unscented Kalman filter.');

clc;
disp('The filtering results of a UKF are now displayed.')
disp(' ');
disp('Push any button to display the smoothing results.');
pause
  
Y_m_su = feval(h_func, U_SM);
plot(x,Y,'.', x, Y_real,'--',x,Y_m_su);
legend('Measurements','Real signal', 'Filtered estimate');
xlim([0 ceil(max(x))]);
title('Estimating a random Sine signal with unscented Kalman smoother (RTS).');

clc;
disp('The smoothing results of a ERTS smoother are now displayed');

UKF_EMM_Y = sum((Y_m_u-Y_real).^2)/n;
URTS_EMM_Y = sum((Y_m_su-Y_real).^2)/n;

disp(' ');
disp('Mean square errors of all estimates:');
fprintf('EKF1-MSE = %.4f\n',sqrt(EMM_Y));
fprintf('ERTS-MSE = %.4f\n',sqrt(ESM1_Y));
fprintf('ETF-MSE = %.4f\n',sqrt(ESM2_Y));
fprintf('EKF2-MSE = %.4f\n',sqrt(EMM_Y2));
fprintf('ERTS2-MSE = %.4f\n',sqrt(ESM1_2_Y));
fprintf('ETF2-MSE = %.4f\n',sqrt(ESM2_2_Y));
fprintf('UKF-RMSE = %.4f\n',sqrt(UKF_EMM_Y));
fprintf('URTS-RMSE = %.4f\n',sqrt(URTS_EMM_Y));

