% Demonstration for Kalman filter and smoother using a 2D CWPA model
%
% Copyright (C) 2007 Jouni Hartikainen
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function kf_cwpa_demo (varargin)
doprint = false; % whether we print figures of results
if nargin >= 1
  doprint = true;
end

% Transition matrix for the continous-time system.
F = [0 0 1 0 0 0;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 0 0 0 0 0;
     0 0 0 0 0 0];
nF = size (F, 1);

% Noise effect matrix for the continous-time system.
L = [0 0;
     0 0;
     0 0;
     0 0;
     1 0;
     0 1];

% Stepsize
dt = 0.5;

% Process noise variance
q  = 0.2;
Qc = diag ([q q]);

% Discretization of the continous-time system.
[A,Q] = lti_disc(F,L,Qc,dt);

% Measurement model.
H  = [1 0 0 0 0 0;
      0 1 0 0 0 0];
nH = size (H, 1);

% Variance in the measurements.
r1 = 10;
r2 = 5;
R  = diag ([r1 r1]);

% Generate the data.
n        = 50;
X_r      = zeros (nF, n);
X_r(:,1) = [0 0 0 0 0 0]';
tmp      = zeros (nF, 1);
for i = 2:n
   X_r(:,i) = A * X_r(:,i-1) + gauss_rnd (tmp, Q);
end

% Generate the measurements.
Y   = zeros (nH, n);
tmp = zeros (nH, 1);
for i = 1:n
   Y(:,i) = H * X_r(:,i) + gauss_rnd (tmp, R);
end

clf; clc;
msg_txt = ['This is a demonstration program for the classical Kalman filter.\n', ...
           'KF is used here to estimate the position of a moving object, whos\n',...
           'dynamics follow the CWPA-model described in the documentation\n',...
           'provided with the toolbox.\n', ...
           'We get noisy measurements from the objects position and velocity\n',...
           'on discrete time steps. The real position of the object and the\n',...
           'measurements made of them are displayed now. The blue line is the\n',...
           'real path of the object and the green dots represents the\n',...
           'measurements. The red circle represents the starting point\n',...
           'of the object.\n\n'];
fprintf (msg_txt);
fprintf('\nFiltering with KF... ');

plot (X_r(1,:),X_r(2,:),Y(1,:),Y(2,:),'.',X_r(1,1),...
     X_r(2,1),'ro','MarkerSize',12);
legend ('Real trajectory', 'Measurements');
title ('Position');
axis tight

if doprint% save an image
 print -dpsc demo1_f1.ps
end

% Initial guesses for the state mean and covariance.
m = [0 0 0 0 0 0]';
P = diag ([0.1 0.1 0.1 0.1 0.5 0.5]);

%% Space for the estimates.
MM = zeros (size(m,1), size(Y,2));
PP = zeros (size(m,1), size(m,1), size(Y,2));

% Filtering steps.
for i = 1:size(Y,2)
   [m,P]     = kf_predict (m, P, A, Q);
   [m,P]     = kf_update (m, P, Y(:,i), H, R);
   MM(:,i)   = m;
   PP(:,:,i) = P;
end

% Smoothing step.
[SM,SP]   = rts_smooth (MM,PP,A,Q);
[SM2,SP2] = tf_smooth  (MM,PP,Y,A,Q,H,R,1);

fprintf('ready.\n\n');
input ('<push any button to see the results>');

subplot (1,2,1);
plot (X_r(1,:), X_r(2,:),'--', MM(1,:), MM(2,:),X_r(1,1),X_r(2,1),...
     'o','MarkerSize',12)
legend ('Real trajectory', 'Filtered');
title ('Position estimation with Kalman filter.');
xlabel ('x');
ylabel ('y');

subplot (1,2,2);
plot (X_r(3,:), X_r(4,:),'--', MM(3,:), MM(4,:),X_r(3,1),...
     X_r(4,1),'ro','MarkerSize',12);
legend ('Real velocity', 'Filtered');
title ('Velocity estimation with Kalman filter.');
xlabel ('x^.');
ylabel ('y^.');

if doprint  % save an image
 print -dpsc demo1_f2.ps
end

clc;
msg_txt = ['The filtering results are displayed now. In the left panel the\n',...
          'estimated path is shown, and, for comparison, in the right panel\n',...
          'the estimated velocity is shown.\n\n'];
fprintf (msg_txt);
input('<push any key to see the smoothing results of a RTS smoother>');

subplot (1,2,1);
plot (X_r(1,:), X_r(2,:),'--', SM(1,:), SM(2,:),X_r(1,1),...
     X_r(2,1),'ro','MarkerSize',12);
legend ('Real trajectory', 'Smoothed');
title ('Position estimation with RTS smoother.');
xlabel ('x');
ylabel ('y');

subplot (1,2,2);
plot (X_r(3,:), X_r(4,:),'--', SM(3,:), SM(4,:),X_r(3,1),...
     X_r(4,1),'ro','MarkerSize',12);
legend ('Real velocity', 'Smoothed');
title ('Velocity estimation with RTS smoother.');
xlabel ('x^.');
ylabel ('y^.');
axis tight

if doprint % save an image
 print -dpsc demo1_f3.ps
end

clc;
msg_txt = ['The smoothing results are displayed now. In the left panel the\n',...
           'estimated path is shown, and, for comparison, in the right panel\n',...
           'the estimated velocity is shown.\n\n'];
fprintf (msg_txt)
input ('<push any key to see the filtering results sequentially>');

subplot (1,2,1);
plot (X_r(1,:), X_r(2,:),'--', SM2(1,:), SM2(2,:),X_r(1,1),...
     X_r(2,1),'ro','MarkerSize',12);
legend ('Real trajectory', 'Smoothed');
title ('Position estimation with TF smoother.');
xlabel ('x');
ylabel ('y');
axis tight

subplot (1,2,2);
plot (X_r(3,:), X_r(4,:),'--', SM2(3,:), SM2(4,:),X_r(3,1),...
     X_r(4,1),'ro','MarkerSize',12);
legend ('Real velocity', 'Smoothed');
title ('Velocity estimation with TF smoother.');
xlabel ('x^.');
ylabel ('y^.');
axis tight

% Track and animate the filtering result.
clc
msg_txt = ['Filtering result is now displayed sequentially on every time step.\n', ...
           'The observations are displayed as green dots and\n',...
           'the real signal as solid blue line. The filtering\n',...
           'result of previous steps is plotted as dashed black\n',...
           'line. The covariance in each step is displayed as a\n',...
           'green ellipse around the mean (black circle) on each\n',...
           'step. The predicted position on next step is displayed\n',...
           'as a red circle.\n\n'];
fprintf (msg_txt);
input ('<push any key to proceed to next step>');

clf
EST = [];
tt  = (0:0.01:1) *2 * pi;
nt  = length (tt);
CC  = [cos(tt);sin(tt)];
for k = 1:size(Y,2)
    M   = MM(:,k);
    P   = PP(:,:,k);
    EST = [EST M];
    M_pred = kf_predict(M,P,A,Q);

    % Confidence ellipse
    cc = repmat (M(1:2),1, nt) + 2 * chol (P(1:2,1:2))' * CC;

    % Animate
    if k == 1
      hndl = plot(X_r(1,:),X_r(2,:),'-',...
             Y(1,:), Y(2,:), '.',...
             M(1),M(2),'ko',...
             M_pred(1),M_pred(2),'ro',...
             EST(1,:),EST(2,:),'k--',...
             cc(1,:),cc(2,:),'g-');
      drawnow;
      axis tight
    else
      set (hndl(3), 'xdata', M(1), 'ydata', M(2));
      set (hndl(4), 'xdata', M_pred(1), 'ydata', M_pred(2));
      set (hndl(5), 'xdata', EST(1,:), 'ydata', EST(2,:));
      set (hndl(6), 'xdata', cc(1,:), 'ydata', cc(2,:));
    end
    pause;
end

% Smoothing result.
clf; clc;
msg_txt = ['Smoothing (using RTS smoother) result is now displayed\n',...
           'sequentially on every time step in reverse order.\n', ...
           'The observations are displayed as green dots and the\n',...
           'real signal as solid blue line. The smoothing result\n',...
           'of previous steps is plotted as dashed black line.\n',...
           'The covariance in each step is displayed as a green\n',...
           'ellipse around the mean (black circle) on each step.\n\n']
fprintf (msg_txt)
input ('<push any key to proceed to next step>');

EST = [];
for k = size(Y,2):-1:1
    M   = SM(:,k);
    P   = SP(:,:,k);
    EST = [EST M];

    % Confidence ellipse
    cc = repmat (M(1:2), 1, nt) + 2 * chol (P(1:2,1:2))' * CC;

    % Animate
    if k == size(Y,2)
      hndl = plot(X_r(1,:),X_r(2,:),'-',...
             Y(1,:), Y(2,:), '.',...
             M(1),M(2),'ko',...
             EST(1,:),EST(2,:),'k--',...
             cc(1,:),cc(2,:),'g-');
      drawnow;
      axis tight
    else
      set (hndl(3), 'xdata', M(1), 'ydata', M(2));
      set (hndl(4), 'xdata', EST(1,:), 'ydata', EST(2,:));
      set (hndl(5), 'xdata', cc(1,:), 'ydata', cc(2,:));
    end
    pause;
end

% MSEs of position estimates
MSE_KF1 = mean ((X_r(1,:) - MM(1,:)).^2);
MSE_KF2 = mean ((X_r(2,:) - MM(2,:)).^2);
MSE_KF  = 1/2 * (MSE_KF1 + MSE_KF2);

MSE_RTS1 = mean ((X_r(1,:) - SM(1,:)).^2);
MSE_RTS2 = mean ((X_r(2,:) - SM(2,:)).^2);
MSE_RTS  = 1/2 * (MSE_RTS1 + MSE_RTS2);

MSE_TF1 = mean ((X_r(1,:) - SM2(1,:)).^2);
MSE_TF2 = mean ((X_r(2,:) - SM2(2,:)).^2);
MSE_TF  = 1/2 * (MSE_TF1 + MSE_TF2);

clc;
fprintf('Mean square errors of position estimates:\n');
fprintf('KF-RMSE = %.4f\n',MSE_KF);
fprintf('RTS-RMSE = %.4f\n',MSE_RTS);
fprintf('TF-RMSE = %.4f\n\n',MSE_TF);
