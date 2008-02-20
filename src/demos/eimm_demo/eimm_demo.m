% Simple demonstration for IMM using velocity and acceleration models

print_figures = 1;
save_figures = 1;

% Dimensionality of the state space
dims = 7;

nmodels = 3;

% Step size
dt = 0.1;

a_func = {};
ia_func = {};
a_param = {};
h_func = {};
h_param = {};
dh_dx_func = {};

a_func{1} = [];
a_func{2} = [];
a_func{3} = @f_turn;
ia_func{1} = [];
ia_func{2} = [];
ia_func{3} = @f_turn_inv;
a_param{1} = [];
a_param{2} = [];
a_param{3} = {dt};
h_func{1} = @bot_h;
h_func{2} = @bot_h;   
h_func{3} = @bot_h;
dh_dx_func{1} = @bot_dh_dx;
dh_dx_func{2} = @bot_dh_dx;
dh_dx_func{3} = @bot_dh_dx;



% Space for models
ind = cell(1,nmodels);
F   = cell(1,nmodels);
L   = cell(1,nmodels);
Qc  = cell(1,nmodels);
A   = cell(1,nmodels);
Q   = cell(1,nmodels);
H   = cell(1,nmodels);
R   = cell(1,nmodels);

ind{2} = [1 2 3 4 5 6]';

% Transition matrix for the continous-time acceleration model.
F{2} = [0 0 1 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        0 0 0 0 0 0;
        0 0 0 0 0 0];

% Noise effect matrix for the continous-time system.
L{2} = [0 0;
        0 0;
        0 0;
        0 0;
        1 0;
        0 1];

% Process noise variance
q2 = .1;
Qc{2} = diag([q2 q2]);

% Discretization of the continous-time system.
[A{2},Q{2}] = lti_disc(F{2},L{2},Qc{2},dt);

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
q1 = .01;
Qc{1} = diag([q1 q1]);

% Discretization of the continous-time system.
[A{1},Q{1}] = lti_disc(F{1},L{1},Qc{1},dt);

% Specification of the turning model

% System components. 7th parameter is the turning rate 
ind{3} = [1 2 3 4 7]';
% Dynamic function 
A{3} = @f_turn_dx;
% Process noise for the turning rate
Qc{3} = 0.001;

% Noise effect matrix
L{3} = [0 0 0 0 1]';

% Process noise covariance
Q{3} = L{3}*Qc{3}*L{3}'*dt;

% Check the derivatives
der_check(a_func{3}, A{3}, 1, [1 1 1 1 0.00001]',{dt});

%mu_ip = [0.90 0.05 0.05];
mu_ip = [0.8 0.1 0.1];
mu_0j = mu_ip;
%p_ij = [0.65 0.35;
%        0.10 0.90];

p_ij = [0.96 0.02 0.02;
        0.02 0.96 0.02;
        0.02 0.02 0.96];

% Generate the data.
n = 200;
X_r = zeros(dims,n);
X_r(:,1) = [0 0 0 -1 0 0 1]';
X_r(7,70) = randn*0.1+0.5;
mstate = zeros(1,n);

% Generate state transitions
mstate(1) = 1;
for i = 2:n
    r = rand;
    for j = 1:size(p_ij,2);
        if r < p_ij(mstate(i-1),j);
            mstate(i) = j;
            break
        else 
            r = r - p_ij(mstate(i-1),j);
        end        
    end
    if mstate(i) == 0
        mstate(i) = size(p_ij,1);
    end
end

% $$$ % Force model transitions
% $$$ mstate(1:50) = 1;
% $$$ mstate(51:70) = 2;
% $$$ mstate(71:120) = 3;
% $$$ mstate(121:150) = 2;
% $$$ mstate(151:170) = 3;
% $$$ mstate(171:200) = 1;

% Generate state values
for i = 2:n
   st = mstate(i);
   if isstr(a_func{st}) | strcmp(class(a_func{st}),'function_handle')
       X_r(ind{st},i) = feval(a_func{st},X_r(ind{st},i-1),a_param{st}) + L{st} * gauss_rnd(zeros(size(Qc{st},2)),Qc{st});
   else 
       X_r(ind{st},i) = A{st}*X_r(ind{st},i-1) + gauss_rnd(zeros(size(A{st},1),1), Q{st});
   end
end

% $$$ S1 = [-10;-10];
% $$$ S2 = [10;10];
% $$$ S3 = [-10;10];
% $$$ S4 = [10;-10];
% $$$ 
% $$$ s = [S1 S2 S3 S4];

x_min = min(X_r(1,:));
x_max = max(X_r(1,:));
y_min = min(X_r(2,:));
y_max = max(X_r(2,:));
offset = 5;

% Positions of sensors
% $$$ S1 = [ 50;-50];
% $$$ S2 = [ 50; 50];
% $$$ S3 = [-50; 50];
% $$$ S4 = [-50;-50];

S1 = [x_max + offset ; y_min - offset];
S2 = [x_max + offset ; y_max + offset];
S3 = [x_min - offset ; y_max + offset];
S4 = [x_min - offset ; y_min - offset];

s = [S1 S2 S3 S4];




% $$$ nsensors = 4;
% $$$ s = zeros(2,nsensors);
% $$$ for i = 1:nsensors
% $$$     s(:,i) = [100*rand-50;100*rand-50];
% $$$ end
% $$$ plot(s(1,:),s(2,:),'.')


h_param{1} = s;
h_param{2} = s;
h_param{3} = s;

% Standard deviation
%sd = 0.05;
%sd = 0.10*rand(1,size(s,2));
sd = 0.01*ones(1,size(s,2));
R = {};
R{1} = diag(sd);
R{2} = R{1};
R{3} = R{1};
% Generate measurements
Y = bot_h(X_r,s);
% Add noise
for i = 1:size(Y,2)
    st = mstate(i);
    for j = 1:size(Y,1)
        Y(j,i) = Y(j,i) + R{st}(j,j) * randn; 
    end
end

ac = floor(n/2)+1:floor(n/2)+2;
clf; clc;
%fprintf('Filtering with KF...');

plot(X_r(1,:),X_r(2,:),...
     X_r(1,1),X_r(2,1),'ro','MarkerSize',12);
legend('Real trajectory', 'Measurements');
title('Position');
pause


m = [0 0 0 -1 0 0 0]';
P = diag([10.1 10.1 1.1 1.1 0.5 0.5 1]);

%% Space for the estimates.

% EKF with model 1
MM1 = zeros(size(A{1},1), size(Y,2));
PP1 = zeros(size(A{1},1), size(A{1},1), size(Y,2));

% EKF with model 2
MM2 = zeros(size(A{2},1), size(Y,2));
PP2 = zeros(size(A{2},1), size(A{2},1), size(Y,2));

% IMM
MM = zeros(size(m,1), size(Y,2));
PP = zeros(size(m,1), size(m,1), size(Y,2));
MM_i = cell(nmodels,n);
PP_i = cell(nmodels,n);
MU = zeros(nmodels,size(Y,2));

%%% Initial estimates %%%

% EKF with model 1
M1 = [0 0 0 -1]';
P1 = diag([1.1 1.1 0.1 0.1]);

% EKF with model 2
M2 = [0 0 0 -1 0 0]';
P2 = diag([1.1 1.1 0.1 0.1 0.5 0.5]);

% IMM
x_ip{1} = [0 0 0 -1]';
x_ip{2} = [0 0 0 -1 0 0]';
x_ip{3} = [0 0 0 -1 1]';

P_ip{1} = diag([1.1 1.1 0.1 0.1]);
P_ip{2} = diag([1.1 1.1 0.1 0.1 0.5 0.5]);
P_ip{3} = diag([1.1 1.1 0.1 0.1 1]);


% Filtering steps.
for i = 1:size(Y,2)
    % EKF with model 1
    [M1,P1] = ekf_predict1(M1,P1,A{1},Q{1});
    [M1,P1] = ekf_update1(M1,P1,Y(:,i),dh_dx_func{1},R{1},h_func{1},[],h_param{1});
    MM1(:,i)   = M1;
    PP1(:,:,i) = P1;
    
    % EKF with model 2
    [M2,P2] = ekf_predict1(M2,P2,A{2},Q{2});
    [M2,P2] = ekf_update1(M2,P2,Y(:,i),dh_dx_func{2},R{2},h_func{2},[],h_param{2});
    MM2(:,i)   = M2;
    PP2(:,:,i) = P2;
        
    %[x_ip,P_ip,mu_ip,m,P] = imm_filter(x_ip,P_ip,mu_ip,p_ij,ind,dims,A,Q,Y(:,i),H,R);
    %[x_p,P_p,c_j] = imm_predict(x_ip,P_ip,mu_ip,p_ij,ind,dims,A,Q);
    [x_p,P_p,c_j] = eimm_predict(x_ip,P_ip,mu_ip,p_ij,ind,dims,A,a_func,a_param,Q);

    [x_ip,P_ip,mu_ip,m,P] = eimm_update(x_p,P_p,c_j,ind,dims,Y(:,i),dh_dx_func,h_func,R,h_param);
    %[x_ip,P_ip,mu_ip,m,P] = uimm_update(x_p,P_p,c_j,ind,dims,Y(:,i),h_func,R,h_param);
    
    MM(:,i)   = m;
    PP(:,:,i) = P;
    MU(:,i)   = mu_ip';
    MM_i(:,i) = x_ip';
    PP_i(:,i) = P_ip';
    plot(MM1(1,1:i),MM1(2,1:i),'y-',...
         MM2(1,1:i),MM2(2,1:i),'b-',...
         MM(1,1:i),MM(2,1:i),'r-',...
         X_r(1,1:i),X_r(2,1:i),'g-');
    %plot(X_r(1,1:i),X_r(2,1:i),'g-');
   
   % Measurement directions

   %x_temp = repmat(X_r(1:2,i),1,size(s,2));
   %len = sqrt(sum((x_temp-s).^2,1));
   
   hold on
   for k = 1:size(s,2)
       len = sqrt(sum((X_r(1:2,i)-s(:,k)).^2,1));
       %len = 20.5;
       dx = len*cos(Y(k,i));
       dy = len*sin(Y(k,i));
   
       plot([s(1,k);s(1,k)+dx], [s(2,k);s(2,k)+dy], 'k--')
       %xlim([-60 60]) 
       %ylim([-60 60]) 
   end
   hold off
   drawnow
   %pause
end

[SM3,SP3,SM3_i,SP3_i,MU_S] = eimm_smooth(MM,PP,MM_i,PP_i,MU,p_ij,mu_0j,ind,dims,A,ia_func,a_param,Q,R,dh_dx_func,h_func,h_param,Y);

[SM1,SP1] = erts_smooth1(MM1,PP1,A{1},Q{1});
[SM2,SP2] = erts_smooth1(MM2,PP2,A{2},Q{2});

%[SM2,SP2] = erts_smooth1(MM(ind{2},:),PP(ind{2},ind{2},:),A{2},Q{2});

% $$$ 
% $$$ plot(X_r(1,:),X_r(2,:),'g-',SM1(1,:),SM1(2,:),'c-',...
% $$$      SM3(1,:),SM3(2,:),'m-',MM(1,:),MM(2,:),'r-');

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
    h = plot(X_r(1,:),X_r(2,:),'g-',...         
             MM1(1,:),MM1(2,:),'r-',...
             SM1(1,:),SM1(2,:),'b-');
    legend('True trajectory',...
           'Filtered',...
           'Smoothed');
    title('Estimates produced by Kalman filter using the model 1.');
    set(h,'markersize',2);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',6);
    if save_figures
        print('-depsc','eimm1.eps');
    end
    pause
    
    h = plot(X_r(1,:),X_r(2,:),'g-',...         
             MM2(1,:),MM2(2,:),'r-',...
             SM2(1,:),SM2(2,:),'b-');
    legend('True trajectory',...
           'Filtered',...
           'Smoothed');
    title('Estimates produced by Kalman filter using the model 2.');
    set(h,'markersize',2);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',6);
    if save_figures
        print('-depsc','eimm2.eps');
    end
    pause
    
    h = plot(X_r(1,:),X_r(2,:),'g-',...         
             MM(1,:),MM(2,:),'r-',...
             SM3(1,:),SM3(2,:),'b-');
    legend('True trajectory',...
           'Filtered',...
           'Smoothed');
    title('Estimates produced by IMM-filter.')
    set(h,'markersize',2);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',6);
    if save_figures
        print('-depsc','eimm3.eps');
    end
    pause
    
    % Determine the real model probabilities
    p_models = zeros(nmodels,n);
    I1 = find(mstate == 1);
    p_models(1,I1) = 1;
    I2 = find(mstate == 2);
    p_models(2,I2) = 1;
    I3 = find(mstate == 3);
    p_models(3,I3) = 1;
 
    
% $$$     h = plot(1:n, p_models(1,:),'r',...
% $$$              1:n,MU(1,:)','r-',...
% $$$              1:n, p_models(2,:),'b',...
% $$$              1:n,MU(2,:)','b-',...
% $$$              1:n, p_models(3,:),'g',...
% $$$              1:n,MU(3,:)','g-')

    % Plot model 1 probability
    h = plot(1:n, p_models(1,:),'g--',...
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
    pause

    % Plot model 2 probability
    h = plot(1:n, p_models(2,:),'g--',...
             1:n,MU(2,:)','r-',...
             1:n,MU_S(2,:)','b-');
    legend('True',...
           'Filtered',...
           'Smoothed');
    title('Probability of model 2');
    ylim([-0.1,1.1]);
    set(h,'markersize',2);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',6);
    pause
    
    % Plot model 3 probability
    h = plot(1:n, p_models(3,:),'g--',...
             1:n,MU(3,:)','r-',...
             1:n,MU_S(3,:)','b-');
    legend('True',...
           'Filtered',...
           'Smoothed');
    title('Probability of model 3');
    ylim([-0.1,1.1]);
    set(h,'markersize',2);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',6);


end