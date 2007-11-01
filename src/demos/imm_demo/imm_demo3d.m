% Simple demonstration for IMM using velocity and acceleration models

% Dimensionality of the state space
dims = 9;

nmodels = 2;

% Space for models
ind = cell(1,2);
A = cell(1,2);
Q = cell(1,2);
H = cell(1,2);
R = cell(1,2);

% Stepsize
dt = 0.1;

ind{1} = [1 2 3 4 5 6 7 8 9]';

% Transition matrix for the continous-time acceleration model.
F1 = [0 0 0 1 0 0 0 0 0;
      0 0 0 0 1 0 0 0 0;
      0 0 0 0 0 1 0 0 0;
      0 0 0 0 0 0 1 0 0;
      0 0 0 0 0 0 0 1 0;
      0 0 0 0 0 0 0 0 1;
      0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0];

% Noise effect matrix for the continous-time system.
L1 = [0 0 0;
      0 0 0;
      0 0 0;
      0 0 0;
      0 0 0;
      0 0 0;
      1 0 0;
      0 1 0;
      0 0 1];

% Process noise variance
q1 = 2;
Qc1 = diag([q1 q1 q1/5]);

% Discretization of the continous-time system.
[A{1},Q{1}] = lti_disc(F1,L1,Qc1,dt);

ind{2} = [1 2 3 4 5 6]';

% Transition matrix for the continous-time velocity model.
F2 = [0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0];
% Noise effect matrix for the continous-time system.
L2 = [0 0 0;
      0 0 0;
      0 0 0;
      1 0 0;
      0 1 0;
      0 0 1];

% Process noise variance
q2 = 0.01;
Qc2 = diag([q2 q2 q2/10]);

% Discretization of the continous-time system.
[A{2},Q{2}] = lti_disc(F2,L2,Qc2,dt);


% Measurement models.
H{1} = [1 0 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0 0];

H{2} = [1 0 0 0 0 0;
        0 1 0 0 0 0;
        0 0 1 0 0 0];
hdims = 3;

% Variance in the measurements.
r1 = .10;
r2 = .10;
r3 = .10;
R{1} = diag([r1 r2 r3]);

r1 = .01;
r2 = .01;
r3 = .01;
R{2} = diag([r1 r2 r3]);


% Generate the data.
n = 1000;
Y = zeros(hdims,n);
X_r = zeros(dims,n);
X_r(:,1) = [0 0 0 0 0 0 0 0 0]';
mstate = zeros(1,n);

mu_ip = [0.1 0.9];
p_ij = [0.65 0.35;
        0.10 0.90];


% $$$ for i = 2:n
% $$$     if rand < mu_ip(1)
% $$$         X_r(:,i) = A{1}*X_r(:,i-1) + gauss_rnd(zeros(d,1), Q{1});
% $$$     else 
% $$$         X_r(:,i) = A{2}*X_r(:,i-1) + gauss_rnd(zeros(d,1), Q{2});
% $$$     end    
% $$$ end

mstate(1) = 2;
for i = 2:n
    r = rand;
    for j = size(p_ij,1)-1;
        if r < p_ij(mstate(i-1),j);
            mstate(i) = j;
        end
    end
    if mstate(i) == 0
        mstate(i) = size(p_ij,1);
    end
end


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

plot(X_r(1,:),X_r(2,:),Y(1,:),Y(2,:),'.',X_r(1,1),...
     X_r(2,1),'ro','MarkerSize',12);
legend('Real trajectory', 'Measurements');
title('Position');

m = [0 0 0 0 0 0 0 0 0]';
P = diag([0.1 0.1 0.1 0.1 0.1 0.1 0.5 0.5 0.5]);

%% Space for the estimates.
MM = zeros(size(m,1), size(Y,2));
PP = zeros(size(m,1), size(m,1), size(Y,2));
MU = zeros(2,size(Y,2));
MM_i = cell(2,n);
PP_i = cell(2,n);

x_ip{1} = [0 0 0 0 0 0 0 0 0]';
x_ip{2} = [0 0 0 0 0 0]';
P_ip{1} = diag([0.1 0.1 0.1 0.1 0.1 0.1 0.5 0.5 0.5]);
P_ip{2} = diag([0.1 0.1 0.1 0.1 0.1 0.1]);


% Filtering steps.
for i = 1:size(Y,2)
   [m,P,x_ip,P_ip,mu_ip] = imm_kf_filter(x_ip,P_ip,mu_ip,p_ij,ind,dims,A,Q,R,H,Y(:,i));
   MM(:,i)   = m;
   PP(:,:,i) = P;
   MU(:,i)   = mu_ip';
   MM_i(:,i) = x_ip';
   PP_i(:,i) = P_ip';
   plot3(MM(1,1:i),MM(2,1:i),MM(3,1:i),'k-',X_r(1,1:i),X_r(2,1:i),X_r(3,1:i),'g-',...
        Y(1,1:i),Y(2,1:i),Y(3,1:i),'.');
   %xlim([-100 100]) 
   %ylim([-100 100]) 
   drawnow
end

%[SM1,SP1,SS1] = rts_smooth(MM,PP,A{1},Q{1});
%[SM2,SP2,SS2] = rts_smooth(MM(ind{2},ind{2}),PP(ind{2},ind{2},:),A{2},Q{2});

[SM3,SP3,SM3_i,SP3_i,MU_S] = imm_smooth(MM,PP,MM_i,PP_i,MU,p_ij,mu_0j,ind,dims,A,Q,R,H,Y);

%plot(MU');
%pause

% $$$ st1 = find(mstate == 1);
% $$$ st2 = find(mstate == 2);
% $$$ plot(MM(1,:),MM(2,:),'k-',X_r(1,:),X_r(2,:),'g-',...
% $$$      Y(1,:),Y(2,:),'.');
% $$$ 
% $$$ plot(SM1(1,:),SM1(2,:),'k-',X_r(1,:),X_r(2,:),'g-',...
% $$$      Y(1,:),Y(2,:),'.',SM2(1,:),SM2(2,:),'y-');

plot3(SM3(1,:),SM3(2,:),SM3(3,:),'k-',X_r(1,:),X_r(2,:),X_r(3,:),'g-',...
      Y(1,:),Y(2,:),Y(3,:),'.',MM(1,:),MM(2,:),MM(3,:),'r-');
   