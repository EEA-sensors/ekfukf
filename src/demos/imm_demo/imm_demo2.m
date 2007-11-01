% Simple demonstration for IMM using velocity and acceleration models

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
dt = 0.5;

ind{1} = [1 2 3 4 5 6]';

% Transition matrix for the continous-time acceleration model.
F1 = [0 0 1 0 0 0;
      0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1;
      0 0 0 0 0 0;
      0 0 0 0 0 0];

% Noise effect matrix for the continous-time system.
L1 = [0 0;
      0 0;
      0 0;
      0 0;
      1 0;
      0 1];

% Process noise variance
q1 = 1.5;
Qc1 = diag([q1 q1]);

% Discretization of the continous-time system.
[A{1},Q{1}] = lti_disc(F1,L1,Qc1,dt);

ind{2} = [1 2 3 4]';

% Transition matrix for the continous-time velocity model.
F2 = [0 0 1 0;
      0 0 0 1;
      0 0 0 0;
      0 0 0 0];
% Noise effect matrix for the continous-time system.
L2 = [0 0;
      0 0;
      1 0;
      0 1];

% Process noise variance
q2 = 0.5;
Qc2 = diag([q2 q2]);

% Discretization of the continous-time system.
[A{2},Q{2}] = lti_disc(F2,L2,Qc2,dt);


% Measurement models.
H{1} = [1 0 0 0 0 0;
        0 1 0 0 0 0];

H{2} = [1 0 0 0;
        0 1 0 0];
hdims = 2;

% Variance in the measurements.
r1 = 10.1;
r2 = 10.1;
R{1} = diag([r1 r2]);

r1 = 10.1;
r2 = 10.1;
R{2} = diag([r1 r2]);


% Generate the data.
% $$$ n = 300;
% $$$ Y = zeros(hdims,n);
% $$$ X_r = zeros(dims,n);
% $$$ X_r(:,1) = [0 0 0 0 0 0]';
% $$$ mstate = zeros(1,n);

mu_ip = [0.1 0.9];
mu_0j = [0.1 0.9];
p_ij = [0.75 0.25;
        0.5 0.5];


% $$$ for i = 2:n
% $$$     if rand < mu_ip(1)
% $$$         X_r(:,i) = A{1}*X_r(:,i-1) + gauss_rnd(zeros(d,1), Q{1});
% $$$     else 
% $$$         X_r(:,i) = A{2}*X_r(:,i-1) + gauss_rnd(zeros(d,1), Q{2});
% $$$     end    
% $$$ end

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


% Acceleration model
% $$$ for i = 2:n
% $$$     st = mstate(i);
% $$$    X_r(ind{st},i) = A{st}*X_r(ind{st},i-1) + gauss_rnd(zeros(size(A{st},1),1), Q{st});
% $$$ end


% $$$ % Generate the measurements.
% $$$ for i = 1:n
% $$$     Y(:,i) = H{mstate(i)}*X_r(ind{mstate(i)},i) + gauss_rnd(zeros(size(Y,1),1), R{mstate(i)});
% $$$ end

% $$$ ac = floor(n/2)+1:floor(n/2)+2;
% $$$ clf; clc;
% $$$ %fprintf('Filtering with KF...');
% $$$ 
% $$$ plot(X_r(1,:),X_r(2,:),Y(1,:),Y(2,:),'.',X_r(1,1),...
% $$$      X_r(2,1),'ro','MarkerSize',12);
% $$$ legend('Real trajectory', 'Measurements');
% $$$ title('Position');

m = [0 0 0 0 0 0]';
P = diag([0.1 0.1 0.1 0.1 0.5 0.5]);

%% Space for the estimates.
% $$$ MM = zeros(size(m,1), size(Y,2));
% $$$ PP = zeros(size(m,1), size(m,1), size(Y,2));
% $$$ MU = zeros(2,size(Y,2));
MM1 = [];
PP1 = [];
MM2 = [];
PP2 = [];
MU = [];
Y = [];
SM = [];
SP = [];

MM1_i = {};
PP1_i = {};
MM2_i = {};
PP2_i = {};


X_r = [];

x_ip{1} = [0 0 0 0 0 0]';
x_ip{2} = [0 0 0 0]';
P_ip{1} = diag([0.1 0.1 0.1 0.1 0.5 0.5]);
P_ip{2} = diag([0.1 0.1 0.1 0.1]);

clf
% Filtering steps.
%for i = 1:size(Y,2)
i = 0;
while(1)
    [x,y,b] = ginput(1);
    if b == 3 && ~isempty(Y)
        MM2_i = cat(2,MM2_i,MM1_i);
        PP2_i = cat(2,PP2_i,PP1_i);
        MM2 = cat(2,MM2,MM1);
        PP2 = cat(3,PP2,PP1);
        MM1 = [];
        PP1 = [];
        [SM,SP,SM_i,SP_i,MU_S] = imm_smooth(MM2,PP2,MM2_i,PP2_i,MU,p_ij,mu_0j,ind,dims,A,Q,R,H,Y);
        plot(SM(1,:),SM(2,:),'r-',X_r(1,:),X_r(2,:),'g-',...
             Y(1,:),Y(2,:),'.');
        xlim([-100 100]) 
        ylim([-100 100]) 
        drawnow
        i = 0;
    else        
        i = i + 1;
    
        % Generate the measurement.
        Y = [Y [[x;y] + gauss_rnd(zeros(2,1), R{1})]];
        X_r = [X_r [x;y]];
    
        [m,P,x_ip,P_ip,mu_ip] = imm_kf_filter(x_ip,P_ip,mu_ip,p_ij,ind,dims,A,Q,R,H,Y(:,end));
        MM1_i(:,i) = x_ip';
        PP1_i(:,i) = P_ip';
        MM1 = [MM1 m];
        PP1(:,:,i) = P;
        MU = [MU mu_ip'];

        plot(MM1(1,:),MM1(2,:),'k-',X_r(1,:),X_r(2,:),'g-',...
             Y(1,:),Y(2,:),'.');
        if ~isempty(SM)
            hold on
            plot(SM(1,:),SM(2,:),'r-');
            hold off
        end
        xlim([-100 100]) 
        ylim([-100 100]) 
        drawnow
    end
end


%plot(MU');
%pause

st1 = find(mstate == 1);
st2 = find(mstate == 2);
plot(MM(1,:),MM(2,:),'k-',X_r(1,:),X_r(2,:),'g-',...
     Y(1,:),Y(2,:),'.');
