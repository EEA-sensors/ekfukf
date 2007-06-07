% Demonstration of univariate nonstationary growth model (UNGM) using
% EKF and UKF.

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


f_func = @ungm_f;
h_func = @ungm_h;
df_func = @ungm_df_dx;
dh_func = @ungm_dh_dx;

n = 500;
x_0 = .1;

Y = zeros(1,n);
Y = zeros(1,n);

u_n = 1;
v_n = 1;

% Generate the true states with process noise
X = zeros(1,n);
X(1) = ungm_f(x_0,1) + gauss_rnd(0,u_n,1);
for i = 2:n
    X(i) = feval(f_func,X(i-1),i) + gauss_rnd(0,u_n,1);
end
    
% Generate the observations with measurement noise
for i = 1:n
    Y(i) = feval(h_func,X(i)) + gauss_rnd(0,v_n,1);
end

Y_r = feval(h_func,X(i));

M = 0;
P = 1;

MM = zeros(size(M,1),size(Y,2));
PP = zeros(size(M,1),size(M,1),size(Y,2));

% Filtering loop for UKF
for k = 1:size(Y,2)
    [M,P,X_s,w] = ukf_predict3(M,P,f_func,u_n,v_n,k);
    [M,P] = ukf_update3(M,P,Y(:,k),h_func,v_n,X_s,w,[]);
    %[M,P] = ukf_predict1(M,P,f_func,u_n,k);
    %[M,P] = ukf_update1(M,P,Y(:,k),h_func,v_n);
   MM(:,k)   = M;
   PP(:,:,k) = P;    
end

params = cell(size(Y,2));
for i = 1:size(Y,2)
   params{i} = i+1; 
end

[MMS1, PPS1] = urts_smooth1(MM,PP,f_func,u_n,params,[],[],[],[],0);

%[MMS2, PPS2] = ufbf_smooth1(MM,PP,Y,[],u_n,[],h_func,v_n,[],[],[],[],0);

% Visualize the estimates
Y_f = zeros(1,n);
for k = 1:size(Y,2)
   Y_f(k) = feval(h_func,MM(:,k),[]);
end

plot(1:100,Y_r(1:100),'-kx',1:100,Y_f(1:100),'--bo')
title('UKF filtering result projected to observation space');
legend('Observations','UKF filtered estimate'); 
xlim([0 100]);
ylim([-5 40]);
print -dpsc ungm_f1.ps

pause

plot(1:100,X(1:100),'-kx',1:100,MM(1:100),'--bo')
title('UKF state filtering result');
xlim([0 100]);
ylim([-40 40]);
legend('Real signal', 'UKF filtered estimate');
print -dpsc ungm_f2.ps

pause

plot(1:100,X(1:100),'-kx',1:100,MMS1(1:100),'--bo')
title('URTS state smoothing result');
xlim([0 100]);
ylim([-40 40]);
legend('Real signal', 'ERTS smoothed estimate');

pause

E = X-MM;
PP1 = squeeze(PP);
plot(1:100,E(1:100),'-kx',1:100,3*sqrt(PP1(1:100))','--r',1:100,-3*sqrt(PP1(1:100))','--r');
title('UKF error');
legend('Prediction error of UKF','3\sigma interval')

pause

UKF_MSE = sum((X-MM).^2)/n;
fprintf('UKF-MSE = %.4f\n',UKF_MSE);
UKS1_MSE = sum((X-MMS1).^2)/n;
fprintf('UKS1-MSE = %.4f\n',UKS1_MSE);


M = x_0;
P = 1;

MM2 = zeros(size(M,1),size(Y,2));
PP2 = zeros(size(M,1),size(M,1),size(Y,2));

% Filtering loop for EKF
for k = 1:size(Y,2)
   [M,P] = ekf_predict1(M,P,df_func,u_n,[],[],k);
   [M,P] = ekf_update1(M,P,Y(:,k),dh_func,v_n,h_func,[],[]);
   MM2(:,k)   = M;
   PP2(:,:,k) = P;    
end

[MMS3,PPS3] = erts_smooth1(MM2,PP2,df_func,u_n,f_func,[],params,0);

% Visualize the estimates
Y_f2 = zeros(1,n);
for k = 1:size(Y,2)
   Y_f2(k) = feval(h_func,MM2(:,k),[]);
end

plot(1:100,Y_r(1:100),'-kx',1:100,Y_f2(1:100),'--bo')
title('EKF filtering result of observations');
legend('Observations','EKF filtered estimate'); 

pause

plot(1:100,X(1:100),'-kx',1:100,MM2(1:100),'--bo')
title('EKF filtering result of states');
legend('Real signal', 'EKF filtered estimate');

pause

plot(1:100,X(1:100),'-kx',1:100,MMS3(1:100),'--bo')
title('ERTS smoothed states');

pause

E = X-MM2;
PP2 = squeeze(PP2);
plot(1:100,E(1:100),'-kx',1:100,3*sqrt(PP2(1:100))','-r',1:100,-3*sqrt(PP2(1:100))','-r');
title('EKF error');
legend('Prediction error of EKF','3\sigma interval');


EKF_MSE = sum((X-MM2).^2)/n;
fprintf('EKF-MSE = %.4f\n',EKF_MSE);
ERTS_MSE = sum((X-MMS3).^2)/n;
fprintf('ERTS-MSE = %.4f\n',ERTS_MSE);


M = x_0;
P = 1;
SX = gauss_rnd(M,P,10000);

MM3 = zeros(size(M,1),size(Y,2));
PP3 = zeros(size(M,1),size(M,1),size(Y,2));

% Filtering loop for bootstrap filter
for k = 1:size(Y,2)
   SX = ungm_f(SX,k) + gauss_rnd(0,1,size(SX,2));
   W  = gauss_pdf(Y(:,k),ungm_h(SX,k),1);
   ind = resampstr(W);
   SX = SX(:,ind);
   M = mean(SX);
   P = var(SX);
   MM3(:,k)   = M;
   PP3(:,:,k) = P;    
end


% Visualize the estimates
Y_f3 = zeros(1,n);
for k = 1:size(Y,2)
   Y_f3(k) = feval(h_func,MM3(:,k),[]);
end

plot(1:100,Y_r(1:100),'-kx',1:100,Y_f3(1:100),'--bo')
title('Bootstrap filtering result of observations')
legend('Original observations', 'filtered estimate')

pause

plot(1:100,X(1:100),'-kx',1:100,MM3(1:100),'--bo')
title('Bootstrap filtering result of states')
legend('Real states','Filtered estimates')

pause

E = X-MM3;
PP3 = squeeze(PP3);
plot(1:100,E(1:100),'-kx',1:100,3*sqrt(PP3(1:100))','-r',1:100,-3*sqrt(PP3(1:100))','-r');
title('Bootstrap filtering error');
legend('Prediction error','3\sigma interval')

BS_MSE = sum((X-MM3).^2)/n;
fprintf('BS-MSE = %.4f\n',BS_MSE);

