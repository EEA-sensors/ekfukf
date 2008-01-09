
runs = 100;
E_EKF = zeros(1,runs);
E_EKS = zeros(1,runs);
E_UKF = zeros(1,runs);
E_UKS = zeros(1,runs);
E_IMM1 = zeros(1,runs);
E_IMMS1 = zeros(1,runs);
E_IMM2 = zeros(1,runs);
E_IMMS2 = zeros(1,runs);

for run = 1:runs
    run
    eimm_demo2;
    E_EKF(run) = MSE_EKF1;
    E_EKS(run) = MSE_EKS1;
    E_UKF(run) = MSE_UKF1;
    E_UKS(run) = MSE_UKS1;    
    E_IMM1(run) = MSE_EIMM;
    E_IMMS1(run) = MSE_EIMMS;
    E_IMM2(run) = MSE_UIMM;
    E_IMMS2(run) = MSE_UIMMS;

end

mean(E_EKF)
mean(E_EKS)
mean(E_UKF)
mean(E_UKS)
mean(E_IMM1)
mean(E_IMMS1)
mean(E_IMM2)
mean(E_IMMS2)


