
runs = 1000;
E_KF1 = zeros(1,runs);
E_KS1 = zeros(1,runs);
E_KF2 = zeros(1,runs);
E_KS2 = zeros(1,runs);
E_IMM = zeros(1,runs);
E_IMMS = zeros(1,runs);

for run = 1:runs
    run
    imm_demo;
    E_KF1(run) = MSE_KF1;
    E_KS1(run) = MSE_KS1;
    E_KF2(run) = MSE_KF2;
    E_KS2(run) = MSE_KS2;
    E_IMM(run) = MSE_IMM;
    E_IMMS(run) = MSE_IMMS;
end

mean(E_KF1)
mean(E_KS1)
mean(E_KF2)
mean(E_KS2)
mean(E_IMM)
mean(E_IMMS)
