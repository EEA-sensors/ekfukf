

%% Run the UNGM demo ns times to obtain some reliable results

  % Number of loops
  ns = 100;
  
  % Allocate space for results
    
  mse_results = zeros(11,ns);
  
  for loop = 1:ns

      % Run UNGM demo in silent mode
      ungm_demo
      
      mse_results(:,loop) = [UKF1_MSE; UKS1_MSE; UKF2_MSE; UKS2_MSE; ...
                             EKF_MSE; ERTS_MSE; BS_MSE; ...
                             GHKF_MSE; GHRTS_MSE; CKF_MSE; CRTS_MSE];

  end

  clc
  
  % Calculate means
  means = mean(mse_results,2);
  UKF1_MSE  = means(1);
  UKS1_MSE  = means(2);
  UKF2_MSE  = means(3);
  UKS2_MSE  = means(4);
  EKF_MSE   = means(5);
  ERTS_MSE  = means(6);
  BS_MSE    = means(7);
  GHKF_MSE  = means(8);
  GHRTS_MSE = means(9);
  CKF_MSE   = means(10);
  CRTS_MSE  = means(11);
      
  % Show average results
  fprintf('Average MSE results over %i Monte Carlo runs\n',ns);
  fprintf('UKF1-MSE  = %.4f\n',UKF1_MSE);
  fprintf('UKS1-MSE  = %.4f\n',UKS1_MSE);
  fprintf('UKF2-MSE  = %.4f\n',UKF2_MSE);
  fprintf('UKS2-MSE  = %.4f\n',UKS2_MSE);
  fprintf('EKF-MSE   = %.4f\n',EKF_MSE);
  fprintf('ERTS-MSE  = %.4f\n',ERTS_MSE);
  fprintf('BS-MSE    = %.4f\n',BS_MSE);
  fprintf('GHKF-MSE  = %.4f\n',GHKF_MSE);
  fprintf('GHRTS-MSE = %.4f\n',GHRTS_MSE);
  fprintf('CKF-MSE   = %.4f\n',CKF_MSE);
  fprintf('CRTS-MSE  = %.4f\n',CRTS_MSE);