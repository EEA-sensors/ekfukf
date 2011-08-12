
%% Run the bot demo ns times to obtain some reliable results

  % Number of loops
  ns = 100;
  
  % Allocate space for results
  rmse_results = zeros(10,ns);
  
  % Show waitbar
  handle = waitbar(0,'Please wait...');  
  loopstart = tic;
  
  for loop = 1:ns

      figure(1)
      
      % Run the demos (with silen=1, keep_trajectory=1 - not in first)
      ekfs_bot_demo
      ukfs_bot_demo
      ghkfs_bot_demo
      ckfs_bot_demo
      
      rmse_results(:,loop) = [ekf_rmse;   % EKF
                              eks_rmse1;  % ERTS
                              eks_rmse2;  % ETF
                              ukf_rmse;   % UKF
                              uks_rmse1;  % URTS
                              uks_rmse2;  % UTF
                              ghkf_rmse;  % GHKF
                              ghks_rmse1; % GHRTS
                              ckf_rmse;   % CHKF
                              cks_rmse1;];% CRTS
                          
      % Update waitbar and show time remaining
      if rem(k,20)==0 
          secondsleft = min(toc(loopstart),Inf)/loop*(ns-loop);
          waitbar(loop/ns,handle,sprintf('Monte Carlo runs\nTime left: %.0f min %.0f s.', ...
          floor(secondsleft/60),rem(secondsleft,60)))
      end
  end

  % Get rid of the waitbar
  close(handle)
  
  clc
  
  % Calculate means
  means = mean(rmse_results,2);
  
  ekf_rmse   = means(1);   % EKF
  eks_rmse1  = means(2);   % ERTS
  eks_rmse2  = means(3);   % ETF
  ukf_rmse   = means(4);   % UKF
  uks_rmse1  = means(5);   % URTS
  uks_rmse2  = means(6);   % UTF
  ghkf_rmse  = means(7);   % GHKF
  ghrts_rmse = means(8);   % GHRTS
  ckf_rmse   = means(9);   % CKF
  crts_rmse  = means(10);   % CRTS
  

      
  % Show average results
  fprintf('Average RMSE results over %i Monte Carlo runs\n',ns);
  fprintf('  EKF1-RMSE  = %.4f\n',ekf_rmse);
  fprintf('  EKS1-RMSE  = %.4f\n',eks_rmse1);
  fprintf('  ETF1-RMSE  = %.4f\n',eks_rmse2);
  fprintf('  UKF1-RMSE  = %.4f\n',ukf_rmse);
  fprintf('  URTS-RMSE  = %.4f\n',uks_rmse1);
  fprintf('  UTF-RMSE   = %.4f\n',uks_rmse2);
  fprintf('  GHKF-RMSE  = %.4f\n',ghkf_rmse);
  fprintf('  GHRTS-RMSE = %.4f\n',ghrts_rmse);
  fprintf('  CKF-RMSE   = %.4f\n',ckf_rmse);
  fprintf('  CRTS-RMSE  = %.4f\n',crts_rmse);
  
  
  
  ekf_rmse;   % EKF
  eks_rmse1;  % ERTS
  eks_rmse2;  % ETF
  ukf_rmse;   % UKF
  uks_rmse1;  % URTS
  uks_rmse2;  % UTF
  ckf_rmse;   % CKF
  crts_rmse;  % CRTS
  ghkf_rmse;  % GHKF
  ghrts_rmse; % GHRTS
  
  
  
  