
%% Run the reentry demo ns times to obtain some reliable results

  % Number of loops
  ns = 100;
  
  % Allocate space for results
  rmse_results = zeros(11,ns);
  
  % Show waitbar
  handle = waitbar(0,'Please wait...');  
  loopstart = tic;
  
  for loop = 1:ns

      figure(1)
      
      % Run the demos (with silen=1, keep_trajectory=1 - not in first)
      reentry_demo;
      
      rmse_results(:,loop) = [ekf_rmse;      % EKF
                              eks_rmse1;     % ERTS
                              eks_rmse2;     % ETF
                              ukf_rmse;      % UKF
                              uks_rmse1;     % URTS 1
                              uks_rmse1b;    % URTS 2
                              uks_rmse2;     % UTF
                              ghkf_rmse;     % GHKF
                              ghrts_rmse1;   % GHRTS
                              ckf_rmse;      % CHKF
                              crts_rmse1;]; % CRTS
                          
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
  uks_rmse1  = means(5);   % URTS 1
  uks_rmse1b = means(6);   % URTS 2
  uks_rmse2  = means(7);   % UTF
  ghkf_rmse  = means(8);   % GHKF
  ghrts_rmse = means(9);   % GHRTS
  ckf_rmse   = means(10);   % CKF
  crts_rmse  = means(11);   % CRTS
  

      
  % Show average results
  fprintf('Average RMSE results over %i Monte Carlo runs\n',ns);
  fprintf('EKF-RMSE   = %.6f\n',ekf_rmse);
  fprintf('ERTS-RMSE  = %.6f\n',eks_rmse1);
  fprintf('ETF-RMSE   = %.6f\n',eks_rmse2);
  fprintf('UKF-RMSE   = %.6f\n',ukf_rmse);
  fprintf('URTS1-RMSE = %.6f\n',uks_rmse1);
  fprintf('URTS2-RMSE = %.6f\n',uks_rmse1b);
  fprintf('UTF-RMSE   = %.6f\n',uks_rmse2);
  fprintf('GHKF-RMSE  = %.6f\n',ghkf_rmse);
  fprintf('GHRTS-RMSE = %.6f\n',ghrts_rmse1);
  fprintf('CKF-RMSE   = %.6f\n',ckf_rmse);
  fprintf('CRTS-RMSE  = %.6f\n',crts_rmse1);
  
  