%EIMM_SMOOTH  EKF based fixed-interval IMM smoother using two IMM-EKF filters.
%
% Syntax:
%   [X_S,P_S,X_IS,P_IS,MU_S] = EIMM_SMOOTH(MM,PP,MM_i,PP_i,MU,p_ij,mu_0j,ind,dims,A,a,a_param,Q,R,H,h,h_param,Y)
%
% In:
%   MM    - Means of forward-time IMM-filter on each time step
%   PP    - Covariances of forward-time IMM-filter on each time step
%   MM_i  - Model-conditional means of forward-time IMM-filter on each time step 
%   PP_i  - Model-conditional covariances of forward-time IMM-filter on each time step
%   MU    - Model probabilities of forward-time IMM-filter on each time step 
%   p_ij  - Model transition probability matrix
%   ind   - Indices of state components for each model as a cell array
%   dims  - Total number of different state components in the combined system
%   A     - Dynamic model matrices for each linear model and Jacobians of each
%           non-linear model's measurement model function as a cell array
%   a     - Cell array containing function handles for dynamic functions
%           for each model having non-linear dynamics
%   a_param - Parameters of a as a cell array.
%   Q     - Process noise matrices for each model as a cell array.
%   R     - Measurement noise matrices for each model as a cell array.
%   H     - Measurement matrices for each linear model and Jacobians of each
%           non-linear model's measurement model function as a cell array
%   h     - Cell array containing function handles for measurement functions
%           for each model having non-linear measurements
%   h_param - Parameters of h as a cell array.
%   Y     - Measurement sequence
%
% Out:
%   X_S  - Smoothed state means for each time step
%   P_S  - Smoothed state covariances for each time step
%   X_IS - Model-conditioned smoothed state means for each time step
%   P_IS - Model-conditioned smoothed state covariances for each time step
%   MU_S - Smoothed model probabilities for each time step
%   
% Description:
%   EKF based two-filter fixed-interval IMM smoother.
%
% See also:
%   EIMM_UPDATE, EIMM_PREDICT

% History:
%   09.01.2008 JH The first official version.
%
% Copyright (C) 2007,2008 Jouni Hartikainen
%
% $Id: imm_update.m 111 2007-11-01 12:09:23Z jmjharti $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.
 
function [x_sk,P_sk,x_sik,P_sik,mu_sk] = eimm_smooth(MM,PP,MM_i,PP_i,MU,p_ij,mu_0j,ind,dims,A,a,a_param,Q,R,H,h,h_param,Y)
    % Default values for mean and covariance
    MM_def = zeros(dims,1);
    PP_def = diag(ones(dims,1));

    % Number of models
    m = length(A);
    
    % Number of measurements
    n = size(Y,2);
    
    % The prior model probabilities for each step
    p_jk = zeros(m,n);
    p_jk(:,1) = mu_0j;
    for i1 = 2:n
        for i2 = 1:m
            p_jk(i2,i1) = sum(p_ij(:,i2).*p_jk(:,i1-1));
        end
    end
     
    % Backward-time transition probabilities
    p_ijb = cell(1,n);
    for k = 1:n
        for i1 = 1:m
            % Normalizing constant
            b_i = sum(p_ij(:,i1).*p_jk(:,k));
            for j = 1:m
                p_ijb{k}(i1,j) = 1/b_i.*p_ij(j,i1).*p_jk(j,k);
            end
        end
    end
    
    % Space for overall smoothed estimates
    x_sk = zeros(dims,n);
    P_sk = zeros(dims,dims,n);
    mu_sk = zeros(m,n);
    
    % Values of smoothed estimates at the last time step.
    x_sk(:,end)   = MM(:,end);
    P_sk(:,:,end) = PP(:,:,end);
    mu_sk(:,end)  = MU(:,end);
    
    % Space for model-conditioned smoothed estimates
    x_sik = cell(m,n);
    P_sik = cell(m,n);
    
    % Values for last time step
    x_sik(:,end) = MM_i(:,end);
    P_sik(:,end) = PP_i(:,end);
    
    % Backward-time estimated model probabilities
    mu_bp = MU(:,end);
    
    % Space for model-conditioned backward-time updated means and covariances
    x_bki = cell(1,m);
    P_bki = cell(1,m);
    
    % Initialize with default values
    for i1 = 1:m
       x_bki{i1} = MM_def;
       x_bki{i1}(ind{i1}) = MM_i{i1,end};
       P_bki{i1} = PP_def;
       P_bki{i1}(ind{i1},ind{i1}) = PP_i{i1,end};
    end
    
    % Space for model-conditioned backward-time predicted means and covariances
    x_kp = cell(1,m);
    P_kp = cell(1,m);
    
    % Initialize with default values
    for i1 = 1:m
       x_kp{i1} = MM_def;
       P_kp{i1} = PP_def;
    end

    for k = n-1:-1:1
        % Space for normalizing constants and conditional model probabilities
        a_j = zeros(1,m);
        mu_bijp = zeros(m,m);
        
        for i2 = 1:m
            % Normalizing constant
            a_j(i2) = sum(p_ijb{k}(:,i2).*mu_bp(:));
             % Conditional model probability
            mu_bijp(:,i2) = 1/a_j(i2).*p_ijb{k}(:,i2).*mu_bp(:); 
            
            % Retrieve the transition matrix or the Jacobian of the dynamic model
            if isnumeric(A{i2})
                A2 = A{i2};
            elseif isstr(A{i2}) | strcmp(class(A{i2}),'function_handle')
                A2 = feval(A{i2},x_bki{i2}(ind{i2}),a_param{i2});
            else
                A2 = A{i2}(x_bki{i2}(ind{i2}),a_param{i2});
            end
            
            % Backward-time EKF prediction step
            [x_kp{i2}(ind{i2}), P_kp{i2}(ind{i2},ind{i2})] = ekf_predict1(x_bki{i2}(ind{i2}),...
                                                              P_bki{i2}(ind{i2},ind{i2}),...
                                                              inv(A2),Q{i2},a{i2},[],a_param{i2});
            
        end 
        
        % Space for mixed predicted mean and covariance
        x_kp0 = cell(1,m);
        P_kp0 = cell(1,m);
        
        % Space for measurement likelihoods
        lhood_j = zeros(1,m);
        
        for i2 = 1:m
            % Initialize with default values
            x_kp0{i2} = MM_def;
            P_kp0{i2} = PP_def;            
            P_kp0{i2}(ind{i2},ind{i2}) = zeros(length(ind{i2}),length(ind{i2}));
            
            % Mix the mean
            for i1 = 1:m
                x_kp0{i2}(ind{i2}) = x_kp0{i2}(ind{i2}) + mu_bijp(i1,i2)*x_kp{i1}(ind{i2});
            end
            
            % Mix the covariance 
            for i1 = 1:m
                P_kp0{i2}(ind{i2},ind{i2}) = P_kp0{i2}(ind{i2},ind{i2}) + mu_bijp(i1,i2)*(P_kp{i1}(ind{i2},ind{i2})+(x_kp{i1}(ind{i2})-x_kp0{i2}(ind{i2}))*(x_kp{i1}(ind{i2})-x_kp0{i2}(ind{i2}))'); 
            end

            % Backward-time EKF update 
            %
            % If the measurement model is linear don't pass h and h_param to ekf_update1
            if isempty(h) | isempty(h{i2})
                [x_bki{i2}(ind{i2}), P_bki{i2}(ind{i2},ind{i2}),K,MUP,S,lhood_j(i2)] = ekf_update1(x_kp0{i2}(ind{i2}),P_kp0{i2}(ind{i2},ind{i2}),Y(:,k),H{i2},R{i2},[],[],[]);
            else
                [x_bki{i2}(ind{i2}), P_bki{i2}(ind{i2},ind{i2}),K,MUP,S,lhood_j(i2)] = ekf_update1(x_kp0{i2}(ind{i2}),P_kp0{i2}(ind{i2},ind{i2}),Y(:,k),H{i2},R{i2},h{i2},[],h_param{i2});
            end
        end
        
        % Normalizing constant
        a_s = sum(lhood_j.*a_j);
        % Updated model probabilities
        mu_bp = 1/a_s.*a_j.*lhood_j;        
        
        % Space for conditional measurement likelihoods
        lhood_ji = zeros(m,m);
        for i1 = 1:m
            for i2 = 1:m
                d_ijk = MM_def;
                D_ijk = PP_def;
                d_ijk = d_ijk + x_kp{i1};
                d_ijk(ind{i2}) = d_ijk(ind{i2}) - MM_i{i2,k};
                PP2 = zeros(dims,dims);
                PP2(ind{i2},ind{i2}) = PP_i{i2,k};
                D_ijk = P_kp{i1} + PP2;

                % Calculate the (approximate) conditional measurement likelihoods
                lhood_ji(i2,i1) = gauss_pdf(d_ijk,0,D_ijk);                
            end
        end
        
        d_j = zeros(m,1);
        for i2 = 1:m
           d_j(i2) = sum(p_ij(i2,:).*lhood_ji(i2,:)); 
        end
        d = sum(d_j.*MU(:,k));
        
        mu_ijsp = zeros(m,m);
        for i1 = 1:m
            for i2 = 1:m
                mu_ijsp(i1,i2) = 1./d_j(i2)*p_ij(i2,i1)*lhood_ji(i2,i1);
            end
        end
                
        mu_sk(:,k) = 1/d.*d_j.*MU(:,k);
        
        % Space for two-step conditional smoothing distributions p(x_k^j|m_{k+1}^i,y_{1:N}),
        % which are a products of two Gaussians
        x_jis = cell(m,m);
        P_jis = cell(m,m);
        for i2 = 1:m
            for i1 = 1:m
                MM1 = MM_def;
                MM1(ind{i2}) = MM_i{i2,k};
                
                PP1 = PP_def;
                PP1(ind{i2},ind{i2}) = PP_i{i2,k};

                iPP1 = inv(PP1);
                iPP2 = inv(P_kp{i1});
                
                % Covariance of the Gaussian product
                P_jis{i2,i1} = inv(iPP1+iPP2);
                % Mean of the Gaussian product
                x_jis{i2,i1} = P_jis{i2,i1}*(iPP1*MM1 + iPP2*x_kp{i1});
            end
        end
        
        % Mix the two-step conditional distributions to yield model-conditioned
        % smoothing distributions.
        for i2 = 1:m
            % Initialize with default values
            x_sik{i2,k} = MM_def;
            P_sik{i2,k} = PP_def;
            P_sik{i2,k}(ind{i2},ind{i2}) = zeros(length(ind{i2}),length(ind{i2}));
            
            % Mixed mean
            for i1 = 1:m
                x_sik{i2,k} = x_sik{i2,k} + mu_ijsp(i1,i2)*x_jis{i2,i1};
            end
            
            % Mixed covariance
            for i1 = 1:m
                P_sik{i2,k} = P_sik{i2,k} + mu_ijsp(i1,i2)*(P_jis{i2,i1} + (x_jis{i2,i1}-x_sik{i2,k})*(x_jis{i2,i1}-x_sik{i2,k})'); 
            end
        end
        
        % Mix the overall smoothed mean
        for i1 = 1:m
            x_sk(:,k) = x_sk(:,k) + mu_sk(i1,k)*x_sik{i1,k};
        end
        
        % Mix the overall smoothed covariance
        for i1 = 1:m
            P_sk(:,:,k) = P_sk(:,:,k) + mu_sk(i1,k)*(P_sik{i1,k} + (x_sik{i1,k}-x_sk(:,k))*(x_sik{i1,k}-x_sk(:,k))');
        end
        
    end
    
