% Fixed-interval IMM smoother using two IMM-filters.
% 
% 
function [x_sk,P_sk,x_sik,P_sik,mu_sk] = imm_smooth(MM,PP,MM_i,PP_i,MU,p_ij,mu_0j,ind,dims,A,Q,R,H,Y)

    MM_def = zeros(dims,1);
    PP_def = diag(100*ones(dims,1));
    %PP_def = zeros(dims,dims);
    
    m = length(A);
    n = size(Y,2);
    
    p_jk = zeros(m,n);
    p_jk(:,1) = mu_0j;
    for i1 = 2:n
        for i2 = 1:m
            p_jk(i2,i1) = sum(p_ij(:,i2).*p_jk(:,i1-1));
        end
    end
    
    p_ijb = cell(1,n);
    for k = 1:n
        for i1 = 1:m
            b_i = sum(p_ij(:,i1).*p_jk(:,k));
            for j = 1:m
                p_ijb{k}(i1,j) = 1/b_i.*p_ij(j,i1).*p_jk(j,k);
            end
        end
    end
    
    x_sk = zeros(dims,n);
    P_sk = zeros(dims,dims,n);
    mu_sk = zeros(m,n);
    
    x_sk(:,end)   = MM(:,end);
    P_sk(:,:,end) = PP(:,:,end);
    mu_sk(:,end)  = MU(:,end);
    
    x_sik = cell(m,n);
    P_sik = cell(m,n);
    
    x_sik(:,end) = MM_i(:,end);
    P_sik(:,end) = PP_i(:,end);
    
    mu_bp = MU(:,end);
    
    x_bki = MM_i(:,end);
    P_bki = PP_i(:,end);
    
    x_kp = cell(1,m);
    P_kp = cell(1,m);
    for i1 = 1:m
       x_kp{i1} = MM_def;
       P_kp{i1} = PP_def;
    end
    for k = n-1:-1:1
        a_j = zeros(1,m);
        mu_bijp = zeros(m,m);
        
        for i2 = 1:m
            a_j(i2) = sum(p_ijb{k}(:,i2).*mu_bp(:));
            mu_bijp(:,i2) = 1/a_j(i2).*p_ijb{k}(:,i2).*mu_bp(:); 

            % Backward prediction
            x_kp{i2}(ind{i2}) = A{i2}\x_bki{i2};
            P_kp{i2}(ind{i2},ind{i2}) = (A{i2}\(P_bki{i2}+Q{i2}))/(A{i2}');
        end

        x_kp0 = cell(1,m);
        P_kp0 = cell(1,m);
        lhood_j = zeros(1,m);
        for i2 = 1:m
            x_kp0{i2} = MM_def;
            P_kp0{i2} = PP_def;            
            P_kp0{i2}(ind{i2},ind{i2}) = zeros(length(ind{i2}),length(ind{i2}));
            for i1 = 1:m
                x_kp0{i2}(ind{i2}) = x_kp0{i2}(ind{i2}) + mu_bijp(i1,i2)*x_kp{i1}(ind{i2});
            end
            for i1 = 1:m
                P_kp0{i2}(ind{i2},ind{i2}) = P_kp0{i2}(ind{i2},ind{i2}) + mu_bijp(i1,i2)*(P_kp{i1}(ind{i2},ind{i2})+(x_kp{i1}(ind{i2})-x_kp0{i2}(ind{i2}))*(x_kp{i1}(ind{i2})-x_kp0{i2}(ind{i2}))'); 
            end
            z = Y(:,k) - H{i2} * x_kp0{i2}(ind{i2});
            S = H{i2} * P_kp0{i2}(ind{i2},ind{i2}) * (H{i2}') + R{i2};
            K = P_kp0{i2}(ind{i2},ind{i2}) * ((H{i2}') / S);

            x_bki{i2}(ind{i2}) = x_kp0{i2}(ind{i2}) + K * z;
            P_bki{i2}(ind{i2},ind{i2}) = (eye(size(H{i2},2),size(H{i2},2))-K*H{i2})*P_kp0{i2}(ind{i2},ind{i2});
            lhood_j(i2) = gauss_pdf(z,0,S);
        end
        a = sum(lhood_j.*a_j);
        mu_bp = 1/a.*a_j.*lhood_j;        
        
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
                %D_ijk =  D_ijk + P_kp{i1};
                %D_ijk(ind{i2},ind{i2}) =  D_ijk(ind{i2},ind{i2}) - PP_i{i2,k};
                lhood_ji(i2,i1) = gauss_pdf(d_ijk,0,D_ijk);                
            end
        end
        
        d_j = zeros(m,1);
        for i2 = 1:m
           d_j(i2) = sum(p_ij(i2,:).*lhood_ji(i2,:)); 
        end
        d = sum(d_j.*MU(:,k));
        
        mu_ijsp = zeros(m,m);
% $$$         for i1 = 1:m
% $$$             mu_ijsp(i1,:) = 1./d_j'.*p_ij(i1,:).*lhood_ji(:,i1)';
% $$$         end

        for i1 = 1:m
            for i2 = 1:m
                mu_ijsp(i1,i2) = 1./d_j(i2)*p_ij(i2,i1)*lhood_ji(i2,i1);
            end
        end
                
        mu_sk(:,k) = 1/d.*d_j.*MU(:,k);

        x_jis = cell(m,m);
        P_jis = cell(m,m);
        for i2 = 1:m
            for i1 = 1:m
                MM1 = MM_def;
                MM1(ind{i2}) = MM_i{i2,k};
                
                PP1 = PP_def;
                PP1(ind{i2},ind{i2}) = PP_i{i2,k};

                %iPP1 = PP_def;
                %iPP2 = PP_def;
                %iPP1(ind{i2},ind{i2}) = inv(PP_i{i2,k});
                %iPP2(ind{i1},ind{i1}) = inv(P_kp{i1}(ind{i1},ind{i1}));

                iPP1 = inv(PP1);
                iPP2 = inv(P_kp{i1});
                
                P_jis{i2,i1} = inv(iPP1+iPP2);
                x_jis{i2,i1} = P_jis{i2,i1}*(iPP1*MM1 + iPP2*x_kp{i1});
            end
        end
        
        for i2 = 1:m
            x_sik{i2,k} = MM_def;
            P_sik{i2,k} = PP_def;
            P_sik{i2,k}(ind{i2},ind{i2}) = zeros(length(ind{i2}),length(ind{i2}));
            for i1 = 1:m
                x_sik{i2,k} = x_sik{i2,k} + mu_ijsp(i1,i2)*x_jis{i2,i1};
            end
            for i1 = 1:m
                P_sik{i2,k} = P_sik{i2,k} + mu_ijsp(i1,i2)*(P_jis{i2,i1} + (x_jis{i2,i1}-x_sik{i2,k})*(x_jis{i2,i1}-x_sik{i2,k})'); 
            end
        end

        for i1 = 1:m
            x_sk(:,k) = x_sk(:,k) + mu_sk(i1,k)*x_sik{i1,k};
        end
        for i1 = 1:m
            P_sk(:,:,k) = P_sk(:,:,k) + mu_sk(i1,k)*(P_sik{i1,k} + (x_sik{i1,k}-x_sk(:,k))*(x_sik{i1,k}-x_sk(:,k))');
        end
        
    end
    
