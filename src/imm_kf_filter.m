
function [x,P,x_i,P_i,mu] = imm_kf_update(x_ip,P_ip,mu_ip,p_ij,ind,dims,A,Q,R,H,Y)
    m = length(x_ip);
    c_j = zeros(1,m);
        
    MM_def = zeros(dims,1);
    PP_def = diag(20*ones(dims,1));

    for j = 1:m
        c_j(j) = 0;
        for i = 1:m
            c_j(j) = c_j(j) + p_ij(i,j).*mu_ip(i);
        end
    end
    
    mu_ij = zeros(m,m);
    for i = 1:m
        for j = 1:m
            mu_ij(i,j) = p_ij(i,j) * mu_ip(i) / c_j(j);
        end
    end
    
    x_0j = cell(1,m);
    for j = 1:m
        x_0j{j} = zeros(dims,1);
        for i = 1:m
            x_0j{j}(ind{i}) = x_0j{j}(ind{i}) + x_ip{i}*mu_ij(i,j);
        end
    end
    
    P_0j = cell(1,m);
    for j = 1:m
        P_0j{j} = zeros(dims,dims);
        %P_0j{j} = PP_def;
        %P_0j{j}(ind{i},ind{i}) = zeros(length(ind{i}),length(ind{i}));
        for i = 1:m
            P_0j{j}(ind{i},ind{i}) = P_0j{j}(ind{i},ind{i}) + mu_ij(i,j)*(P_ip{i} + (x_ip{i}-x_0j{j}(ind{i}))*(x_ip{i}-x_0j{j}(ind{i}))');
        end
    end
    
    x_p = cell(1,m);
    P_p = cell(1,m);
    K_p = cell(1,m);
    x_i = cell(1,m);
    P_i = cell(1,m);
    lambda = zeros(1,m);

    for i = 1:m
        [x_p{i}, P_p{i}] = kf_predict(x_0j{i}(ind{i}),P_0j{i}(ind{i},ind{i}),A{i},Q{i});
        [x_i{i}, P_i{i}, K, IM, IS] = kf_update(x_p{i},P_p{i},Y,H{i},R{i});
        
        lambda(i) = kf_lhood(x_p{i},P_i{i},Y,H{i},R{i});
        %temp = (Y - IM);
        %lambda(i) = exp(-0.5*(temp')*(IS\temp)) / (sqrt(2*pi*det(IS)));        
    end
    
    mu = zeros(1,m); 
    c = sum(lambda.*c_j);
    mu = c_j.*lambda/c;
    
    x = zeros(dims,1);
    P = zeros(dims,dims);
    %P = PP_def;
    for i = 1:m
        x(ind{i}) = x(ind{i}) + mu(i)*x_i{i};
        P(ind{i},ind{i}) = P(ind{i},ind{i}) + mu(i)*(P_i{i} + (x_i{i}-x(ind{i}))*(x_i{i}-x(ind{i}))');
    end
    
    