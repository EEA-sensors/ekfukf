  %
  % Compute condition numbers of transition matrices
  % for re-entry simulation with no noise.
  %
  fprintf('[Generating data...]\n');
  
  nsteps = round(200/dt);
  x = true_m0;
  X = zeros(size(x,1),nsteps);
  Y = zeros(2,nsteps);
  T = zeros(1,nsteps);

  C = zeros(1,nsteps);
  t = 0;
  for k=1:nsteps
    ddt = dt / sim_iter;
    for i=1:sim_iter
      A = dreentry_da(x,{ddt,b0,H0,Gm0,R0});
      x = dreentry_a(x,{ddt,b0,H0,Gm0,R0});
    end
    c = cond(A)
    
    t = t + dt;
    C(:,k) = c;
    T(k) = t;
  end
