  %
  % Compute condition numbers of transition matrices
  % for re-entry simulation with no noise.
  %
  
  % Copyright (C) 2005-2006 Simo Särkkä
  %
  % This software is distributed under the GNU General Public 
  % Licence (version 2 or later); please refer to the file 
  % Licence.txt, included with the software, for details.
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
      A = reentry_df_dx(x,{ddt,b0,H0,Gm0,R0});
      x = reentry_f(x,{ddt,b0,H0,Gm0,R0});
    end
    c = cond(A)
    
    t = t + dt;
    C(:,k) = c;
    T(k) = t;
  end
