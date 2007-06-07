function dY = ekf_demo1_dh_dx(x,param)
  f = x(1,:);
  w = x(2,:);
  a = x(3,:);
  
  dY = [(a.*cos(f))' zeros(size(f,2),1) (sin(f))'];