% Jacobian of the measurement model function in the random sine signal demo

function dY = ekf_sine_dh_dx(x,param)
f = x(1,:);
w = x(2,:);
a = x(3,:);
  
dY = [(a.*cos(f))' zeros(size(f,2),1) (sin(f))'];