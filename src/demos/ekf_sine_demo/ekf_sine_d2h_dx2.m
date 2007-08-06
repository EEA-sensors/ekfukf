% Hessian of the measurement model function in the random sine signal demo

function df = ekf_sine_d2h_dx2(x,param)
f = x(1);
a = x(3);

df = zeros(1,3,3);
df(1,:,:) = [-a*sin(f) 0 cos(f);
             0         0      0;
             cos(f)    0      0];
