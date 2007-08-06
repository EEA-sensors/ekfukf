% Measurement model function for the random sine signal demo

function Y = ekf_sine_h(x,param)
   f = x(1,:);
   a = x(3,:);

   Y = a.*sin(f);
   if size(x,1) == 7
      Y = Y + x(7,:);
   end
