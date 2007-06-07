function x = dreentry_ia(x,param)
% x = dreentry_ia(x,param)

  y = dreentry_a(x,param);
  x = 2 * x(1:size(y,1)) - y;

