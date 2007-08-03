function y = dreentry_h2(x,param)
% y = dreentry_h(x,param)

  xr = param{1};
  yr = param{2};
  y = zeros(2,size(x,2));
  y(1,:) = sqrt((x(1,:) - xr).^2 + (x(2,:) - yr).^2);
  y(2,:) = atan2(x(2,:) - yr, x(1,:) - xr);
  
  if size(x,1) == 10
    y(1,:) = y(1,:) + x(9,:);
    y(2,:) = y(2,:) + x(10,:);
  elseif size(x,1) == 7
    y(1,:) = y(1,:) + x(6,:);
    y(2,:) = y(2,:) + x(7,:);
  end  