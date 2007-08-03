function dy = dreentry_dh(x,param)
% dy = dreentry_dh(x,param)

  xr = param{1};
  yr = param{2};

  y = zeros(2,size(x,2));
  y(1,:) = sqrt((x(1,:) - xr).^2 + (x(2,:) - yr).^2);
  y(2,:) = atan2(x(2,:) - yr, x(1,:) - xr);
  
  dy = zeros(2,size(x,1));

  dy(1,1) = (x(1) - xr) / y(1);
  dy(1,2) = (x(2) - yr) / y(1);
  dy(2,1) = -(x(2)-yr) / ((x(1) - xr)^2 + (x(2) - yr)^2);
  dy(2,2) = (x(1)-xr) / ((x(1) - xr)^2 + (x(2) - yr)^2);

