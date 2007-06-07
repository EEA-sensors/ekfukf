function Y = ekf_demo1_h(x,param)
f = x(1,:);
a = x(3,:);

Y = a.*sin(f);
if size(x,1) == 7
    Y = Y + x(7,:);
end
