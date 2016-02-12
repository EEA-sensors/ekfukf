% Cenerate data of discrete-time re-entry mechanics

% Copyright (C) 2005-2006 Simo Särkkä
%               2007      Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

fprintf('[Generating data...]\n');

nsteps = round(200/dt);
x = true_m0 + sqrt(true_P0) * randn(5,1);
X = zeros(size(x,1),nsteps);
Y = zeros(2,nsteps);
T = zeros(1,nsteps);

t = 0;
for k=1:nsteps
    ddt = dt / sim_iter;
    for i=1:sim_iter
        x = reentry_f(x,{ddt,b0,H0,Gm0,R0});
        x = x + L * sqrt(ddt * true_Qc) * randn(3,1);
    end
    
    y = reentry_h(x,{xr,yr}) + diag([sqrt(vr) sqrt(va)]) * randn(2,1);
    t = t + dt;
    X(:,k) = x;
    Y(:,k) = y;
    T(k) = t;
end

aa = 0.02*(-1:0.1:4);
cx = R0 * cos(aa);
cy = R0 * sin(aa);
plot(xr,yr,'o',X(1,:),X(2,:),'--',cx,cy,'-');

