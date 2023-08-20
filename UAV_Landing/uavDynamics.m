function dv = uavDynamics(v,u)
%%  name of states
dx = v(7);
dy = v(8);
dz = v(9);
dphi = v(10);
dtheta = v(11);
dpsi = v(12);

%%  constants
K1 = 1e-3;
K2 = 1e-3;
K3 = 1e-3;
K4 = 1e-3;
K5 = 1e-3;
K6 = 1e-3;
m = 2.85;
Ix = 0.0552;
Iy = 0.0552;
Iz = 0.1104;

%%  dynamics
%   f & g
f_x = -K1/m*dx;
g_x = 1;
%   f & g
f_y = -K2/m*dy;
g_y = 1;
%   f & g
f_z = -K3/m*dz;
g_z = 1;
%   f & g
f_phi = -K4/Ix*dphi + (Iy - Iz)/Ix*dtheta*dpsi;
g_phi = 1/Ix;
%   f & g
f_theta = -K5/Iy*dtheta + (Iz - Ix)/Iy*dphi*dpsi;
g_theta = 1/Iy;
%   f & g
f_psi = -K6/Iz*dpsi + (Ix - Iy)/Iz*dphi*dtheta;
g_psi = 1/Iz;

%%  state eqution
dv = [v(7:12);f_x+g_x*u(1);f_y+g_y*u(2);f_z+g_z*u(3);...
    f_phi+g_phi*u(4);f_theta+g_theta*u(5);f_psi+g_psi*u(6)];
end