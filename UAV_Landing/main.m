%%  Be name Khoda  %%

%%  clear everything
clc, clear, close all, warning off

%%  load desired theta file
load thetad

%%  sampling time
Ts = 0.1;

%%  platform velocity
pv = 3;

%%  wind speed
Vw = 5;

%%  disturbance type
dist_type = 1;

%%  system linearization
%  name of states
v = sym('v',[6,1]);
dv = sym('dv',[6,1]);
u = sym('u',[6,1]);
dx = dv(1);
dy = dv(2);
dz = dv(3);
dphi = dv(4);
dtheta = dv(5);
dpsi = dv(6);
%  constants
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
%  dynamics
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
%  state eqution
ddv = [f_x+g_x*u(1);f_y+g_y*u(2);f_z+g_z*u(3);...
    f_phi+g_phi*u(4);f_theta+g_theta*u(5);f_psi+g_psi*u(6)];
f = [dv;ddv]';
%   state space matrices
A = eval(subs(jacobian(f,[v;dv]),[v;dv;u],zeros(18,1)));
B = eval(subs(jacobian(f,u),[v;dv;u],zeros(18,1)));
C = [eye(6),zeros(6,6)];
D = zeros(6,6);
%   discrete plant
plant = c2d(ss(A,B,C,D),Ts);

%%  linear MPC
p = 30; % Prediction horizon
m = 3; % Control horizon
mpcobj = mpc(plant,Ts,p,m); clc

%%  nonlinear MPC
nx = 12;
ny = 6;
nu = 6;
nlobj = nlmpc(nx,ny,nu);
nlobj.Ts = Ts;
nlobj.PredictionHorizon = 20;
nlobj.ControlHorizon = 3;
nlobj.Model.StateFcn = 'uavDynamics';
nlobj.Model.OutputFcn = @(x,u,Ts) x(1:6);
nlmdl = 'sim_nmpc_2020a';
nlobj.States(3).Min = 0.1;
nlobj.States(7).Min = -5;
nlobj.States(8).Min = -5;
nlobj.States(9).Min = -5;
nlobj.States(7).Max = 5;
nlobj.States(8).Max = 5;
nlobj.States(9).Max = 5;

%%  RL
load Agent_RL

%%  MPC-RL
load Agent_MPC_RL

%%  run simulations
s4 = sim('sim_nmpc_2020a',30); clc
s3 = sim('sim_mpc_2020a',30); clc
s2 = sim('sim_mpc_rl_2020a',30); clc
s1 = sim('sim_rl_2020a',30); clc

[~,idx] = max([max(s1.t),max(s2.t),max(s3.t),max(s4.t)]);
switch idx
    case 1
        T = s1.t; R = s1.r;
    case 2
        T = s2.t; R = s2.r;
    case 3
        T = s3.t; R = s3.r;
    case 4
        T = s4.t; R = s4.r;
end

%%  results
xmin = min(s1.r(:,1))-5; xmax = max(s1.r(:,1))+5;
ymin = min(s1.r(:,2))-5; ymax = max(s1.r(:,2))+5;
zmin = min(s1.r(:,3))-1; zmax = max(s1.r(:,3))+3;
figure(221)
for k = 1:length(s1.t)
   clf
   hold on
   quad1(:,1:2) =  plotCircle(s1.p(k,1),s1.p(k,2),0.1);
   quad1(:,3) = s1.p(k,3)*ones(size(s1.p(k,1)));
   plat1(:,1:2) =  plotCircle(s1.r(k,1),s1.r(k,2),0.3);
   plat1(:,3) =  1*ones(size(s1.p(k,1)));
   fill3(plat1(:,1),plat1(:,2),plat1(:,3),'r-')
   fill3(quad1(:,1),quad1(:,2),quad1(:,3),'b-')
   xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
   title('MPC'), legend('Platform','Quadcopter')
   xlim([xmin,xmax]), ylim([ymin,ymax]), zlim([zmin,zmax])
   view(-45,50)
   pause(0.001)
end
pause(1)
figure(222)
for k = 1:length(s2.t)
   clf
   hold on
   quad2(:,1:2) =  plotCircle(s2.p(k,1),s2.p(k,2),0.1);
   quad2(:,3) = s2.p(k,3)*ones(size(s2.p(k,1)));
   plat2(:,1:2) =  plotCircle(s2.r(k,1),s2.r(k,2),0.3);
   plat2(:,3) =  1*ones(size(s2.p(k,1)));
   fill3(plat2(:,1),plat2(:,2),plat2(:,3),'r-')
   fill3(quad2(:,1),quad2(:,2),quad2(:,3),'b-')
   xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
   title('NMPC'), legend('Platform','Quadcopter')
   xlim([xmin,xmax]), ylim([ymin,ymax]), zlim([zmin,zmax])
   view(-45,50)
   pause(0.001)
end
pause(1)
figure(223)
for k = 1:length(s3.t)
   clf
   hold on
   quad3(:,1:2) =  plotCircle(s3.p(k,1),s3.p(k,2),0.1);
   quad3(:,3) = s3.p(k,3)*ones(size(s3.p(k,1)));
   plat3(:,1:2) =  plotCircle(s3.r(k,1),s3.r(k,2),0.3);
   plat3(:,3) =  1*ones(size(s3.p(k,1)));
   fill3(plat3(:,1),plat3(:,2),plat3(:,3),'r-')
   fill3(quad3(:,1),quad3(:,2),quad3(:,3),'b-')
   xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
   title('RL'), legend('Platform','Quadcopter')
   xlim([xmin,xmax]), ylim([ymin,ymax]), zlim([zmin,zmax])
   view(-45,50)
   pause(0.001)
end
pause(1)
figure(224)
for k = 1:length(s4.t)
   clf
   hold on
   quad4(:,1:2) =  plotCircle(s4.p(k,1),s4.p(k,2),0.1);
   quad4(:,3) = s4.p(k,3)*ones(size(s4.p(k,1)));
   plat4(:,1:2) =  plotCircle(s4.r(k,1),s4.r(k,2),0.3);
   plat4(:,3) =  1*ones(size(s4.p(k,1)));
   fill3(plat4(:,1),plat4(:,2),plat4(:,3),'r-')
   fill3(quad4(:,1),quad4(:,2),quad4(:,3),'b-')
   xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
   title('MPC-RL'), legend('Platform','Quadcopter')
   xlim([xmin,xmax]), ylim([ymin,ymax]), zlim([zmin,zmax])
   view(-45,50)
   pause(0.001)
end
pause(1)

%%  results
% Create a new figure
figure, hold on

% X position plot
subplot(3,2,1), hold on
plot(T,R(:,1),'linewidth',1.5)
plot(s1.t,s1.p(:,1),'linewidth',1.5)
plot(s2.t,s2.p(:,1),'linewidth',1.5)
plot(s3.t,s3.p(:,1),'linewidth',1.5)
plot(s4.t,s4.p(:,1),'linewidth',1.5)
xlabel('Time (s)')
ylabel('X position (m)')
legend('Ref','MPC','NLMPC','RL','MPC-RL','location','best')

% Y position plot
subplot(3,2,2), hold on
plot(T,R(:,2),'linewidth',1.5)
plot(s1.t,s1.p(:,2),'linewidth',1.5)
plot(s2.t,s2.p(:,2),'linewidth',1.5)
plot(s3.t,s3.p(:,2),'linewidth',1.5)
plot(s4.t,s4.p(:,2),'linewidth',1.5)
xlabel('Time (s)')
ylabel('Y position (m)')
legend('Ref','MPC','NLMPC','RL','MPC-RL','location','best')

% Z position plot
subplot(3,2,3), hold on
plot(T,R(:,3),'linewidth',1.5)
plot(s1.t,s1.p(:,3),'linewidth',1.5)
plot(s2.t,s2.p(:,3),'linewidth',1.5)
plot(s3.t,s3.p(:,3),'linewidth',1.5)
plot(s4.t,s4.p(:,3),'linewidth',1.5)
xlabel('Time (s)')
ylabel('Z position (m)')
legend('Ref','MPC','NLMPC','RL','MPC-RL','location','best')

% X velocity plot
subplot(3,2,4), hold on
plot(s1.t,s1.v(:,1),'linewidth',1.5)
plot(s2.t,s2.v(:,1),'linewidth',1.5)
plot(s3.t,s3.v(:,1),'linewidth',1.5)
plot(s4.t,s4.v(:,1),'linewidth',1.5)
xlabel('Time (s)')
ylabel('X Velocity (m/s)')
legend('MPC','NLMPC','RL','MPC-RL','location','best')

% Y velocity plot
subplot(3,2,5), hold on
plot(s1.t,s1.v(:,2),'linewidth',1.5)
plot(s2.t,s2.v(:,2),'linewidth',1.5)
plot(s3.t,s3.v(:,2),'linewidth',1.5)
plot(s4.t,s4.v(:,2),'linewidth',1.5)
xlabel('Time (s)')
ylabel('Y Velocity (m/s)')
legend('MPC','NLMPC','RL','MPC-RL','location','best')

% Z velocity plot
subplot(3,2,6), hold on
plot(s1.t,s1.v(:,3),'linewidth',1.5)
plot(s2.t,s2.v(:,3),'linewidth',1.5)
plot(s3.t,s3.v(:,3),'linewidth',1.5)
plot(s4.t,s4.v(:,3),'linewidth',1.5)
xlabel('Time (s)')
ylabel('Z Velocity (m/s)')
legend('MPC','NLMPC','RL','MPC-RL','location','best')
