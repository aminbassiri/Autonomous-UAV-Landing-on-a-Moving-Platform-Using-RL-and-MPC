%% clear everything
clc, clear, close all, warning off

%% load desired theta file
load thetad

%% sampling time
Ts = 0.1;

%% platform velocity
pv = 10;

%% wind speed
Vw = 3;

%% disturbance type
dist_type = 1;

%% system linearization
% name of states
v = sym('v',[6,1]);
dv = sym('dv',[6,1]);
u = sym('u',[6,1]);
dx = dv(1);
dy = dv(2);
dz = dv(3);
dphi = dv(4);
dtheta = dv(5);
dpsi = dv(6);

% constants
K1 = 1e-3; K2 = 1e-3; K3 = 1e-3; K4 = 1e-3; K5 = 1e-3; K6 = 1e-3;
m = 2.85;
Ix = 0.0552; Iy = 0.0552; Iz = 0.1104;

% dynamics: f & g
f_x = -K1/m*dx; g_x = 1;
f_y = -K2/m*dy; g_y = 1;
f_z = -K3/m*dz; g_z = 1;
f_phi = -K4/Ix*dphi + (Iy - Iz)/Ix*dtheta*dpsi; g_phi = 1/Ix;
f_theta = -K5/Iy*dtheta + (Iz - Ix)/Iy*dphi*dpsi; g_theta = 1/Iy;
f_psi = -K6/Iz*dpsi + (Ix - Iy)/Iz*dphi*dtheta; g_psi = 1/Iz;

% state equation
ddv = [f_x+g_x*u(1); f_y+g_y*u(2); f_z+g_z*u(3); ...
       f_phi+g_phi*u(4); f_theta+g_theta*u(5); f_psi+g_psi*u(6)];
f = [dv;ddv]';

% state space matrices
A = eval(subs(jacobian(f,[v;dv]),[v;dv;u],zeros(18,1)));
B = eval(subs(jacobian(f,u),[v;dv;u],zeros(18,1)));
C = [eye(6), zeros(6,6)];
D = zeros(6,6);

% discrete plant
plant = c2d(ss(A,B,C,D),Ts);

%% linear MPC
p = 30; % Prediction horizon
m = 3; % Control horizon
mpcobj = mpc(plant,Ts,p,m); clc

%% nonlinear MPC
nx = 12; ny = 6; nu = 6;
nlobj = nlmpc(nx,ny,nu);
nlobj.Ts = Ts;
nlobj.PredictionHorizon = 20;
nlobj.ControlHorizon = 3;
nlobj.Model.StateFcn = 'uavDynamics';
nlobj.Model.OutputFcn = @(x,u,Ts) x(1:6);
nlmdl = 'sim_nmpc_2020a';
nlobj.States(3).Min = 0.1;
nlobj.States(7).Min = -5; nlobj.States(8).Min = -5; nlobj.States(9).Min = -5;
nlobj.States(7).Max = 5; nlobj.States(8).Max = 5; nlobj.States(9).Max = 5;

%% RL
load Agent_RL

%% MPC-RL
load Agent_MPC_RL

%% run simulations
s4 = sim('sim_nmpc_2020a',30); clc
s3 = sim('sim_mpc_2020a',30); clc
s2 = sim('sim_mpc_rl_2020a',30); clc
s1 = sim('sim_rl_2020a',30); clc

[~,idx] = max([max(s1.t), max(s2.t), max(s3.t), max(s4.t)]);
switch idx
    case 1, T = s1.t; R = s1.r;
    case 2, T = s2.t; R = s2.r;
    case 3, T = s3.t; R = s3.r;
    case 4, T = s4.t; R = s4.r;
end

%% Plot all simulations in one figure with legends
figure('Position', [100, 100, 1200, 800])
colors = {'b', 'b', 'b', 'b'}; % Define trajectory colors for quadcopter

xmin = min([min(s1.r(:,1)), min(s2.r(:,1)), min(s3.r(:,1)), min(s4.r(:,1))])-5;
xmax = max([max(s1.r(:,1)), max(s2.r(:,1)), max(s3.r(:,1)), max(s4.r(:,1))])+5;
ymin = min([min(s1.r(:,2)), min(s2.r(:,2)), min(s3.r(:,2)), min(s4.r(:,2))])-5;
ymax = max([max(s1.r(:,2)), max(s2.r(:,2)), max(s3.r(:,2)), max(s4.r(:,2))])+5;
zmin = min([min(s1.r(:,3)), min(s2.r(:,3)), min(s3.r(:,3)), min(s4.r(:,3))])-1;
zmax = max([max(s1.r(:,3)), max(s2.r(:,3)), max(s3.r(:,3)), max(s4.r(:,3))])+3;

% MPC Simulation
subplot(2,2,1), hold on
h_platform = drawPlatform(s1.r(end,1), s1.r(end,2), 1);
h_quad = drawQuadcopter(s1.p(end,1), s1.p(end,2), s1.p(end,3));
h_traj_quad = plot3(s1.p(:,1), s1.p(:,2), s1.p(:,3), [colors{1} '--'], 'LineWidth', 1.5);
h_traj_platform = plot3(s1.r(:,1), s1.r(:,2), s1.r(:,3), 'r:', 'LineWidth', 2); % Platform trajectory
xlabel('x (m)', 'FontSize', 25), ylabel('y (m)', 'FontSize', 25), zlabel('z (m)', 'FontSize', 25)
set(gca, 'FontSize', 25)  % Set font size for axes
title('MPC', 'FontSize', 25)
grid on, view(-45,50)
xlim([xmin, xmax]), ylim([ymin, ymax]), zlim([zmin, zmax])
customLegend(h_quad, h_platform, h_traj_quad, h_traj_platform);

% NMPC Simulation
subplot(2,2,2), hold on
h_platform = drawPlatform(s2.r(end,1), s2.r(end,2), 1);
h_quad = drawQuadcopter(s2.p(end,1), s2.p(end,2), s2.p(end,3));
h_traj_quad = plot3(s2.p(:,1), s2.p(:,2), s2.p(:,3), [colors{2} '--'], 'LineWidth', 1.5);
h_traj_platform = plot3(s2.r(:,1), s2.r(:,2), s2.r(:,3), 'r:', 'LineWidth', 2); % Platform trajectory
xlabel('x (m)', 'FontSize', 25), ylabel('y (m)', 'FontSize', 25), zlabel('z (m)', 'FontSize', 25)
set(gca, 'FontSize', 25)  % Set font size for axes
title('NMPC', 'FontSize', 25)
grid on, view(-45,50)
xlim([xmin, xmax]), ylim([ymin, ymax]), zlim([zmin, zmax])
customLegend(h_quad, h_platform, h_traj_quad, h_traj_platform);

% RL Simulation
subplot(2,2,3), hold on
h_platform = drawPlatform(s3.r(end,1), s3.r(end,2), 1);
h_quad = drawQuadcopter(s3.p(end,1), s3.p(end,2), s3.p(end,3));
h_traj_quad = plot3(s3.p(:,1), s3.p(:,2), s3.p(:,3), [colors{3} '--'], 'LineWidth', 1.5);
h_traj_platform = plot3(s3.r(:,1), s3.r(:,2), s3.r(:,3), 'r:', 'LineWidth', 2); % Platform trajectory
xlabel('x (m)', 'FontSize', 25), ylabel('y (m)', 'FontSize', 25), zlabel('z (m)', 'FontSize', 25)
set(gca, 'FontSize', 25)  % Set font size for axes
title('RL', 'FontSize', 25)
grid on, view(-45,50)
xlim([xmin, xmax]), ylim([ymin, ymax]), zlim([zmin, zmax])
customLegend(h_quad, h_platform, h_traj_quad, h_traj_platform);

% MPC-RL Simulation
subplot(2,2,4), hold on
h_platform = drawPlatform(s4.r(end,1), s4.r(end,2), 1);
h_quad = drawQuadcopter(s4.p(end,1), s4.p(end,2), s4.p(end,3));
h_traj_quad = plot3(s4.p(:,1), s4.p(:,2), s4.p(:,3), [colors{4} '--'], 'LineWidth', 1.5);
h_traj_platform = plot3(s4.r(:,1), s4.r(:,2), s4.r(:,3), 'r:', 'LineWidth', 2); % Platform trajectory
xlabel('x (m)', 'FontSize', 25), ylabel('y (m)', 'FontSize', 25), zlabel('z (m)', 'FontSize', 25)
set(gca, 'FontSize', 25)  % Set font size for axes
title('MPC-RL', 'FontSize', 25)
grid on, view(-45,50)
xlim([xmin, xmax]), ylim([ymin, ymax]), zlim([zmin, zmax])
customLegend(h_quad, h_platform, h_traj_quad, h_traj_platform);

%% Save High-Quality Figures for LaTeX
print('quad_trajectory_simulation', '-dpdf', '-r300');  % High-resolution PDF for LaTeX
print('quad_trajectory_simulation', '-depsc', '-r300');  % High-resolution EPS for LaTeX

%% Results Comparison Plot (Position and Velocity)
figure('Position', [100, 100, 1200, 700]) % Larger figure window for clarity

% X position plot
subplot(3,2,1), hold on
plot(T,R(:,1),'k--','LineWidth',2) % Reference line
plot(s1.t,s1.p(:,1),'b','LineWidth',2) % MPC - Blue
plot(s2.t,s2.p(:,1),'g','LineWidth',2) % NMPC - Green
plot(s3.t,s3.p(:,1),'r','LineWidth',2) % RL - Red
plot(s4.t,s4.p(:,1),'m','LineWidth',2) % MPC-RL - Magenta
xlabel('Time (s)', 'FontSize', 30), ylabel('X position (m)', 'FontSize', 25)
legend('Ref','MPC','NMPC','RL','MPC-RL', 'FontSize', 16)
grid on

% Y position plot
subplot(3,2,2), hold on
plot(T,R(:,2),'k--','LineWidth',2)
plot(s1.t,s1.p(:,2),'b','LineWidth',2) % MPC - Blue
plot(s2.t,s2.p(:,2),'g','LineWidth',2) % NMPC - Green
plot(s3.t,s3.p(:,2),'r','LineWidth',2) % RL - Red
plot(s4.t,s4.p(:,2),'m','LineWidth',2) % MPC-RL - Magenta
xlabel('Time (s)', 'FontSize', 30), ylabel('Y position (m)', 'FontSize', 25)
legend('Ref','MPC','NMPC','RL','MPC-RL', 'FontSize', 16)
grid on

% Z position plot
subplot(3,2,3), hold on
plot(T,R(:,3),'k--','LineWidth',2)
plot(s1.t,s1.p(:,3),'b','LineWidth',2) % MPC - Blue
plot(s2.t,s2.p(:,3),'g','LineWidth',2) % NMPC - Green
plot(s3.t,s3.p(:,3),'r','LineWidth',2) % RL - Red
plot(s4.t,s4.p(:,3),'m','LineWidth',2) % MPC-RL - Magenta
xlabel('Time (s)', 'FontSize', 30), ylabel('Z position (m)', 'FontSize', 25)
legend('Ref','MPC','NMPC','RL','MPC-RL', 'FontSize', 16)
grid on

%% Save Position and Velocity Plots for LaTeX
% High-quality PDF and EPS export
print('position_velocity_comparison', '-dpdf', '-r300');
print('position_velocity_comparison', '-depsc', '-r300');

% X velocity plot
subplot(3,2,4), hold on
plot(s1.t,s1.v(:,1),'b','LineWidth',2) % MPC - Blue
plot(s2.t,s2.v(:,1),'g','LineWidth',2) % NMPC - Green
plot(s3.t,s3.v(:,1),'r','LineWidth',2) % RL - Red
plot(s4.t,s4.v(:,1),'m','LineWidth',2) % MPC-RL - Magenta
xlabel('Time (s)', 'FontSize', 30), ylabel('X Velocity (m/s)', 'FontSize', 25)
legend('MPC','NMPC','RL','MPC-RL', 'FontSize', 16)
grid on

% Y velocity plot
subplot(3,2,5), hold on
plot(s1.t,s1.v(:,2),'b','LineWidth',2) % MPC - Blue
plot(s2.t,s2.v(:,2),'g','LineWidth',2) % NMPC - Green
plot(s3.t,s3.v(:,2),'r','LineWidth',2) % RL - Red
plot(s4.t,s4.v(:,2),'m','LineWidth',2) % MPC-RL - Magenta
xlabel('Time (s)', 'FontSize', 30), ylabel('Y Velocity (m/s)', 'FontSize', 25)
legend('MPC','NMPC','RL','MPC-RL', 'FontSize', 16)
grid on

% Z velocity plot
subplot(3,2,6), hold on
plot(s1.t,s1.v(:,3),'b','LineWidth',2) % MPC - Blue
plot(s2.t,s2.v(:,3),'g','LineWidth',2) % NMPC - Green
plot(s3.t,s3.v(:,3),'r','LineWidth',2) % RL - Red
plot(s4.t,s4.v(:,3),'m','LineWidth',2) % MPC-RL - Magenta
xlabel('Time (s)', 'FontSize', 30), ylabel('Z Velocity (m/s)', 'FontSize', 25)
legend('MPC','NMPC','RL','MPC-RL', 'FontSize', 16)
grid on

% Save the figure
print('velocity_comparison', '-dpdf', '-r300');  % Save velocity comparison

%% Custom Legend Function
function customLegend(h_quad, h_platform, h_traj_quad, h_traj_platform)
    % Create custom patches for the platform and quadcopter
    patch_platform = plot(NaN, NaN, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Platform (red circle)
    line_quadcopter_marker = plot(NaN, NaN, 'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b');  % Small quadcopter
    line_traj_quad = plot(NaN, NaN, 'b--', 'LineWidth', 1.5);  % Quadcopter trajectory
    line_traj_platform = plot(NaN, NaN, 'r:', 'LineWidth', 2); % Platform trajectory
    
    legend([line_quadcopter_marker, patch_platform, line_traj_quad, line_traj_platform], ...
           'Quadcopter', 'Platform', 'Quadcopter Trajectory', 'Platform Trajectory', 'Location', 'best', 'FontSize', 20);
end

%% Function to draw a simple 3D quadcopter
function h_quad = drawQuadcopter(x, y, z)
    armLength = 0.5; % Length of quadcopter arms
    propellerRadius = 0.2; % Radius of propellers
    h1 = plot3([x - armLength, x + armLength], [y, y], [z, z], 'k-', 'LineWidth', 2); % X-axis arm
    plot3([x, x], [y - armLength, y + armLength], [z, z], 'k-', 'LineWidth', 2); % Y-axis arm
    theta = linspace(0, 2*pi, 100);
    for i = [-1, 1]
        for j = [-1, 1]
            xp = propellerRadius * cos(theta) + x + i*armLength/2;
            yp = propellerRadius * sin(theta) + y + j*armLength/2;
            zp = z * ones(size(xp));
            fill3(xp, yp, zp, 'b', 'FaceAlpha', 1.0); % Propellers
        end
    end
    h_quad = h1; % Return a handle for the quadcopter arm to use in the legend
end

%% Function to draw a circular platform
function h_platform = drawPlatform(x, y, z)
    platformRadius = 0.9; % Platform radius
    theta = linspace(0, 2*pi, 100);
    xp = platformRadius * cos(theta) + x;
    yp = platformRadius * sin(theta) + y;
    zp = z * ones(size(xp));
    h_platform = fill3(xp, yp, zp, 'r', 'FaceAlpha', 0.4); % Platform as a filled circle
end
