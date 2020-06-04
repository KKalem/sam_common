%% Optimal control for hydrobatics: Turbo-turn on the surface

%3DOF model, optimal control on the surface for turbo turn using propulsion
%and thrust vectoring subsystems
% Solution using PMP, MPC 
% Sriharsha Bhat, 05.05.2020

close all
clear all 
clc

addpath([pwd '/mfunctions']);
addpath([pwd '/mfunctions/utils']);
addpath([pwd '/mfunctions/dynamics']);
addpath([pwd '/workfiles']);
addpath([pwd '/workfiles/data']);
addpath([pwd '/workfiles/cad']);


%% Initialize variables for nonlinear model

m= 15.4; % Mass(kg)
Izz = 10; %Moment of inertia 

x_g=0;
y_g=0;
z_g=0;


%Hydrodynamics
Xuu= 1;
Yvv= 100;
Nrr= 150; %15;
x_cp = 0.1;
y_cp = 0;


K_T = 0.1;
rpm = 1000;
d_r = -0.1;


%% Setup and solve optimal control problem


%% Perform time integration to observe the effect of the control

% variable s= [x y psi u v r]
r_g= [x_g, y_g, z_g];

y0= [0 0 0 0 0 0];
tspan = [0 200];

[t,s]= ode23s(@(t,s) Fossen3DOF(t,s,m,Izz,r_g, Xuu, Yvv, Nrr, x_cp, y_cp, K_T,rpm, d_r), tspan, y0);

%% Plot results

x = s(:,1);
y = s(:,2);
psi = s(:,3);
u = s(:,4);
v = s(:,5);
r = s(:,6);

figure(1)
    subplot(2,3,1)
        plot(x,y)
        hold on
        grid on
        xlabel('x(m)')
        ylabel('y(m)')
        title('Positions')

    subplot(2,3,2)
        hold on
        grid on
        plot(t,rad2deg(psi))
        xlabel('t(s)')
        ylabel('psi(°)')
        title('Rotations')
    
    subplot(2,3,3)
        hold on
        grid on
        plot(t,u, t, v)
        xlabel('t(s)')
        ylabel('vel(m/s)')
        title('Linear Vel')
        legend('u','v')
        
    subplot(2,3,4)
        hold on
        grid on
        plot(t,r)
        xlabel('t(s)')
        ylabel('r(rad/s)')
        title('Angular Vel')
        legend('r')

    subplot(2,3,5)
        hold on
        grid on
        plot(t,rad2deg(d_r),'o')
        xlabel('t(s)')
        ylabel('d_r(°)')
        title('rudder')
        
    subplot(2,3,6)
        hold on
        grid on
        plot(t,rpm, 'o')
        xlabel('t(s)')
        ylabel('rpm')
        title('Prop rpm')
    
%% Symbolic model 
    
% Mass and inertia
%syms m Izz

% Cg position
% syms x_g y_g z_g

% States
% syms x y psi u v r
% syms Xuu Yvv Nrr x_cp y_cp

% Controls
% syms K_T rpm d_r  
% 
% % Mass and inertia matrix
% M= [m 0 -m*y_g;
%     0 m m*x_g;
%     -m*y_g m*x_g Izz];
% 
% % Coriolis and centripetal matrix
% C_RB= [0 0 -m*(x_g*r+v);
%        0 0 -m*(y_g*r-u);
%        m*(x_g*r+v) m*(y_g*r-u) 0 ];
% 
% % Damping matrix
% D= [Xuu*abs(u) 0 0;
%     0 Yvv*abs(v) 0;
%     -y_cp*Xuu*u x_cp*Yvv*v Nrr*abs(r)];
% 
% % buoyancy (0 in 2d case on surface)
% g_eta= [0;
%         0;
%         0];
% 
% % controls    
% F_T= K_T*rpm;
% tau_c = [ F_T* cos(d_r);
%           -F_T* sin(d_r);
%           0];
% 
% % state matrices      
% eta = [x y psi]';
% nu = [u v r]';
% 
% %rotational transform between body and NED
% J_eta= [cos(psi) -sin(psi) 0; 
%         sin(psi) cos(psi) 0;
%         0 0 1]; 
% 
% % Kinematics
% etadot= J_eta*nu
% 
% % Dynamics
% nudot= inv(M)*(tau_c-(C_RB+D)*nu-g_eta)