%% Optimal control for hydrobatics: 6dOF maneuvering

%6DOF fossen model, 
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

m = 15.4; % Mass(kg)
Ixx = 10;
Iyy = 10;
Izz = 10; %Moment of inertia 
I_o = [Ixx 0 0;
       0 Iyy 0;
       0 0 Izz];

x_g = 0;
y_g = 0;
z_g = 0;
r_g = [x_g, y_g, z_g];

x_b = 0;
y_b = 0;
z_b = 0;
r_b = [x_b, y_b, z_b];

%Weight and buoyancy
W = m*9.81;
B = W;

%Hydrodynamics
Xuu = 1;
Yvv = 100;
Zww = 100;
Kpp = 100;
Mqq = 100;
Nrr = 15;

x_cp = 0.1;
y_cp = 0.0;
z_cp = 0.0;
r_cp = [x_cp,y_cp,z_cp];

% Control actuators
K_T = [0.1 0.1];
Q_T = [0.01 -0.01];

rpm1 = 1000;
rpm2 = 1000;

d_r = 0.1;
d_e = 0.1;


%% Setup and solve optimal control problem


%% Perform time integration to observe the effect of the control

% variable s= [x y z theta phi psi u v w p q r]


y0= [0 0 0 0 0 0 0 0 0 0 0 0];
tspan = [0 20];

[t,s]= ode23s(@(t,s) Fossen6DOF(t,s,m,I_o,r_g,r_b,r_cp, Xuu,Yvv,Zww,Kpp,Mqq,Nrr, W,B, K_T,Q_T, rpm1,rpm2, d_r, d_e), tspan, y0);

%% Plot results

x   = s(:,1);
y   = s(:,2);
z   = s(:,3);
theta = s(:,4);
phi = s(:,5);
psi = s(:,6);
u = s(:,7);
v = s(:,8);
w = s(:,9);
p = s(:,10);
q = s(:,11);
r = s(:,12);

figure(1)
    subplot(2,3,1)
        grid on
        plot3(x,y,z)
        xlabel('x(m)')
        ylabel('y(m)')
        zlabel('z(m)')
        title('Trajectory')

    subplot(2,3,2)
        hold on
        grid on
        plot(t,rad2deg(theta),t,rad2deg(phi),t,rad2deg(psi))
        xlabel('t(s)')
        ylabel('Rotations(°)')
        legend('roll','pitch','yaw')
        title('Orientations')
    
    subplot(2,3,3)
        hold on
        grid on
        plot(t,u,t,v, t,w)
        xlabel('t(s)')
        ylabel('vel(m/s)')
        legend('u','v','w')
        title('linear velocities')

    subplot(2,3,4)
        hold on
        grid on
        plot(t,p,t,q, t,r)
        xlabel('t(s)')
        ylabel('angular vel(m/s)')
        legend('p','q','r')
        title('angular velocities')

    subplot(2,3,5)
        hold on
        grid on
        plot(t,d_e,'o',t,d_r,'o')
        xlabel('t(s)')
        ylabel('thrust vector(°)')
        legend('d_r','d_r')
        title('thrust vector')
        
     subplot(2,3,6)
        hold on
        grid on
        plot(t,rpm1,'o',t,rpm2,'o')
        xlabel('t(s)')
        ylabel('rpm')
        legend('prop1','prop2')
        title('propeller rpm')
    
figure(3)
    scale = 100;
    step =1;
    fontsize= 12;
    set(gca,'fontsize', fontsize);
    %trajectory3(x,-y,-z,phi,-theta,-psi,scale/800,step*4,'Samassembly');
    trajectory3(x,y,z,-phi,-theta,-psi,scale/400,step*4,'Samassembly');
