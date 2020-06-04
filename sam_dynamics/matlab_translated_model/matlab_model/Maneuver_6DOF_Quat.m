%% Optimal control for hydrobatics: 6dOF maneuvering

%6DOF fossen model with unit quaternion kinematics  
% Solution using PMP, MPC 
% Sriharsha Bhat, 05.05.2020

close all
clear all 
clc

%% Initialize workspace
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
Nrr = 15;%150;

x_cp = 0.1;
y_cp = 0.00;
z_cp = 0.00;
r_cp = [x_cp,y_cp,z_cp];

% Control actuators
K_T = [0.1 0.1];
Q_T = [0.001 -0.001];

rpm1 = 1000;
rpm2 = 1000;

d_r = 0.1;
d_e = 0.1;


%% Setup and solve optimal control problem


%% Perform time integration to observe the effect of the control

% variable s= [x y z theta phi psi u v w p q r]

ixyz = [0 0 0]; %x, y, z in m
iquat = [1 0 0 0]; % eta0, eps1, eps2, eps3 in unit quaternions
iuvw = [0 0 0]; % u, v, w in m/s
ipqr= [0 0 0]; % p,q,r in rad/s

s0 = [ixyz iquat iuvw ipqr]; % combine initial conditions

tspan = [0 20]; % timespan for ODE integration

[t,s]= ode23s(@(t,s) Fossen6DOF_Quat(t,s,m,I_o,r_g,r_b,r_cp, Xuu,Yvv,Zww,Kpp,Mqq,Nrr, W,B, K_T,Q_T, rpm1,rpm2, d_r, d_e), tspan, s0);

%% Plot results

x   = s(:,1);
y   = s(:,2);
z   = s(:,3);
eta0 = s(:,4);
eps1 = s(:,5);
eps2 = s(:,6);
eps3 = s(:,7);
u = s(:,8);
v = s(:,9);
w = s(:,10);
p = s(:,11);
q = s(:,12);
r = s(:,13);
% theta = s(:,14);
% phi = s(:,15);
% psi = s(:,16);

% Rotations
    quatrot= [eta0 eps1 eps2 eps3];
    eulrot = quat2eul(quatrot, 'XYZ');
    theta= eulrot(:,1);
    phi= eulrot(:,2);
    psi= eulrot(:,3); 

figure(1)
    subplot(2,3,1)
        hold on
        grid on
        plot(t,eta0, t,eps1, t,eps2, t,eps3 )
        xlabel('t(s)')
        ylabel('Rotation quaternions')
        legend('eta0','eps1','eps2','eps3')
        title('Orientations (Quat)')
        
    subplot(2,3,2)
        hold on
        grid on
        plot(t,rad2deg(theta),t,rad2deg(phi),t,rad2deg(psi))
        xlabel('t(s)')
        ylabel('Rotations(°)')
        legend('roll','pitch','yaw')
        title('Orientations (Euler)')
    
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
    
figure(2)
scale = 100;
step =1;
fontsize= 12;
set(gca,'fontsize', fontsize);
%trajectory3(x,-y,-z,phi,-theta,-psi,scale/800,step*4,'Samassembly');
trajectory3(x,y,z,-phi,-theta,-psi,scale/400,step*4,'Samassembly');
