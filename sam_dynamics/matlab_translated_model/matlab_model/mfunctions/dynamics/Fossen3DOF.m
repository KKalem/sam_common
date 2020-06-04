function sdot= Fossen3DOF(t,s,m,Izz,r_g, Xuu, Yvv, Nrr, x_cp, y_cp, K_T,rpm, d_r)

    x= s(1);
    y= s(2);
    psi= s(3);
    u= s(4);
    v= s(5);
    r= s(6);
    
    eta= [x y psi]';
    nu=[u v r]'; 

    %cg position
    x_g= r_g(1);
    y_g= r_g(2);
    z_g= r_g(3);
    
    
    % Mass and inertia matrix
    M= [m 0 -m*y_g;
        0 m m*x_g;
        -m*y_g m*x_g Izz];

    % Coriolis and centripetal matrix
    C_RB= [0 0 -m*(x_g*r+v);
           0 0 -m*(y_g*r-u);
           m*(x_g*r+v) m*(y_g*r-u) 0 ];

    % Damping matrix (can be later updated with the lookup-tables)
    D= [Xuu*abs(u) 0 0;
        0 Yvv*abs(v) 0;
        -y_cp*Xuu*abs(u) x_cp*Yvv*abs(v) Nrr*abs(r)];
    
%         D= [Xuu*abs(u) 0 0;
%         0 Yvv*abs(v) 0;
%         -y_cp*Xuu*u x_cp*Yvv*v Nrr*abs(r)];
%     
%     D= [Xuu*abs(u) 0 0;
%         0 Yvv*abs(v) 0;
%         0 0 Nrr*abs(r)];


    % buoyancy (0 in 2d case on surface)
    g_eta= [0;
            0;
            0];

    % controls    
    F_T= K_T*rpm;
    tau_c = [ F_T* cos(d_r);
              -F_T* sin(d_r);
              0];

          
    %rotational transform between body and NED
    J_eta= [cos(psi) -sin(psi) 0; 
            sin(psi) cos(psi) 0;
            0 0 1]; 

% Kinematics
etadot= J_eta*nu;

% Dynamics
nudot= inv(M)*(tau_c-(C_RB+D)*nu-g_eta);

sdot= [etadot;
        nudot];




