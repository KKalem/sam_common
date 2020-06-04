function sdot= Fossen6DOF(t, s, m, I_o, r_g, r_b, r_cp, Xuu, Yvv, Zww, Kpp, Mqq, Nrr, W,B, K_T, Q_T, rpm1,rpm2, d_r, d_e)

    x   = s(1);
    y   = s(2);
    z   = s(3);
    theta = s(4);
    phi = s(5);
    psi = s(6);
    u = s(7);
    v = s(8);
    w = s(9);
    p = s(10);
    q = s(11);
    r = s(12);
    
    eta= [x y z theta phi psi]';
    nu=[u v w p q r]'; 

    %cg position
    x_g= r_g(1);
    y_g= r_g(2);
    z_g= r_g(3);
    
    %cb position
    x_b= r_b(1);
    y_b= r_b(2);
    z_b= r_b(3);
    
    %cp position
    x_cp= r_cp(1);
    y_cp= r_cp(2);
    z_cp= r_cp(3);

      
    % Mass and inertia matrix
    M = [m*eye(3,3), -m*Skew(r_g);
        m*Skew(r_g), I_o];

    % Coriolis and centripetal matrix
    nu1 = [u v w]';
    nu2 = [p q r]';
    C_RB = [zeros(3,3), -m*Skew(nu1)-m*Skew(nu2)*Skew(r_g);
            -m*Skew(nu1)+m*Skew(r_g)*Skew(nu2), -Skew(I_o*nu2)];

    % Damping matrix (can be later updated with the lookup-tables)
    D= [Xuu*abs(u) 0 0 0 0 0;
        0 Yvv*abs(v) 0 0 0 0;
        0 0 Zww*abs(w) 0 0 0;
        0 -z_cp*Yvv*v y_cp*Zww*w Kpp*abs(p) 0 0;
        z_cp*Xuu*u 0 -x_cp*Zww*w 0 Mqq*abs(q) 0;
        -y_cp*Xuu*u x_cp*Yvv*v 0 0 0 Nrr*abs(r)];
        
    % currently in Euler angles, convert to Quaternions!
    
    %rotational transform between body and NED
    R_euler = [cos(psi)*cos(theta)   -sin(psi)*cos(phi)+cos(psi)*sin(theta)*sin(phi)   sin(psi)*sin(phi)+cos(psi)*cos(phi)*sin(theta);
           sin(psi)*cos(theta)   cos(psi)*cos(phi)+sin(phi)*sin(theta)*sin(psi)    -cos(psi)*sin(phi)+sin(theta)*sin(psi)*cos(phi);
           -sin(theta)           cos(theta)*sin(phi)                               cos(theta)*cos(phi)];
    T_euler = [1   sin(phi)*tan(theta)   cos(phi)*tan(theta);
           0   cos(phi)              -sin(phi);
           0   sin(phi)/cos(theta)   cos(phi)/cos(theta)];

    J_eta = [R_euler zeros(3,3);
         zeros(3,3) T_euler];
    
    % buoyancy (0 in 2d case on surface)
%     g_eta = [(W-B)*sin(theta);
%              -(W-B)*cos(theta)*sin(phi);
%              -(W-B)*cos(theta)*cos(phi);
%              -(y_g*W-y_b*B)*cos(theta)*cos(phi) + (z_g*W-z_b*B)*cos(theta)*sin(phi);
%              (z_g*W-z_b*B)*sin(theta) + (x_g*W-x_b*B)*cos(theta)*cos(phi);
%              -(x_g*W-x_b*B)*cos(theta)*sin(phi) - (y_g*W-y_b*B)*sin(theta)];
% 
    f_g = [0; 0; W];
    f_b = [0; 0; -B];
    g_eta = [inv(R_euler)*(f_g+f_b);
            Skew(r_g)*inv(R_euler)*f_g + Skew(r_b)*inv(R_euler)*f_b];

    % controls    
    F_T= K_T*[rpm1; rpm2];
    M_T= Q_T*[rpm1; rpm2];
    tau_c = [ F_T* cos(d_e)*cos(d_r);
              -F_T* sin(d_r);
              F_T* sin(d_e)*cos(d_r);
              M_T* cos(d_e)*cos(d_r);
              -M_T* sin(d_r);
              M_T* sin(d_e)*cos(d_r)];

% Kinematics
etadot= J_eta*nu;

% Dynamics
nudot= inv(M)*(tau_c-(C_RB+D)*nu-g_eta);

sdot= [etadot;
        nudot];




