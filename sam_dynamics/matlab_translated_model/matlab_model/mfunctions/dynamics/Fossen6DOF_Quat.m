function sdot= Fossen6DOF_Quat(t, s, m, I_o, r_g, r_b, r_cp, Xuu, Yvv, Zww, Kpp, Mqq, Nrr, W,B, K_T, Q_T, rpm1,rpm2, d_r, d_e)

    x   = s(1);
    y   = s(2);
    z   = s(3);
    eta0 = s(4);
    eps1 = s(5);
    eps2 = s(6);
    eps3 = s(7);
    u = s(8);
    v = s(9);
    w = s(10);
    p = s(11);
    q = s(12);
    r = s(13);
    
    eta= [x y z eta0 eps1 eps2 eps3]';
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
        0 -z_cp*Yvv*abs(v) y_cp*Zww*abs(w) Kpp*abs(p) 0 0;
        z_cp*Xuu*abs(u) 0 -x_cp*Zww*abs(w) 0 Mqq*abs(q) 0;
        -y_cp*Xuu*abs(u) x_cp*Yvv*abs(v) 0 0 0 Nrr*abs(r)];
    
%      D= [Xuu*abs(u) 0 0 0 0 0;
%         0 Yvv*abs(v) 0 0 0 0;
%         0 0 Zww*abs(w) 0 0 0;
%         0 -z_cp*Yvv*v y_cp*Zww*w Kpp*abs(p) 0 0;
%         z_cp*Xuu*u 0 -x_cp*Zww*w 0 Mqq*abs(q) 0;
%         -y_cp*Xuu*u x_cp*Yvv*v 0 0 0 Nrr*abs(r)];
    
    %rotational transform between body and NED in quaternions
    T_q = 0.5*[-eps1 -eps2 -eps3;
                eta0 -eps3 eps2;
                eps3 eta0 -eps1;
                -eps2 eps1 eta0];
        
    R_q = [1-2*(eps2^2+eps3^2)   2*(eps1*eps2-eps3*eta0)   2*(eps1*eps3+eps2*eta0);
           2*(eps1*eps2+eps3*eta0)   1-2*(eps1^2+eps3^2)   2*(eps2*eps3-eps1*eta0);
           2*(eps1*eps3-eps2*eta0)   2*(eps2*eps3+eps1*eta0)   1-2*(eps1^2+eps2^2)];

    J_eta = [R_q zeros(3,3);
            zeros(4,3) T_q];
    
    % buoyancy in quaternions
    f_g = [0; 0; W];
    f_b = [0; 0; -B];
    g_eta = [inv(R_q)*(f_g+f_b);
            Skew(r_g)*inv(R_q)*f_g + Skew(r_b)*inv(R_q)*f_b];

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
        nudot;];





