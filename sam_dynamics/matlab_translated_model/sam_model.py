#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8


import numpy as np
from scipy.integrate import solve_ivp


def skew(l):
    l = l.flatten()
    l1 = l[1-1]
    l2 = l[2-1]
    l3 = l[3-1]
    return np.array([[0, -l3, l2],
                     [l3, 0, -l1],
                     [-l2, l1, 0]])


def Fossen6DOF_Quat(t, s, m, I_o, r_g, r_b, r_cp, Xuu, Yvv, Zww, Kpp, Mqq, Nrr, W, B, K_T, Q_T, rpm1, rpm2, d_r, d_e):

    x = s[1-1]
    y = s[2-1]
    z = s[3-1]
    eta0 = s[4-1]
    eps1 = s[5-1]
    eps2 = s[6-1]
    eps3 = s[7-1]
    u = s[8-1]
    v = s[9-1]
    w = s[10-1]
    p = s[11-1]
    q = s[12-1]
    r = s[13-1]

    eta= np.array([[x], [y], [z], [eta0], [eps1], [eps2], [eps3]])
    nu= np.array([[u], [v], [w], [p], [q], [r]])

    #cg position
    x_g= r_g[1-1]
    y_g= r_g[2-1]
    z_g= r_g[3-1]

    #cb position
    x_b= r_b[1-1]
    y_b= r_b[2-1]
    z_b= r_b[3-1]

    #cp position
    x_cp= r_cp[1-1]
    y_cp= r_cp[2-1]
    z_cp= r_cp[3-1]

    # Mass and inertia matrix
    M = np.block([[m*np.eye(3,3), -m*skew(r_g)],
                  [m*skew(r_g), I_o]])

    assert M.shape == (6,6), M

    # Coriolis and centripetal matrix
    nu1 = np.array([[u], [v], [w]])
    nu2 = np.array([[p], [q], [r]])
    top_right = -m*skew(nu1)-m*skew(nu2)*skew(r_g)
    bottom_left =-m*skew(nu1) + m*skew(r_g)*skew(nu2)
    Ionu2 = I_o.dot(nu2)
    bottom_right = -skew(Ionu2)
    C_RB = np.block([[np.zeros([3,3]), top_right],
                     [bottom_left, bottom_right]])

    assert C_RB.shape == (6,6), C_RB

    # Damping matrix (can be later updated with the lookup-tables)

    D= np.array([[Xuu*abs(u), 0, 0, 0, 0, 0],
                 [0, Yvv*abs(v), 0, 0, 0, 0],
                 [0, 0, Zww*abs(w), 0, 0, 0],
                 [0, -z_cp*Yvv*abs(v), y_cp*Zww*abs(w), Kpp*abs(p), 0, 0],
                 [z_cp*Xuu*abs(u), 0, -x_cp*Zww*abs(w), 0, Mqq*abs(q), 0],
                 [-y_cp*Xuu*abs(u), x_cp*Yvv*abs(v), 0, 0, 0, Nrr*abs(r)]])

    assert D.shape == (6,6), D

    #rotational transform between body and NED in quaternions
    T_q = 0.5*np.array([[-eps1, -eps2, -eps3],
                        [eta0, -eps3, eps2],
                        [eps3, eta0, -eps1],
                        [-eps2, eps1, eta0]])

    assert T_q.shape == (4,3), T_q

    R_q = np.array([[1-2*(eps2**2+eps3**2),  2*(eps1*eps2-eps3*eta0),  2*(eps1*eps3+eps2*eta0)],
                    [2*(eps1*eps2+eps3*eta0),  1-2*(eps1**2+eps3**2),  2*(eps2*eps3-eps1*eta0)],
                    [2*(eps1*eps3-eps2*eta0),  2*(eps2*eps3+eps1*eta0),  1-2*(eps1**2+eps2**2)]])

    assert R_q.shape == (3,3), R_q

    J_eta = np.block([[R_q, np.zeros([3,3])],
                      [np.zeros([4,3]), T_q]])

    assert J_eta.shape == (7,6), J_eta


    # buoyancy in quaternions
    f_g = np.array([[0], [0], [W]])
    f_b = np.array([[0], [0], [-B]])
    row1 = [np.linalg.inv(R_q).dot((f_g+f_b))]
    row2 = [skew(r_g).dot(np.linalg.inv(R_q)).dot(f_g) + skew(r_b).dot(np.linalg.inv(R_q)).dot(f_b)]
    g_eta = np.block([row1,
                      row2])

    assert g_eta.shape == (6,1), g_eta


    # controls
    F_T= K_T.dot(np.array([[rpm1], [rpm2]]))
    M_T= Q_T.dot(np.array([[rpm1], [rpm2]]))
    tau_c = np.block([ [F_T*np.cos(d_e)*np.cos(d_r)],
                       [-F_T*np.sin(d_r)],
                       [F_T*np.sin(d_e)*np.cos(d_r)],
                       [M_T*np.cos(d_e)*np.cos(d_r)],
                       [-M_T*np.sin(d_r)],
                       [M_T*np.sin(d_e)*np.cos(d_r)] ])

    assert tau_c.shape == (6,1), tau_c

    # Kinematics
    etadot = np.block([J_eta.dot(nu)])

    assert etadot.shape == (7,1), etadot

    # Dynamics
    invM = np.linalg.inv(M)
    nugeta = nu-g_eta
    crbd = C_RB+D
    other = crbd.dot(nugeta)
    other2 = tau_c-other
    nudot = invM.dot(other2)

    assert nudot.shape == (6,1), nudot

    sdot= np.block([ [etadot],
                     [nudot] ])


    return sdot.flatten()




class AUV(object):
    def __init__(self,
                 position,
                 orientation,
                 linear_vel,
                 angular_vel,
                 rpm1=0,
                 rpm2=0,
                 d_r=0,
                 d_e=0
                 ):
        # control variables
        self.rpm1 = rpm1
        self.rpm2 = rpm2
        self.d_r = d_r
        self.d_e = d_e


        assert len(position) == 3, "Position must be a list of [x,y,z]"
        assert len(orientation) == 4, "Orientation must be a list of [w,x,y,z]"
        assert len(linear_vel) == 3, "Linear vel must be a list of [u,v,w]"
        assert len(angular_vel) == 3, "Linear vel must be a list of [p,q,r]"

        self.position = position
        self.orientation = orientation
        self.linear_vel = linear_vel
        self.angular_vel = angular_vel

        # physical desc. of the auv, do not touch unless you know what you're doing
        self._m = 15.4
        self._Ixx = 10
        self._Iyy = 10
        self._Izz = 10
        self._I_o = np.array([[self._Ixx, 0, 0],
                              [0, self._Iyy, 0],
                              [0, 0, self._Izz]])

        self._x_g = 0
        self._y_g = 0
        self._z_g = 0
        self._r_g = np.array([self._x_g, self._y_g, self._z_g])

        self._x_b = 0
        self._y_b = 0
        self._z_b = 0
        self._r_b = np.array([self._x_b, self._y_b, self._z_b])

        #Weight and buoyancy
        self._W = self._m*9.81
        self._B = self._W

        #Hydrodynamics
        self._Xuu = 1
        self._Yvv = 100
        self._Zww = 100
        self._Kpp = 100
        self._Mqq = 100
        self._Nrr = 150

        self._x_cp = 0.1
        self._y_cp = 0.00
        self._z_cp = 0.00
        self._r_cp = np.array([self._x_cp, self._y_cp, self._z_cp])

        # Control actuators
        self._K_T = np.array([0.1, 0.1])
        self._Q_T = np.array([0.001, -0.001])



    def simulate_for(self, seconds, rpm1, rpm2, d_r, d_e):
        assert all( [rpm1 <= 2000, rpm1 >= -2000, rpm2 <= 2000, rpm2 >= -2000] ), "RPMs out of bounds"
        assert all( [d_r <= 0.15, d_r >= -0.15, d_e <= 0.15, d_e >= -0.15] ), "d_r or d_e out of bounds"
        self.rpm1 = rpm1
        self.rpm2 = rpm2
        self.d_r = d_r
        self.d_e = d_e

        tspan = np.array([0, seconds]) # timespan for ODE integration
        ## Perform time integration to observe the effect of the control
        # variable s= [x y z theta phi psi u v w p q r]
        s0 = np.hstack([self.position,
                        self.orientation,
                        self.linear_vel,
                        self.angular_vel]) # combine initial conditions
        print("Solving IVP for s0:{}".format(s0))
        ret = solve_ivp(self.fossen_ts, tspan, s0)
        y = ret['y']

        # just set the _final_ values as out internal state
        self.position = y[0:3]
        self.orientation = y[3:7]
        self.linear_vel = y[7:10]
        self.angular_vel = y[10:13]

        # but return the whole trajectory for interested parties
        return y



    def fossen_ts(self, t, s):
        return Fossen6DOF_Quat(t,
                               s,
                               self._m,
                               self._I_o,
                               self._r_g,
                               self._r_b,
                               self._r_cp,
                               self._Xuu,
                               self._Yvv,
                               self._Zww,
                               self._Kpp,
                               self._Mqq,
                               self._Nrr,
                               self._W,
                               self._B,
                               self._K_T,
                               self._Q_T,
                               self.rpm1,
                               self.rpm2,
                               self.d_r,
                               self.d_e)






pos = [0,0,0]
ori = [1, 0,0,0]
lvel = [0,0,0]
avel = [0,0,0]


import time
t0 = time.time()
auv = AUV(pos, ori, lvel, avel)
y = auv.simulate_for(100, 1000, 1000, 0.1, 0.1)
t1 = time.time()
print('single',t1-t0)



ys =[]

for rpm in np.linspace(100.,2000.,10):
    for d in np.linspace(0.,0.15,3):
        rpm1 = rpm
        rpm2 = rpm
        d_r = d
        d_e = d
        auv = AUV(pos, ori, lvel, avel)
        y = auv.simulate_for(20, rpm1, rpm2, d_r, d_e)
        ys.append(y)
t1 = time.time()
print(t1-t0)

import matplotlib.pyplot as plt
plt.ion()

for y in ys:
    plt.plot(y[0,:], y[1,:])

