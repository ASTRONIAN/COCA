#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  19 20:45:08 2023

@author: Subramanian
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
#plt.style.use('dark_background')

import planetary_data as pd
import tools as t


class orbpro:
    def __init__(self,state0,tspan,dt,coes=False,deg=True,cb=pd.earth):

        if coes:
            self.r0,self.v0,_ = t.coes2rv(state0,deg=deg,mu=cb['mu'])
        else:
            self.r0=state0[:3]
            self.v0=state0[3:]

        self.y0=self.r0.tolist()+self.v0.tolist()
#        self.y0=np.array(self.r0.tolist()+self.v0.tolist())
        self.tspan=tspan
        self.dt=dt
        self.cb=cb

        self.n_steps = int(np.ceil(self.tspan/self.dt))

        self.ts = np.zeros((self.n_steps,1))
        self.y  = np.zeros((self.n_steps,6))
#        self.y0 = self.r0 + self.v0
        self.ts[0] = 0
        self.y[0,:] = np.array(self.y0)
        self.step = 1

        self.solver = ode(self.diffy_q)
        self.solver.set_integrator('lsoda')
        self.solver.set_initial_value(self.y0,0)

        self.propagate_orbit()

    def propagate_orbit(self):

        while self.solver.successful() and self.step<self.n_steps:
            self.solver.integrate(self.solver.t+self.dt)
            self.ts[self.step]=self.solver.t
            self.y[self.step]=self.solver.y
            self.step+=1

        self.rs=self.y[:,:3]
        self.vs=self.y[:,3:]
#        plot(rs)

    def diffy_q(self,t,y):
        rx,ry,rz,vx,vy,vz=y
        r=np.array([rx,ry,rz])
        v=np.array([vx,vy,vz])

        norm_r=np.linalg.norm(r)

        ax,ay,az=-r*self.cb["mu"]/norm_r**3

        return [vx,vy,vz,ax,ay,az]

    def plot_3d(self,show_plot=False,save_plot=False,title='Traj'):
        fig = plt.figure(figsize=(16,8))
        ax = fig.add_subplot(111,projection='3d')

        ax.plot(self.rs[:,0],self.rs[:,1],self.rs[:,2],'w',label= 'Trajectory')
        ax.plot([self.rs[0,0]],[self.rs[0,1]],[self.rs[0,2]],'wo',label = 'Initial Position')

        _u,_v = np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
        _x = self.cb['radius']*np.cos(_u)*np.sin(_v)
        _y = self.cb['radius']*np.sin(_u)*np.sin(_v)
        _z = self.cb['radius']*np.cos(_v)
        ax.plot_surface(_x,_y,_z, cmap='Blues')

        l = self.cb['radius']*2
        x,y,z = [[0,0,0],[0,0,0],[0,0,0]]
        u,v,w = [[1,0,0],[0,1,0],[0,0,1]]
        ax.quiver(x,y,z,u,v,w,color='k')

        max_val = np.max(np.abs(self.rs))

        ax.set_xlim([-max_val, max_val])
        ax.set_ylim([-max_val, max_val])
        ax.set_zlim([-max_val, max_val])

        ax.set_xlabel(['x (km)'])
        ax.set_ylabel(['Y (km)'])
        ax.set_zlabel(['Z (km)'])

        ax.set_aspect('equal')

        ax.set_title(title)
        plt.legend()

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)
