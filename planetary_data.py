#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 20:22:56 2020

@author: astronian
"""
from numpy import array

G_meters = 6.67408e-11
G = G_meters*10**-9

day2sec =  24*3600.0

sun = {
       'name':'Sun',
       'mass':1.989e24,
       'mu':1.32712e11,
       'radius':695700.0
       }

atm= array([[63.096,2.059e-4],[251.189,5.909e-11],[1000.0,3.561e-15]])
earth = {
        'name':'Earth',
        'mass':5.972e24,
        'mu':5.972e24*G,
        'radius':6378.0,
        'J2' : -1.082635854e-3,
        'spice_file':'..\\..\\spice_data\\solar_system\\solar_system.mk',
        'zs':atm[:,0],
        'rhos':atm[:,1]*10**8,
        'atm_rot_vector' :array([0.0,0.0,72.9211e-6])
        }