#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 16:16:31 2022

@author: alecchurch
"""
# %% Imports
import time
import numpy as np
import matplotlib.pyplot as plt
import requests
from random import randrange

# %% CAD API Call Function
def API(CADparams):
    url_CAD =  "https://ssd-api.jpl.nasa.gov/cad.api"
    out = requests.get(url_CAD,params=CADparams).json()
    # Edge case catch, no returned items from API "count = 0"
    if out['count'] == '0':
        return 0,0
    return out['data']
    #return out

# %% Function to call API Function and process data
def processData():
    start = time.time()
    # Constants & Preallocation
    bodies = ['Merc','Venus','Earth','Mars']
    grav = [3.7, 8.9, 9.8, 3.7]
    densities = [5429, 5243, 5514, 3932]
    v_esc = [4.3, 10.4, 11.2, 5.0]
    Dist = []
    V = []
    dia = []
    D_c = []
    p = 0.1
    
    # For each body, call API data, for each asteroid calc & save distance, diameter, & velocity.
    for i,body in enumerate(bodies):
        
        # API request parameters
        CADparams = {
        'date-min':'2002-12-01',
        'limit':100000,
        'body':body,
        }
        
        # NASA JPL API Function Call
        data = API(CADparams)
        #data = out['data']
        
        # Preallocate lists
        diameter = []
        dist = []
        v = []
        d_c = []
        
        # For each asteroid, if it has a measured H_mag, calculate diameter and save both distance and velocity.
        for asteroid in data:
            if type(asteroid[10]) != type(None):
                crater_d = crater(float(asteroid[7]),grav[i],densities[i],((1329/np.sqrt(p))*10**(-0.2*float(asteroid[-1]))),v_esc[i])
                diametr = (1329/np.sqrt(p))*10**(-0.2*float(asteroid[-1]))
                if diametr > 1:
                    v.append(float(asteroid[7]))
                    diameter.append(diametr)
                    d_c.append(crater_d)
        
        # Save all of the asteroid data into a sublist that relates to the body it passes by.
        V.append(v)
        Dist.append(dist)
        dia.append(diameter)
        D_c.append(d_c)
    
    # Calculate and print runtime of API call and data saving section of main function.
    endAPI = time.time()
    print("API Program execution time:", (endAPI-start), "s.")
    return V,Dist,dia,D_c

# %% Crater diameter calculation
def crater(v0,g,rho_g,d,v_esc):
    # Average asteroid density [kg/m^3]
    rho_i = 2000
    # Impact Velocity [km/s]
    v_i = np.sqrt(v0**2 + v_esc**2)
    # Initial velocity, gravity, and diameter to calculate
    D_c = (v_i**2/g)**(1/4) * (rho_g/rho_i)**(1/4) * d**(3/4)
    return D_c

# %% Size-Frequency Distribution Constants Calculation
def constants(d_c):
    b = 3
    time = 20
    q = b+1
    C_i = -b * (len(d_c)/time) / (np.max(d_c)**(-b) - np.min(d_c)**(-b))
    C_c = C_i / b
    return b,q,C_i,C_c

# %% SFD Function
def sfd(d,C_i,q):
    N_i = C_i * d**(-q)
    return N_i

# %% SFD Function for log plotting
def sfd2(D,C_i):
    b = 3
    q = b + 1
    N = C_i*D**(-q)
    return N

# %% Cratering Simulation Function
def simulation(T,craters,d_c):
    start = time.time()
    partial = False
    crater_points = []
    area_covered = 0
    
    # Calc constants
    b,q,C_i,C_c = constants(d_c)
    print('C_i = {}'.format(C_i))
    # Find formation frequency at different crater sizes
    dt_impact = []
    for d in craters:
        dt_impact.append([sfd(d,C_i,q),d])
    impact_events = np.ones(len(dt_impact))
    crater_count = 0
    
    # Simulation
    for t in T:
        print(t)
        for i,Nd in enumerate(dt_impact):
            if (Nd[0]*t) >= impact_events[i]:
                impact_events[i] += 1
                crater_count += 1
                d = Nd[1]
                #print('At time {} [years] from simulation start, a crater of diameter {} [km] will impact.'.format(t,d))
                # Random crater impact point
                x0 = randrange(1,500,1)
                y0 = randrange(1,500,1)
                r = int(d/2)
                del_ind = []
                # Create circle equation:
                theta = np.linspace(0,2*np.pi,100)
                x = r * np.cos(theta) + x0
                y = r * np.sin(theta) + y0
                area_covered += np.pi * r**2
                # Iterate through all old craters
                for k,cp in enumerate(crater_points): 
                    inside_index = []
                    distance = np.sqrt(((x0-cp[2])**2 + (y0-cp[3])**2))
                    if (cp[4]+r) > distance:
                        if partial == True:
                            for i in range(len(cp[0])):
                                p_dist = np.sqrt((x0-cp[0][i])**2+(y0-cp[1][i])**2) # Distance from new center to all points on old crater
                                if r >= p_dist: # If the distance to old point is less than new radius, the crater will destroy that point
                                    inside_index.append(i)
                            replace_x = np.delete(cp[0],inside_index)
                            replace_y = np.delete(cp[1],inside_index)
                            crater_points[k][0] = replace_x
                            crater_points[k][1] = replace_y
                        elif partial == False:
                            del_ind.append(k)
                for ind in sorted(del_ind, reverse=True):
                    del crater_points[ind]
                crater_points.append([x,y,x0,y0,r])
    endT = time.time()
    print('Simulation call time to complete: {} s'.format(endT-start))
    print('Beginning plotting...')
    plotting(crater_points)
    return

# %% Plotting Function
def plotting(crater_points):
    figure,axes=plt.subplots()
    axes.set_aspect(1)
    axes.set(xlim=(0,500),ylim=(0,500))
    for cp in crater_points:
        circ_x = cp[0]
        circ_y = cp[1]
        axes.scatter(circ_x,circ_y,s=.06,c='black',marker='.')
    plt.show()
    print('Plotting complete...')
    return

# %% Main function
def Main():
    bodies = ['Merc','Venus','Earth','Mars']
    
    Nm = []
    Nv =[]
    Ne=[]
    Nma =[]
    craters = np.linspace(10,200,96)
    for c in craters:
        Nm.append(sfd2(c,8.984e3))
        Nv.append(sfd2(c,6.475e3))
        Ne.append(sfd2(c,6.157e3))
        Nma.append(sfd2(c,2.33e3))
    
    plt.figure(0)
    plt.loglog(craters,Nm,'red',label='Mercury')
    plt.loglog(craters,Nv,'orange',label='Venus')
    plt.loglog(craters,Ne,'green',label='Earth')
    plt.loglog(craters,Nma,'blue',label='Mars')
    plt.title('Logarithmic Size-Frequency Distribution for Planets')
    plt.xlabel(r'Crater Diameter [$km$]')
    plt.ylabel(r'Frequency [$Years^{-1}$]')
    plt.legend()
    plt.show()
    
    # Call Data Acquisition & Processing Function
    V,Dist,dia,D_c = processData()
    print('API call complete...')
    
    # All data from each planetary body has been collected into V, Dist, dia, & crater diameter for each asteroid.
    T = np.linspace(1000,10**6,500)
    craters = np.linspace(1,200,96)
    for i,d_c in enumerate(D_c):
        print('Calling simulation for body: {}'.format(bodies[i]))
        simulation(T,craters,d_c)
        print('Simulation and plotting complete for lifetime cratering on body: {}'.format(bodies[i]))


# %% Main Function Begin
if __name__ == '__main__':
    Main()



