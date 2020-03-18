#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Descriptif du fichier
"""

# Importation des librairies
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import numpy as np
import scipy.integrate as scp
import matplotlib.pyplot as plt

from astropy.time import Time
import astropy.coordinates as asc
import astropy.units as u
# Definition des fonctions

class point:

    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
        self.all = np.array([x,y,z])

    def modif(self,array):
        """Attribuer les nouvelles coordonnees contenues dans array a self"""
        self.x = array[0]
        self.y = array[1]
        self.z = array[2]
        self.all = array

class SolarSystem:

    def __init__(self,name,position,velocity,mass):
        """Initiation des objets de classe"""
        self.name = name
        self.position = position
        self.velocity = velocity
        self.mass = mass

    def norm(self,other):
        """Norme de du vecteur (self_other)"""
        return(np.sqrt((other.position.x - self.position.x)**2 + (other.position.y - self.position.y)**2 + (other.position.z - self.position.z)**2))

    def acceleration(self,other,G):
        """Calcule les composantes de l'acceleration de self """
        acceleration = point(0,0,0)
        if len(other) > 1:
            for planet in other:
                # print(planet.name)
                if planet.name != self.name:
                    norm=SolarSystem.norm(self,planet)
                    acceleration.x += G * planet.mass * ( planet.position.x - self.position.x) / norm**3
                    acceleration.y += G * planet.mass * ( planet.position.y - self.position.y) / norm**3
                    acceleration.z += G * planet.mass * ( planet.position.z - self.position.z) / norm**3
                    acceleration.modif([acceleration.x,acceleration.y,acceleration.z])
        else:
            norm = SolarSystem.norm(self,other)
            acceleration.x += G * other.mass * ( other.positon.x - self.position.x) / norm**3
            acceleration.y += G * other.mass * ( other.position.y - self.position.y) / norm**3
            acceleration.z += G * other.mass * ( other.position.z - self.position.z) / norm**3
            acceleration.modif([acceleration.x,acceleration.y,acceleration.z])
        return acceleration

def deriv(X,t,system,G):
    """Calcule la derivee X (position et vitesse) au temps t"""
    X_copy = np.copy(X).reshape(len(system),6)
    k=0
    for body in system:
        body.position.modif(X_copy[k,:3])
        body.velocity.modif(X_copy[k,3:]) ; k +=1
    X_new = np.zeros(len(X))
    X_new[::6]=X[3::6]
    X_new[1::6]=X[4::6]
    X_new[2::6]=X[5::6]
    k=0
    for body in system:
        X_new[6*k+3:6*k+6] = body.acceleration(system,G).all; k += 1
    return X_new


# Programme principal
if __name__ == "__main__":

    t = Time("2020-02-10 23:00")
    jup = asc.get_body_barycentric_posvel('Earth', t)
    pos=np.array(jup[0].xyz.to(u.m))
    vel=np.sum(jup[1].xyz.to(u.m/u.s).value**2)

    AU = 149597870700
    d = 24*3600
    run = False
    save = False
    load = True
    plot = True
    t_sim = 1e9
    n_points = int(1e4)


    #Constantes
    G = 6.67e-11

    #Corps du systeme
    Sun = SolarSystem("Sun",point(0,0,0),point(0,0,0),1.989e30)
    Mercury = SolarSystem("Mercury",point(56e9,0,0),point(0,0,0),3.3e23)
    Venus = SolarSystem("Venus",point(108e9,0,0),point(0,0,0),4.9e24)
    Earth = SolarSystem("Earth",point(149.6e9,0,0),point(0,0,0),5.972e24)
    Mars = SolarSystem("Mars",point(228e9,0,0),point(0,0,0),6.4e23)
    Jupiter = SolarSystem("Jupiter",point(778.5e9,0,0),point(0,0,0),1.89e27)
    Saturn = SolarSystem("Saturn",point(1427e9,0,0),point(0,0,0),5.7e26)
    Uranus = SolarSystem("Uranus",point(2871e9,0,0),point(0,0,0),8.7e25)
    Neptune = SolarSystem("Neptune",point(4498e9,0,0),point(0,0,0),1e26)

    system = [Sun,Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune]
    planets = system[1:]                #J'ai ajouté ça pour faire quelques précisions dans la suite


    for body in system:
        pos_vel = asc.get_body_barycentric_posvel(body.name, t)
        pos=pos_vel[0].xyz.to(u.m).value
        vel=pos_vel[1].xyz.to(u.m/u.s).value
        body.position.modif(pos_vel[0].xyz.to(u.m).value)
        body.velocity.modif(pos_vel[1].xyz.to(u.m/u.s).value)


    #Pour calculer plus précisément les vitesses initiales (la traj de la Terre avait une gueule bizarre avant, tu avais mis 35000 m/s à la place d'un peu moins de 30000 m/s, mais celle de Jupiter était assez correcte
    #for pl in planets :
    #    pl.velocity.modif([0, np.sqrt(G*Sun.mass/pl.position.x), 0])

    #Pour calculer les conditions initiales du Soleil afin d'avoir un mouvement stable, et des traj concentriques
    #Sun.position.modif([-np.sum([pl.position.x*pl.mass/Sun.mass for pl in planets]), 0, 0])
    #Sun.velocity.modif([0, -np.sum([pl.velocity.y*pl.mass/Sun.mass for pl in planets]), 0])

    long=len(system)

    X_init=np.zeros((long,6))
    
    k=0
    for body in system:
        X_init[k] = np.concatenate((body.position.all,body.velocity.all)) ; k += 1

    X_init=X_init.reshape(long*6)

    time=np.linspace(0,t_sim,n_points)

    if run == True:
        sol,dic=scp.odeint(deriv,X_init,time,args=(system,G),full_output = True)
        if save == True:
            np.savetxt('system_{:.1e}_{:.1e}.csv'.format(t_sim,n_points),sol,delimiter=';')

    if load == True:

        sol = np.loadtxt('system_{:.1e}_{:.1e}.csv'.format(t_sim,n_points),delimiter=';')
        #sol = np.loadtxt('sol.csv',delimiter = ';')
    if plot == True:


        fig=plt.figure(figsize=(12,6))
        ax=fig.subplots(nrows=1,ncols=3)
        k=0
        for body in system:
            ax[0].plot(sol[:,6*k],sol[:,6*k+1],label=body.name)
            #ax[0].legend(loc="right")
            #ax[0].set_aspect('equal', 'box')
            ax[1].plot(sol[:,6*k],sol[:,6*k+2],label=body.name)
            #ax[1].set_aspect('equal', 'box')
            #ax[1].legend(loc="right")
            ax[2].plot(sol[:,6*k+1],sol[:,6*k+2],label=body.name)
            #ax[2].set_aspect('equal', 'box')
            ax[2].legend(loc="right")
            k += 1
        #ax.set_aspect('equal')
        plt.tight_layout()
        plt.show()
        plt.show()


        #Figure
        fig = plt.figure()
        ax = fig.add_subplot((111), projection='3d')
        #ax = fig.gca(projection='3d')
        k=0
        for body in system:
            print(body.name)
            #plt.plot(sol[:,6*k],sol[:,6*k+1],label=body.name)
            ax.plot(sol[:,6*k],sol[:,6*k+1],sol[:,6*k+2],'+',label=body.name)
            k += 1
        #ax.set_aspect('equal')
        plt.legend(loc="best")
        plt.tight_layout()
        plt.show()
