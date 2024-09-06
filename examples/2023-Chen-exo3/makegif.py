import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib as mpl
from PIL import Image

# Path to your single combined .nc file
filepath = "cart_uranus2j-main.nc"


with Dataset(filepath, mode="r") as nc:
    lat = nc.variables["lat"][:]
    print(lat)
    time = nc.variables["time"][:]
    lon = nc.variables["lon"][:]
    x1 = nc.variables["x1"][:]
    vel1= nc.variables["vel1"][:]
    temp= nc.variables["temp"][:]
    press= nc.variables["press"][:]
    vlat= nc.variables["vlat"][:]
    vlon= nc.variables["vlon"][:]

secperyear = 3.154e+7
steparr = np.array([0.0001,0.001,0.01,0.1,0.1,1,10,100])
a = np.array([[1],[0.5],[0.25],[0.2]])
steparr = steparr*a


def makegif(data,folder,trim,title,short,cmap="RdYlBu",frames=len(time)):
    # Meshgrid of 2D coordinates.
    Lon, Lat = np.meshgrid(lon, lat)
    # Start the figure, make sure it's square and turn off the Axes labels.
    fig, ax = plt.subplots()
    # ax.axis('equal')
    # ax.axis('off')
    vmax = np.percentile(np.mean(data[:frames,:,:,:],axis=1),100-trim)
    vmin = np.percentile(np.mean(data[:frames,:,:,:],axis=1),trim)
    A = np.square(steparr/(vmax-vmin)-0.05)
    n,m = np.unravel_index(A.argmin(), A.shape)
    step = steparr[n,m]
    levels = np.arange(step*int(vmin/step)-step,vmax+step,step)
    global cf
    # Plot the filled and line contours, and the outline of the mask.
    cf = ax.contourf(lon,lat,np.mean(data[0,:,:,:],axis=0), levels,cmap=cmap,extend="both")

    plt.title('Uranus {0} t={1:.2e}y'.format(title,time[0]/secperyear))
    plt.xlabel("Longitude (degrees)")
    plt.ylabel("Latitude (degrees)")
    plt.colorbar(cf,ax=ax,extend="both")

    def animate(i):
        """Set the data for the ith iteration of the animation."""
        global cf
        # Update the plot objects: remove the previous collections to save memory.
        for coll in cf.collections:
            coll.remove()
        cf = ax.contourf(lon,lat,np.mean(data[i,:,:,:],axis=0), levels,cmap=cmap,extend="both")
        plt.title('Uranus {0} t={1:.2e}y'.format(title,time[i]/secperyear))
        return cf

    anim = animation.FuncAnimation(fig, animate, frames=frames, repeat=False)
    anim.save('{}/uranus{}.gif'.format(folder,short), fps=5)

def makeprofgif(data,folder,trim,title,short,cmap="RdYlBu",frames=len(time)):
    # Meshgrid of 2D coordinates.
    Lon, Lat = np.meshgrid(lon, lat)
    # Start the figure, make sure it's square and turn off the Axes labels.
    fig, ax = plt.subplots()
    # ax.axis('equal')
    # ax.axis('off)
    vmax = np.percentile(np.mean(data[:frames,:,:,:],axis=3),100-trim)
    vmin = np.percentile(np.mean(data[:frames,:,:,:],axis=3),trim)
    A = np.square(steparr/(vmax-vmin)-0.05)
    n,m = np.unravel_index(A.argmin(), A.shape)
    step = steparr[n,m]
    levels = np.arange(step*int(vmin/step)-step,vmax+step,step)
    global cf
    # Plot the filled and line contours, and the outline of the mask.
    cf = ax.contourf(lat,(x1-np.min(x1))/1e3,np.mean(data[0,:,:,:],axis=2), levels,cmap=cmap,extend="both")

    plt.title('Uranus {0} Profile t={1:.2e}y'.format(title,time[0]/secperyear))
    plt.xlabel("Latitude (degrees)")
    plt.ylabel("Altitude (km)")
    plt.colorbar(cf,ax=ax,extend="both")

    def animate(i):
        """Set the data for the ith iteration of the animation."""
        global cf
        # Update the plot objects: remove the previous collections to save memory.
        for coll in cf.collections:
            coll.remove()
        cf = ax.contourf(lat,(x1-np.min(x1))/1e3,np.mean(data[i,:,:,:],axis=2), levels,cmap=cmap,extend="both")
        plt.title('Uranus {0} t={1:.2e}y'.format(title,time[i]/secperyear))
        return cf

    anim = animation.FuncAnimation(fig, animate, frames=frames, repeat=False)
    anim.save('{}/uranus{}prof.gif'.format(folder,short), fps=5)


#######
folder = "gifs"
if not(os.path.isdir(folder)):
    os.makedirs(folder)
makegif(np.log10(press),folder,2,"Pressure (log(Pa))","2jpress",frames=32)
makeprofgif(np.log10(press),folder,2,"Pressure Profile (log(Pa))","2jpress",frames=32)
# print(time)
# print(vlat[:32,:,:,:])
#use seismic colour map for temperature