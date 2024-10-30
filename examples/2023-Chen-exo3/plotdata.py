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
filepath = "cart_uranus2p-main.nc"


with Dataset(filepath, mode="r") as nc:
    lat = nc.variables["lat"][:]
    time = nc.variables["time"][:]
    lon = nc.variables["lon"][:]
    vel1= nc.variables["vel1"][:]
    temp= nc.variables["temp"][:]
    press= nc.variables["press"][:]
    vlat= nc.variables["vlat"][:]
    vlon= nc.variables["vlon"][:]
    x1 = nc.variables["x1"][:]
    forc = nc.variables["forcing"][:]
    goal = nc.variables["goalprof"][:]


def makeprofgif(data,folder,title,short,x,y,avgaxis,xlabel,ylabel,frames=len(time),trim=1,cmap="RdYlBu"):
    steparr = np.array([0.0001,0.001,0.01,0.1,0.1,1,10,100])
    a = np.array([[1],[0.5],[0.25],[0.2]])
    steparr = steparr*a
    secperyear = 3.154e+7

    # Meshgrid of 2D coordinates.
    # Start the figure, make sure it's square and turn off the Axes labels.
    fig, ax = plt.subplots()
    # ax.axis('equal')
    # ax.axis('off)
    vmax = np.percentile(np.mean(data[:,:,:,:],axis=avgaxis+1),100-trim)
    vmin = np.percentile(np.mean(data[:,:,:,:],axis=avgaxis+1),trim)
    print(vmin)
    print(vmax)
    A = np.square(steparr/(vmax-vmin)-0.05)
    n,m = np.unravel_index(A.argmin(), A.shape)
    step = steparr[n,m]
    levels = np.arange(step*int(vmin/step)-step,vmax+step,step)
    global cf
    # Plot the filled and line contours, and the outline of the mask.
    cf = ax.contourf(x,y,np.mean(data[0,:,:,:],axis=avgaxis), levels,cmap=cmap,extend="both")

    plt.title('Uranus {0} Profile t={1:.2e}y'.format(title,time[0]/secperyear))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.colorbar(cf,ax=ax,extend="both")

    def animate(i):
        """Set the data for the ith iteration of the animation."""
        global cf
        # Update the plot objects: remove the previous collections to save memory.
        for coll in cf.collections:
            coll.remove()
        cf = ax.contourf(x,y,np.mean(data[i,:,:,:],axis=avgaxis), levels,cmap=cmap,extend="both")
        plt.title('Uranus {0} t={1:.2e}y'.format(title,time[i]/secperyear))
        return cf

    anim = animation.FuncAnimation(fig, animate, frames=frames, repeat=False)
    anim.save('{}/uranus{}prof.gif'.format(folder,short), fps=10)

# def makeprofgif(data,folder,trim,title,short,cmap="RdYlBu",frames=len(time)):
#     # Meshgrid of 2D coordinates.
#     # Start the figure, make sure it's square and turn off the Axes labels.
#     fig, ax = plt.subplots()
#     # ax.axis('equal')
#     # ax.axis('off)
#     vmax = np.percentile(np.mean(data[:frames,:,:,:],axis=3),100-trim)
#     vmin = np.percentile(np.mean(data[:frames,:,:,:],axis=3),0)
#     A = np.square(steparr/(vmax-vmin)-0.05)
#     n,m = np.unravel_index(A.argmin(), A.shape)
#     step = steparr[n,m]
#     levels = np.arange(step*int(vmin/step)-step,vmax+step,step)
#     global cf
#     # Plot the filled and line contours, and the outline of the mask.
#     cf = ax.contourf(lat,press,np.mean(data[0,:,:,:],axis=2), levels,cmap=cmap,extend="both")

#     plt.title('Uranus {0} Profile t={1:.2e}y'.format(title,time[7]/secperyear))
#     plt.xlabel("Latitude (degrees)")
#     plt.ylabel("Pressure (Pa)")
#     plt.yscale("log")
#     plt.colorbar(cf,ax=ax,extend="top")
#     plt.gca().invert_yaxis()
#     plt.savefig('{}/uranus{}prof.png'.format(folder,short))


#######
folder = "gifs4"
if not(os.path.isdir(folder)):
    os.makedirs(folder)


#plot temperature
data= np.copy(temp[0:147,:,:,:])
timerange = [0,147]
heightrange=[0,len(x1)]
latrange = [0,len(lat)]
lonrange = [0,len(lon)]
avgaxis= 2 #0 is height, 1 is lat, 2 is lon
title = "Temperature"
shorttitle="2ptemp"
setframes = len(data[timerange[0]:timerange[1],0,0,0])
x = lat
y = (x1-min(x1))/1e3
xlabel = "Latitude (deg)"
ylabel = "Height (km)"
makeprofgif(data[timerange[0]:timerange[1],heightrange[0]:heightrange[1],latrange[0]:latrange[1],lonrange[0]:lonrange[1]],folder,title,shorttitle,x,y,avgaxis,xlabel,ylabel,frames=setframes,trim=0,cmap="seismic")
#use seismic colour map for temperature

#plot forcing
data= np.copy(forc[0:147,:,:,:])
timerange = [0,147]
heightrange=[0,len(x1)]
latrange = [0,len(lat)]
lonrange = [0,len(lon)]
avgaxis= 2 #0 is height, 1 is lat, 2 is lon
title = "Forcing"
shorttitle="2pforc"
setframes = len(data[timerange[0]:timerange[1],0,0,0])
x = lat
y = (x1-min(x1))/1e3
xlabel = "Latitude (deg)"
ylabel = "Height (km)"
makeprofgif(data[timerange[0]:timerange[1],heightrange[0]:heightrange[1],latrange[0]:latrange[1],lonrange[0]:lonrange[1]],folder,title,shorttitle,x,y,avgaxis,xlabel,ylabel,frames=setframes,trim=0,cmap="seismic")
#use seismic colour map for temperature

#plot goal profile
data= np.copy(goal[0:147,:,:,:])
timerange = [0,147]
heightrange=[0,len(x1)]
latrange = [0,len(lat)]
lonrange = [0,len(lon)]
avgaxis= 2 #0 is height, 1 is lat, 2 is lon
title = "Goal Profile"
shorttitle="2pgoal"
setframes = len(data[timerange[0]:timerange[1],0,0,0])
x = lat
y = (x1-min(x1))/1e3
xlabel = "Latitude (deg)"
ylabel = "Height (km)"
makeprofgif(data[timerange[0]:timerange[1],heightrange[0]:heightrange[1],latrange[0]:latrange[1],lonrange[0]:lonrange[1]],folder,title,shorttitle,x,y,avgaxis,xlabel,ylabel,frames=setframes,trim=0,cmap="seismic")
#use seismic colour map for temperature

#plot zonal winds
data= np.copy(vel1[0:147,:,:,:])
timerange = [0,147]
heightrange=[0,len(x1)]
latrange = [0,len(lat)]
lonrange = [0,len(lon)]
avgaxis= 2 #0 is height, 1 is lat, 2 is lon
title = "Zonal Wind Profile"
shorttitle="2pzw"
setframes = len(data[timerange[0]:timerange[1],0,0,0])
x = lat
y = (x1-min(x1))/1e3
xlabel = "Latitude (deg)"
ylabel = "Height (km)"
makeprofgif(data[timerange[0]:timerange[1],heightrange[0]:heightrange[1],latrange[0]:latrange[1],lonrange[0]:lonrange[1]],folder,title,shorttitle,x,y,avgaxis,xlabel,ylabel,frames=setframes,trim=0)
