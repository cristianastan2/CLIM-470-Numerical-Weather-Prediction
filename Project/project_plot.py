"""
** 1 Nov. 2026 release
Features added
1.  reads "Hs" files for all 3 resolutions - HsNNN.dat
2.  produces a topo map for each resolution
3.  produces a single topo section including all 3 resolutions --
    they overlay each other, but red and green markers help.  see legend
4.  added minor grid lines to the topo section
5.  added datetime to filenames for png plot files
"""

from datetime import datetime
from os import WCONTINUED

import numpy as np
from array import array
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MultipleLocator


"""ds_u_bin = open(dirin+uDfil,'rb')
for t in range(ntime):
    tmp = ds_u_bin.read(recl)
    tmp1 = array('f',tmp)
    #tmp2=tmp1[1:-1] #for fortran sequential
    u_bin[t,:,:]=np.reshape(tmp1,(nlat,nlon))
del ds_u_bin,tmp,tmp1

ds_v_bin = open(dirin+vDfil,'rb')
for t in range(ntime):
    tmp = ds_v_bin.read(recl)
    tmp1 = array('f',tmp)
    #tmp2=tmp1[1:-1] #for fortran sequential
    v_bin[t,:,:]=np.reshape(tmp1,(nlat,nlon))
del ds_v_bin,tmp,tmp1 """

# calculate time mean
#z_bin_mean = z_bin.mean(axis=0)
#u_bin_mean = u_bin.mean(axis=0)
#v_bin_mean = v_bin.mean(axis=0)

#np.shape(z_bin_mean)

def plot_topo(z_bin, grid_cell_size, run_date):
    # plot one time

    # zc = ax.contour(z_bin[:,:])
    # ax.set_title('project surface topography')

    x = np.arange(0,6010, grid_cell_size)
    y = np.arange(0,2010, grid_cell_size)
    fig, ax = plt.subplots(figsize=(8, 5), dpi=150)

    cmap = LinearSegmentedColormap.from_list('topo', ['lightblue', 'blue', 'green', 'yellow'])
    zc = ax.contourf(x, y, z_bin, cmap=cmap)
    # zc = ax.contourf(x,y,z_bin[:,:], cmap='terrain')
    ax.set_title('project topography', pad=22)
    ax.text(0.5, 1.02, "Gary Greenwood (CLIM470)", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=10, color="gray")
    ax.text(0.5, -0.5, "showing topography generated with cell-size: "+ str(grid_cell_size),
            transform=ax.transAxes, ha="center", va="bottom", fontsize=10, color="black")
    fig.colorbar(zc, ax=ax, label='elevation')
    ax.set_aspect('equal')
    ax.grid(True, color='black', alpha=0.4, linestyle='--', linewidth=0.6)
    ax.set_ylabel("km in y direction (latitude)", color="red")
    ax.set_xlabel("km in x direction (longitude)", color="red")

    fig.savefig("./firsttopo" + str(grid_cell_size) + "_" + run_date +".png", bbox_inches="tight")
    plt.clf()   # clear current figure
    plt.cla()   # clear current axes
    plt.close(fig)
    return (1)

def plot_section(z_bin, sizes, hi_resolution, run_date):
    #
    #  this plot_section will overlay values for all resolutions, as indicated in sizes.
    #  z_bin is a dict containing arrays for each resolution
    x = np.arange(0,6010, hi_resolution)    #x-axis will be for highest resolution

    # z = z_bin[0, :]
    z1 = z_bin[sizes[0]][1,:]
    z2 = z_bin[sizes[1]][1,:]
    nlon2 = z2.shape[0]
    x2 = np.arange(nlon2) * sizes[1]
    print(f'x2 number: ', x2)
    z3 = z_bin[sizes[2]][1,:]
    nlon3 = z3.shape[0]
    x3 = np.arange(nlon3) * sizes[2]

    fig, ax = plt.subplots(figsize=(8, 5), dpi=150)
    ax.plot(x, z1, color='blue', linewidth=1, label='125km section',
            marker = "D", markersize = 3, markerfacecolor='white')
    #         color="green",
    #         linestyle="--",
    #         marker="o",
    ax.plot(x2, z2, color='red', linewidth=1, label='250km section',
            marker="o", markersize = 3)
    ax.plot(x3, z3, color='green', linewidth=1, label='500km section',
            linestyle="--")
    # marker values here >> https://matplotlib.org/stable/api/markers_api.html


    # fig.colorbar(zc, ax=ax, label='elevation')
    ax.set_aspect('equal')
    # set up grid lines
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_minor_locator(MultipleLocator(250))
    ax.yaxis.set_major_locator(MultipleLocator(500))
    ax.yaxis.set_minor_locator(MultipleLocator(125))
    ax.grid(True, which='major', color='black', alpha=0.4, linestyle='--', linewidth=0.8)
    ax.grid(True, which='minor', color='black', alpha=0.2, linestyle=':', linewidth=0.6)

    # title and text
    ax.set_title('Project Topography - Section', pad=22)
    ax.text(0.5, 1.02, "Gary Greenwood (CLIM470)", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=10, color="gray")
    ax.text(0.5, -0.5, "showing topo section generated for all grid cell sizes "
            + str(sizes[0]) + ", " + str(sizes[1]) + ", " + str(sizes[2]),
            transform=ax.transAxes, ha="center", va="bottom", fontsize=10, color="black")

    # set up labels and legend
    ax.set_ylabel("z height in meters", color="red")
    ax.set_xlabel("km in x direction (longitude)", color="red")
    ax.get_legend_handles_labels()
    ax.legend(loc = "upper right")

    fig.savefig("./alltoposection_" + run_date +".png", bbox_inches="tight")
    plt.clf()   # clear current figure
    plt.cla()   # clear current axes
    plt.close(fig)
    return(1)

def get_data(sizes, resolution, reclength):
    # no time steps, single variable...
    # ntime=20
    sizes = (125, 250, 500)    # grid cell sizes.
    values = {}
    z_bin = {}

    for i in sizes:
        # dwith open(f"HS_{i}") as f:
        #    values[i] = f.read()   # or however you parse the file

        #heights[i] = getdata(f"height_{i}.dat")

        # first calculate recl based on resolution
        grid_cell_size=i     # will be  125, 250, 500
        nlon=int(6000/grid_cell_size)+1     # +1 includes both edges
        nlat=int(2000/grid_cell_size)+1     # +1 includes both edges
        recl=(nlon*nlat)*4                  # length in bytes, each data point is 4-bytes
        z_bin[i]=np.zeros((nlat,nlon))

        file_open = False

        try:
            newfname = f"Hs{i}.dat"
            ds_z_bin = open(newfname, "rb")
            # ds_z_bin = open(f"Hs{i}.dat", "rb")
        except FileNotFoundError:
            print("The file: ", newfname, " does not exist.")
            return(97, None)
        except Exception:
            print("Other file open error: ", Exception)
            return(99, None)          # don't know what to do at this point
        else:
            file_open = True
            print("File opened successfully.")

        try:
            values[i] = ds_z_bin.read(recl)
        except Exception:
            print("File read error: ", Exception)
            return (98, None)
        z_bin[i] =  np.frombuffer(values[i], dtype=np.float32).reshape(nlat, nlon)
    # z_bin is a dictionary containing numpy arrays.
    #  now create arrays based on grid size.
    #tmp125 = array('f',values[125])
    #tmp250 = array('f', values[250])
    #tmp500 = array('f', values[500])
    #tmp2=tmp1[1:-1] #for fortran sequential
    #z_bin[:,:]=np.reshape(tmp125,(nlat,nlon))

    return (0, values, z_bin)

    # ds_z_bin = open(dirin+zDfil,'rb')
    # for t in range(ntime):


def otherplots(z_bin_mean, nlon, nlat, grid_cell_size):
    """
    This function is included in the project_plt.py to provide examples
    for other ways to plot values
    :param z_bin_mean:
    :param nlon:
    :param nlat:
    :param grid_cell_size:
    :return:
    """
    x = np.arange(0,6010, grid_cell_size)
    y = np.arange(0,2010, grid_cell_size)
    # plot filled contours - for later
    fig, ax = plt.subplots(figsize=(10,4))
    zc = ax.contourf(x, y, z_bin_mean,cmap=cm.coolwarm)
    plt.colorbar(zc, location= 'bottom')

    # plot the contour lines using gray lines
    zl = ax.contour(x, y, z_bin_mean,cmap=cm.gray)
    ax.set_title('gh time mean')


    # Plot gh and wind vectors
    fig, ax = plt.subplots(figsize=(10,3))
    zc = ax.contour(z_bin_mean[36:72,:])
    skip=(slice(None, None, 3), slice(None, None, 3))
    x=np.arange(nlon)
    y=np.arange(nlat)
    # zc = ax.quiver(x[::3],y[0:36:3],u_bin_mean[36:72,:][skip], v_bin_mean[36:72,:][skip],color='grey',
    #               headwidth=2,headlength=2)

    # Plot gh at all longitudes and a random latitude in the Northern Hemisphere
    lons=np.arange(0,360,2.5)
    lats=np.arange(-90,90,2.5)
    fig, ax = plt.subplots(figsize=(10,3))
    z_lat1 = ax.plot(lons, z_bin_mean [65,:],color='b',label=str(lats [65])+' N')
    z_lat2 = ax.plot(lons, z_bin_mean [70,:],color='r',label=str(lats [70])+' N')
    ax.set_title('gh at'+str(lats [65])+' N and '+str (lats [70])+' N')
    ax.set_ylabel('geopotential height, m')
    ax.set_xlabel('longitude')
    plt.legend ()
    return()

def main():
    print(f'starting...')

    now_utc = datetime.now()

    # Get the current UTC time, properly timezone-aware
    run_date =  now_utc.strftime("%Y%m%d_%H%M")

    # we will get the resolution for plotting, but graphing will use all 3
    # resolution = input("enter 500, 250 or 125 to select resolution to plot")
    # ERROR CHECKING needed to verify entered data is correct
    """while resolution.strip() not in ('500', '250', '125'):
        print(f'resolution entered = ', resolution)
        resolution = input("enter 500, 250 or 125 to select resolution to plot")

    for i in (125, 250, 500):
        res_string = str(i)
        dirin='./'
        zDfil='Hs' +res_string + '.dat'      # default filename
        print(f'filename with resolution =', zDfil)

        # these may be used later
        uDfil='uwnd.dat'
        vDfil='vwnd.dat'

        # no time steps, single variable...
        # ntime=20
        grid_cell_size=i     # woll be  125, 250, 500
        nlon=int(6000/grid_cell_size)+1     # +1 includes both edges
        nlat=int(2000/grid_cell_size)+1     # +1 includes both edges
        recl=(nlon*nlat)*4
        z_bin=np.zeros((nlat,nlon))
        #u_bin=np.zeros((ntime,nlat,nlon))
        #v_bin=np.zeros((ntime,nlat,nlon))

        # ds_z_bin = open(dirin+zDfil,'rb')
        # for t in range(ntime):
        """
    sizes = (125, 250, 500)    # grid cell sizes.
    rc, values, z_bin = get_data(sizes, None, None)
    if rc != 0 : exit (rc)
    # need to recalc nlon/nlat to get new recl...
    # rc, tmp_250 = get_data(zDfil, "250", recl)

    #tmp1 = array('f',tmp)
    #tmp2=tmp1[1:-1] #for fortran sequential
    #z_bin[:,:]=np.reshape(tmp1,(nlat,nlon))
    #del ds_z_bin,tmp,tmp1

    for j in sizes:
        topo_fig = plot_topo(z_bin[j], j, run_date)         # single grid size, but loop thru 3
    section_fig = plot_section(z_bin, sizes, 125, run_date)   # arrays for each of 3 gridsizes
    exit(0)

# *************************
# Script Execution (Bottom)
# --------------------------

if __name__ == "__main__":
    print("No args, calling main")
    main()
