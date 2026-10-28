import numpy as np
from array import array
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap


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

def plot_topo(z_bin, grid_cell_size):
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
    ax.grid(True, color='black', alpha=0.4, linestyle='--', linewidth=0.5)
    ax.set_ylabel("km in y direction (latitude)", color="red")
    ax.set_xlabel("km in x direction (longitude)", color="red")

    fig.savefig("./firsttopo" + str(grid_cell_size) + ".png", bbox_inches="tight")
    plt.clf()   # clear current figure
    plt.cla()   # clear current axes
    plt.close(fig)
    return (1)

def plot_section(z_bin, grid_cell_size):
    x = np.arange(0,6010, grid_cell_size)
    z = z_bin[0, :]
    fig, ax = plt.subplots(figsize=(8, 5), dpi=150)
    ax.plot(x, z, color='blue', linewidth=2, label='topo section')

    ax.set_title('Project Topography - Section', pad=22)
    ax.text(0.5, 1.02, "Gary Greenwood (CLIM470)", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=10, color="gray")
    ax.text(0.5, -0.5, "showing topo section generated with cell-size " + str(grid_cell_size),
            transform=ax.transAxes, ha="center", va="bottom", fontsize=10, color="black")
    # fig.colorbar(zc, ax=ax, label='elevation')
    ax.set_aspect('equal')
    ax.grid(True, color='black', alpha=0.4, linestyle='--', linewidth=0.5)
    ax.set_ylabel("z height in meters", color="red")
    ax.set_xlabel("km in x direction (longitude)", color="red")

    fig.savefig("./firsttoposection" + str(grid_cell_size) + ".png", bbox_inches="tight")
    plt.clf()   # clear current figure
    plt.cla()   # clear current axes
    plt.close(fig)
    return(1)

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
    resolution = input("enter 500, 250 or 125 to select resolution to plot")
    # ERROR CHECKING needed to verify entered data is correct

    dirin='./'
    zDfil='Hs' + resolution.strip() + '.dat'
    print(f'filename with resolution =', zDfil)
    uDfil='uwnd.dat'
    vDfil='vwnd.dat'

    # no time steps, single variable...
    # ntime=20
    grid_cell_size=int(resolution)      # could also be 500, 250
    nlon=int(6000/grid_cell_size)+1     # +1 includes both edges
    nlat=int(2000/grid_cell_size)+1     # +1 includes both edges
    recl=(nlon*nlat)*4
    z_bin=np.zeros((nlat,nlon))
    #u_bin=np.zeros((ntime,nlat,nlon))
    #v_bin=np.zeros((ntime,nlat,nlon))

    ds_z_bin = open(dirin+zDfil,'rb')
    # for t in range(ntime):
    tmp = ds_z_bin.read(recl)
    tmp1 = array('f',tmp)
    #tmp2=tmp1[1:-1] #for fortran sequential
    z_bin[:,:]=np.reshape(tmp1,(nlat,nlon))
    #del ds_z_bin,tmp,tmp1

    topo_fig = plot_topo(z_bin, grid_cell_size)
    section_fig = plot_section(z_bin, grid_cell_size)
    exit(2)

# *************************
# Script Execution (Bottom)
# --------------------------

if __name__ == "__main__":
    print("No args, calling main")
    main()
