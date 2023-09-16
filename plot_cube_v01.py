import numpy as np
from ase.io.cube import read_cube_data
import matplotlib.pyplot as plt

importfolder3A= '/scratch/hpc-prf-eprop2d/eprop2d1_Jan/11_MIR_plasmonics/01_50_kV_m_30_THz_1010_atoms_20_cycles/01_dist_3A/' # data folder
importfolder10A= '/scratch/hpc-prf-eprop2d/eprop2d1_Jan/11_MIR_plasmonics/01_50_kV_m_30_THz_1010_atoms_20_cycles/04_dist_10A/' # data folder

timesteps = 20 # number of timesteps
time_array = np.arange(timesteps) * 0.2

## import .cube file for reference electron density
cube_data_ref_3A, cube_atoms_ref_3A = read_cube_data(importfolder3A + 'Propagation_metal_cluster-ELECTRON_DENSITY-1.cube') # import .cube file
cube_data_ref_10A, cube_atoms_ref_10A = read_cube_data(importfolder10A + 'Propagation_metal_cluster-ELECTRON_DENSITY-1.cube') # import .cube file


## import .cube file for reference electric field
cube_data_ref_3A_Ez, cube_atoms_ref_3A_Ez = read_cube_data(importfolder3A + 'Propagation_metal_cluster-efield_z-1.cube') # import .cube file
cube_data_ref_10A_Ez, cube_atoms_ref_10A_Ez = read_cube_data(importfolder10A + 'Propagation_metal_cluster-efield_z-1.cube') # import .cube file

# extract 2D cross-section for reference (background) e-density
cube_data_2d_ref_3A = cube_data_ref_3A[:, len(cube_data_ref_3A[:, 0, 0]) // 2, :]
cube_data_2d_ref_10A = cube_data_ref_10A[:, len(cube_data_ref_10A[:, 0, 0]) // 2, :]

# extract 2D cross-section for reference (background) electric field
cube_data_2d_ref_3A_Ez = cube_data_ref_3A_Ez[:, len(cube_data_ref_3A_Ez[:, 0, 0]) // 2, :]
cube_data_2d_ref_10A_Ez = cube_data_ref_10A_Ez[:, len(cube_data_ref_10A_Ez[:, 0, 0]) // 2, :]

# declare 2D density and Ez slice array
cube_data_2d_3A = np.zeros((len(cube_data_ref_3A[0,:,0]), len(cube_data_ref_3A[0,0,:]), timesteps))
cube_data_2d_3A_Ez = np.zeros((len(cube_data_ref_3A[0,:,0]), len(cube_data_ref_3A[0,0,:]), timesteps))

cube_data_2d_10A = np.zeros((len(cube_data_ref_10A[0,:,0]), len(cube_data_ref_10A[0,0,:]), timesteps))
cube_data_2d_10A_Ez = np.zeros((len(cube_data_ref_10A[0,:,0]), len(cube_data_ref_10A[0,0,:]), timesteps))

# import loop : import .cube data for each timestep, and slice it so we have just a cross section in x-y.
# full 3D data is thrown away to keep memory use down
for i in range(timesteps):#range(5):

    # Load the cube data
    # densities
    cube_data_3A, cube_atoms_3A = read_cube_data(importfolder3A + 'Propagation_metal_cluster-ELECTRON_DENSITY-1_' + str(i+1) + '.cube')
    cube_data_10A, cube_atoms_10A = read_cube_data(importfolder10A + 'Propagation_metal_cluster-ELECTRON_DENSITY-1_' + str(i+1) + '.cube')
    #fields
    cube_data_3A_Ez, cube_atoms_3A_Ez = read_cube_data(importfolder3A + 'Propagation_metal_cluster-efield_z-1_' + str(i+1) + '.cube')
    cube_data_10A_Ez, cube_atoms_10A_Ez = read_cube_data(importfolder10A + 'Propagation_metal_cluster-efield_z-1_' + str(i+1) + '.cube')

    print('imported t='+str(i))


    # Take a cross-section along the x-axis to get a cross-section of the y-z plane.
    #for densities
    cube_data_2d_3A[:,:,i] = cube_data_3A[:,len(cube_data_3A[:,0,0])//2,:] - cube_data_2d_ref_3A
    cube_data_2d_10A[:,:,i] = cube_data_10A[:,len(cube_data_10A[:,0,0])//2,:] - cube_data_2d_ref_10A
    #for fields
    cube_data_2d_3A_Ez[:,:,i] = cube_data_3A_Ez[:,len(cube_data_3A_Ez[:,0,0])//2,:] - cube_data_2d_ref_3A_Ez
    cube_data_2d_10A_Ez[:,:,i] = cube_data_10A_Ez[:,len(cube_data_10A_Ez[:,0,0])//2,:] - cube_data_2d_ref_10A_Ez

# calculate best density scale max and min for video

# for e density
if np.max(np.abs(cube_data_2d_3A)) > np.max(np.abs(cube_data_2d_10A)):
    min_lim = -np.max(np.abs(cube_data_2d_3A))/3.0
    max_lim = np.max(np.abs(cube_data_2d_3A))/3.0
else:
    min_lim = -np.max(np.abs(cube_data_2d_10A))/3.0
    max_lim = np.max(np.abs(cube_data_2d_10A))/3.0

# for Ez
if np.max(np.abs(cube_data_2d_3A_Ez)) > np.max(np.abs(cube_data_2d_10A_Ez)):
    min_lim_Ez = -np.max(np.abs(cube_data_2d_3A_Ez))
    max_lim_Ez = np.max(np.abs(cube_data_2d_3A_Ez))
else:
    min_lim_Ez = -np.max(np.abs(cube_data_2d_10A_Ez))
    max_lim_Ez = np.max(np.abs(cube_data_2d_10A_Ez))

# Get atomic positions and scale them to match the grid
# Transform atomic coordinates to match the grid used for electron density

# THIS COMPRESSES 3D POSITIONS INTO 2D
atomic_positions_3A = cube_atoms_3A.positions
atomic_positions_3A[:, 1] = atomic_positions_3A[:, 1] / cube_atoms_3A.cell[1, 1] * cube_data_3A.shape[1]
atomic_positions_3A[:, 2] = atomic_positions_3A[:, 2] / cube_atoms_3A.cell[2, 2] * cube_data_3A.shape[2]

y_atoms_3A = atomic_positions_3A[:, 1]
z_atoms_3A = atomic_positions_3A[:, 2]

atomic_positions_10A = cube_atoms_10A.positions
atomic_positions_10A[:, 1] = atomic_positions_10A[:, 1] / cube_atoms_10A.cell[1, 1] * cube_data_10A.shape[1]
atomic_positions_10A[:, 2] = atomic_positions_10A[:, 2] / cube_atoms_10A.cell[2, 2] * cube_data_10A.shape[2]

y_atoms_10A = atomic_positions_10A[:, 1]
z_atoms_10A = atomic_positions_10A[:, 2]


# Define the grid for e-density and field maps
y_3A = np.linspace(0, cube_data_3A.shape[1]-1, cube_data_3A.shape[1])
z_3A = np.linspace(0, cube_data_3A.shape[2]-1, cube_data_3A.shape[2])
Y_3A, Z_3A = np.meshgrid(y_3A, z_3A)

y_10A = np.linspace(0, cube_data_10A.shape[1]-1, cube_data_10A.shape[1])
z_10A = np.linspace(0, cube_data_10A.shape[2]-1, cube_data_10A.shape[2])
Y_10A, Z_10A = np.meshgrid(y_10A, z_10A)

# import transients
dipoles_time = np.genfromtxt('Dipoles.dat', skip_header = 1, delimiter=',')
fields_time = np.genfromtxt('Efield.dat', skip_header = 1, delimiter=',')

# plotting loop for creating each frame of the video
for i in range(timesteps):

    fig = plt.figure()
    gs = fig.add_gridspec(3, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    ax5 = fig.add_subplot(gs[2, :])

    fig.set_size_inches(10, 15)

    # Plot the 2D cross-section of electron density
    im_axs0 = ax1.pcolor(Y_3A.T, Z_3A.T, cube_data_2d_3A[:,:,i], cmap='PiYG', vmin=min_lim, vmax=max_lim)
    # Scatter the atoms on top
    ax1.scatter(y_atoms_3A, z_atoms_3A,s=3,  color='yellow', alpha = 1)
    # Add colorbar
    plt.colorbar(im_axs0, ax=ax1, fraction=0.056, pad=0.01)

    # remove ticks (I have no proper length scale anyhow....)
    ax1.tick_params(left = False, right = False , labelleft = False ,
                    labelbottom = False, bottom = False)

    ax1.set_aspect('equal')
    ax1.set_title('Electron density, 3 Angstrom')

    # Plot the 2D cross-section of electron density
    im_axs1 = ax2.pcolor(Y_10A.T, Z_10A.T, cube_data_2d_10A[:, :, i], cmap='PiYG', vmin=min_lim, vmax=max_lim)
    # Scatter the atoms on top
    ax2.scatter(y_atoms_10A, z_atoms_10A, s=3, color='yellow', alpha=1)
    # Add colorbar
    plt.colorbar(im_axs1, ax=ax2, fraction=0.056, pad=0.01)

    # remove ticks (I have no proper length scale anyhow....)
    ax2.tick_params(left=False, right=False, labelleft=False,
                    labelbottom=False, bottom=False)

    ax2.set_aspect('equal')
    ax2.set_title('Electron density, 10 Angstrom')

    # Plot the 2D cross-section of Ez
    im_axs0 = ax3.pcolor(Y_3A.T, Z_3A.T, cube_data_2d_3A_Ez[:,:,i], cmap='seismic', vmin=min_lim_Ez, vmax=max_lim_Ez)
    # Scatter the atoms on top
    ax3.scatter(y_atoms_3A, z_atoms_3A,s=3,  color='yellow', alpha = 1)
    # Add colorbar
    plt.colorbar(im_axs0, ax=ax3, fraction=0.056, pad=0.01)

    # remove ticks (I have no proper length scale anyhow....)
    ax3.tick_params(left = False, right = False , labelleft = False ,
                    labelbottom = False, bottom = False)

    ax3.set_aspect('equal')
    ax3.set_title('Electric field (z), 3 Angstrom')

    # Plot the 2D cross-section of Ez
    im_axs1 = ax4.pcolor(Y_10A.T, Z_10A.T, cube_data_2d_10A_Ez[:, :, i], cmap='seismic', vmin=min_lim_Ez, vmax=max_lim_Ez)
    # Scatter the atoms on top
    ax4.scatter(y_atoms_10A, z_atoms_10A, s=3, color='yellow', alpha=1)
    # Add colorbar
    plt.colorbar(im_axs1, ax=ax4, fraction=0.056, pad=0.01)

    # remove ticks (I have no proper length scale anyhow....)
    ax4.tick_params(left=False, right=False, labelleft=False,
                    labelbottom=False, bottom=False)

    ax4.set_aspect('equal')
    ax4.set_title('Electric field (z), 10 Angstrom')


    # note I flipped the polarity of this data for easy reading.
    ax5.plot(fields_time[:,0], fields_time[:,1]/13.0, label = '$E$')
    ax5.plot(dipoles_time[:,0], dipoles_time[:,1], label = '$p_{10}$')
    ax5.plot(dipoles_time[:,0], dipoles_time[:,2], label = '$p_{3}$')
    ax5.plot(dipoles_time[:,0], dipoles_time[:,3], label = '$p_{3}$ - $p_{10}$')

    ax5.axvline(time_array[i])
    ax5.legend(loc = 'upper left', frameon = False)
    ax5.set_xlim([10, 80])
    ax5.set_xlabel('Time (fs)')
    ax5.set_ylabel('Dipole')

    plt.tight_layout()

    # save each plot to file
    if i<10:
      plt.savefig('01_plotting_20_timesteps/frame=00' + str(i) + '.png')
    elif i<100:
      plt.savefig('01_plotting_20_timesteps/frame=0' + str(i) + '.png')
    else:
      plt.savefig('01_plotting_20_timesteps/frame=' + str(i) + '.png')

    print('plotted t='+str(i))

    plt.close()



# to make the frames into the movie, quicktime works well on mac.
