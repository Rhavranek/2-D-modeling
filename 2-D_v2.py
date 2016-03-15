"""
Created on Tue Mar 15 15:07:30 2016

@author: rachelhavranek
2nd authors: jim mize!
"""

import numpy as np
import landlab as ll
import matplotlib as plt

   
"""initialize"""
#define variables 
k=2.5 #w/(m-k)
rho=2700. #kg/m3
Cp=2000. #J/kg K %rough approximation    
num_rows=201 #z array
num_cols=201 #x array
dx=40. #this will make both my dx and my dz equal in both directions
geotherm=25. #degrees C/km    
Ts=5. #surface temp degrees C
Tmagma=1500.
Qm=0.04 #w/m2
tmax=np.pi*(10**7)*200000. #total run time in seconds
kappa=k/(rho*Cp)
#dt=dx**2/((kappa*2)/4)
dt = 0.2 * (dx**2 / kappa)
print 'kappa:', kappa
print 'dt:', dt

# Create grid
mg = ll.RasterModelGrid(num_rows,num_cols,dx)
#core_nodes = mg.core_n

# Create data fields
z = mg.add_empty('node', 'depth')
x=mg.add_empty ('node', 'distance')
Q=mg.add_empty('node','heat flux') #should this be at every node or at Hedge?
T=mg.add_empty('node','temperature')
dTdz = mg.add_zeros('link', 'temperature_gradient')
#Tgrad=mg.calculate_gradients_at_active_links #can you do mg.add_empty here?
#dQdz=mg.calculate_gradients_at_active_links

#allgrades=mg.calculate_gradients_at_links
#Tgrad=[mg.active_links]

#Set Boundary Conditions
mg.set_closed_boundaries_at_grid_edges(True, True, True, False) #from Emily's code

#initialize arrays
z[:]=mg.node_y
x[:]=mg.node_x
T[:] = Ts + z*Qm/k

#this tells the computer that we have a square of magma at the bottom
for i in range (int(mg.number_of_nodes)):
    if 2000 < int(x[i]):
        if int(x[i]) < 6000:
            if int(z[i]) > 7000:
                T[i] = Tmagma



"""run"""
for j in range (int(tmax/dt)):
    dTdz[mg.active_links]=mg.calculate_gradients_at_active_links(T)
    Q=-k*dTdz
    dQdz= mg.calculate_flux_divergence_at_nodes(Q[mg.active_links])
    T[mg.core_nodes] += (-dQdz[mg.core_nodes])/(rho*Cp) * dt
   #update T
    print 'dTdz:', dTdz[2000]
    print 'Q', Q[2000]
    print 'dQdz:', dQdz[500]
    print 'T:', T[500]
    
    ## this language is from a landlab tutorial 
    # FINALIZE

        # Get a 2D array version of the elevations
        #zr = mg.node_vector_to_raster(z)

        # Create a shaded image
       #pylab.close()  # clear any pre-existing plot
        #im = pylab.imshow(zr, cmap=pylab.cm.RdBu, extent=[0,numcols*dx,0,numrows*dx],
                                         # origin='lower')
        # add contour lines with labels
        #cset = pylab.contour(zr, extent=[0,numcols*dx,numrows*dx,0], hold='on',
                                                # origin='image')
        #pylab.clabel(cset, inline=True, fmt='%1.1f', fontsize=10)

        # add a color bar on the side
       # cb = pylab.colorbar(im)
        #cb.set_label('Elevation in meters')

        # add a title and axis labels
        #pylab.title('Simulated topography with uplift and diffusion')
        #pylab.xlabel('Distance (m)')
        #pylab.ylabel('Distance (m)')

        # Display the plot
        #pylab.show()
        #print('Run time = '+str(time.time()-start_time)+' seconds')

#if __name__ == "__main__":
      
        
