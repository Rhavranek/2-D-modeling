"""
Created on Tue Mar 15 15:07:30 2016

@author: rachelhavranek
2nd authors: jim mize!
"""

import numpy as np
import landlab as ll
import matplotlib.pyplot as plt

   
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
dt = 0.2 * (dx**2 / kappa)
print 'kappa:', kappa
print 'dt:', dt

# Create grid
mg = ll.RasterModelGrid(num_rows,num_cols,dx)


# Create data fields
z = mg.add_empty('node', 'depth')
x=mg.add_empty ('node', 'distance')
Q=mg.add_empty('node','heat flux') #should this be at every node or at Hedge?
T=mg.add_empty('node','temperature')
dTdz = mg.add_zeros('link', 'temperature_gradient')

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
    #print 'dTdz:', dTdz[2000]
    #print 'Q', Q[2000]
    #print 'dQdz:', dQdz[500]
    #print 'T:', T[500]
    
    ## this language is from a landlab tutorial 
    
"""FINALIZE"""

plt.figure(1)
im = plt.imshow(T, cmap=plt.cm.jet, extent=[0,1000,0,1000], origin='lower')  # display a colored image
plt.colorbar(im)
plt.title('Temperature Profile of Cooling Magma Body')

plt.show()

print('Done!')
