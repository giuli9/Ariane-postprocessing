#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:57:30 2023

@author: giuli9
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib.colors import BoundaryNorm
import seawater # using EOS80 for plotting isopycnals

ds1=xr.open_dataset("Lions/ariane_positions_quantitative.nc")
ds2=xr.open_dataset("Tyrrhenian/ariane_positions_quantitative.nc")
ds3=xr.open_dataset("Sicily/ariane_positions_quantitative.nc")

# Name of the initial and final section
# (forward in time)
initial_section1="Gulf of Lions"
initial_section2="Northern Tyrrhenian"
initial_section3="Sicily Strait"
final_section="Gibraltar Strait"

#=====================================================================
# the normalization factor is the cumulative transport of all and only
# particles which have reached final sections (no meandering, lost and 
# evaporated)
#=====================================================================
tot_transp=2883310065.5042796

# definition of variables

# Lions
init_transp1=ds1.init_transp.values
init_salt1=ds1.init_salt.values
init_temp1=ds1.init_temp.values
final_transp1=ds1.final_transp.values
final_salt1=ds1.final_salt.values
final_temp1=ds1.final_temp.values

# Tyrrhenian
init_transp2=ds2.init_transp.values
init_salt2=ds2.init_salt.values
init_temp2=ds2.init_temp.values
final_transp2=ds2.final_transp.values
final_salt2=ds2.final_salt.values
final_temp2=ds2.final_temp.values

# Sicily Strait
init_transp3=ds3.init_transp.values
init_salt3=ds3.init_salt.values
init_temp3=ds3.init_temp.values
final_transp3=ds3.final_transp.values
final_salt3=ds3.final_salt.values
final_temp3=ds3.final_temp.values

init_transp_sum1=np.sum(init_transp1)
final_transp_sum1=np.sum(final_transp1)

init_transp_sum2=np.sum(init_transp2)
final_transp_sum2=np.sum(final_transp2)

init_transp_sum3=np.sum(init_transp3)
final_transp_sum3=np.sum(final_transp3)

bins=100

xmin=35.5
xmax=39.5
ymin=11
ymax=29

delta_t=0.2
delta_s=0.04

def init_final(init_transp, init_salt, init_temp, final_transp, final_salt, final_temp):
	init=np.zeros((bins,bins))
	final=np.zeros((bins,bins))
	parcels=np.size(init_transp)

	for i in range(0,parcels):
		  ind_t=int(((init_temp[i]-ymin )/ delta_t).round())
		  ind_s=int(((init_salt[i]-xmin) / delta_s).round())
    
		  init[ind_t,ind_s]=init[ind_t,ind_s]+(((init_transp[i]/tot_transp))*100.)
		  ind_t=int(((final_temp[i]-ymin) / delta_t).round())
		  ind_s=int(((final_salt[i]-xmin) / delta_s).round())
		  final[ind_t,ind_s]=final[ind_t,ind_s]+(((final_transp[i]/tot_transp))*100.)
    
	return init, final   

lions=init_final(init_transp1, init_salt1, init_temp1, final_transp1, final_salt1, final_temp1)
tyrr=init_final(init_transp2, init_salt2, init_temp2, final_transp2, final_salt2, final_temp2)
sicily=init_final(init_transp3, init_salt3, init_temp3, final_transp3, final_salt3, final_temp3)

x_2D=np.zeros((bins,bins))
y_2D=np.zeros((bins,bins))

for i in range(0,bins):
    for j in range(0,bins):
        x_2D[j,i]= xmin + i*delta_s
        y_2D[j,i]= ymin + j*delta_t
        
x,y=np.meshgrid(x_2D[0,:],y_2D[:,0])

# masking values equal or less than 0.0
init_mask1 = np.ma.masked_less_equal(lions[0],0.0)
final_mask1 = np.ma.masked_less_equal(lions[1],0.0)

init_mask2 = np.ma.masked_less_equal(tyrr[0],0.0)
final_mask2 = np.ma.masked_less_equal(tyrr[1],0.0)

init_mask3 = np.ma.masked_less_equal(sicily[0],0.0)
final_mask3 = np.ma.masked_less_equal(sicily[1],0.0)

# weighted averages of temperature and salinity 
t1i=np.average(init_temp1, weights=init_transp1)
t1f=np.average(final_temp1, weights=final_transp1)

t2i=np.average(init_temp2, weights=init_transp2)
t2f=np.average(final_temp2, weights=final_transp2)

t3i=np.average(init_temp3, weights=init_transp3)
t3f=np.average(final_temp3, weights=final_transp3)

s1i=np.average(init_salt1, weights=init_transp1)
s1f=np.average(final_salt1, weights=final_transp1)

s2i=np.average(init_salt2, weights=init_transp2)
s2f=np.average(final_salt2, weights=final_transp2)

s3i=np.average(init_salt3, weights=init_transp3)
s3f=np.average(final_salt3, weights=final_transp3)


# creating the figure
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(18,29))

# defining my colormap
viridis_big = cm.get_cmap('plasma')
all_cmp = ListedColormap(viridis_big(np.linspace(0.97,0.06, 128)))

# calculating density anomalies sigma_0
levels=np.array([28.2,28.4,28.6,28.8,29,29.2,29.4])
rho_0 = seawater.dens0(x, y)
sigma_0= rho_0-1000.0

#===================================================================
# TS diagrams at Gibraltar (final section) and at each entry section
#===================================================================
def custome(axs,section):
	axs.set_xlim(37.5,39)
	axs.set_ylim(12.5,16)
	if section == "Gibraltar Strait":
		axs.set_title(section+" exit",size=27,fontweight='bold')
	else:
		axs.set_title(section+" entry",size=27,fontweight='bold')
	axs.set_xlabel('Salinity (psu)',size=24)
	axs.set_ylabel(r'$\theta$ (°C)',size=24)
	axs.tick_params(axis='x', labelsize=20,pad=10)
	axs.tick_params(axis='y', labelsize=20)
	positions = [37.5, 37.75, 38., 38.25, 38.5,38.75,39.]
	labels = ['37.5', '37.75', '38', '38.25', '38.5', '38.75','39']
	axs.set_xticks(positions,labels)
	
	return axs

#===================================================================
# final section entering at the Gulf of Lions
#===================================================================
bounds = [ 0, 0.01, 0.1, 1. ,  5, 10 , 15, 20 ,25, 28.5]
norm = BoundaryNorm(bounds, all_cmp.N)
mycmp=cm.ScalarMappable(norm=norm, cmap=all_cmp)

mesh1=axs[0,1].pcolormesh(x,y,init_mask1,cmap=all_cmp,shading='nearest',norm=norm) #vmin=1,vmax=28.5
mesh =axs[0,1].pcolor(x,y,lions[0],cmap="gray",alpha=0.15,facecolor = "none",edgecolor = 'k')
cs =axs[0,1].contour(x,y, sigma_0,levels=levels, colors='gray', zorder=1)
cl=axs[0,1].clabel(cs,fontsize=20,colors='black',inline=True,inline_spacing=-15)    
cbr=plt.colorbar(mesh1,ax=axs[0,1],ticks=[0,round(0.01,2), round(0.1,1), 1,round(6.5,1),round(12),round(17.5,1),round(23),round(28.5,1)],location='bottom')
cbr.ax.set_xticklabels(['{:.0f}'.format(0),'{:.2f}'.format(0.01),'{:.1f}'.format(0.1),'{:.0f}'.format(1),'{:.0f}'.format(5),'{:.0f}'.format(10),'{:.0f}'.format(15),'{:.0f}'.format(23),'{:.1f}'.format(28.5)])

cbr.ax.tick_params(labelsize=20)
# customize axes
axs[0,1]=custome(axs[0,1],final_section)

#==================================================================
# Gulf of Lions entry section
#==================================================================
bounds = [ 0,0.01, 0.1,1., 4.5, 8. ,11.5,15.,18.5,22.]
norm = BoundaryNorm(bounds, all_cmp.N)
mycmp=cm.ScalarMappable(norm=norm, cmap=all_cmp)

mesh2=axs[0,0].pcolormesh(x,y,final_mask1,cmap=all_cmp,shading='nearest',norm=norm) #vmin=1,vmax=22
mesh =axs[0,0].pcolor(x,y,lions[1],cmap='gray',alpha=0.15,facecolor = "none",edgecolor = 'k')
cbr=plt.colorbar(mesh2,ax=axs[0,0],ticks=[0,round(0.01,2),round(0.1,1),1,round(4.5,1),round(8),round(11.5,1),round(15),round(18.5,1),round(22)],location='bottom')
cbr.ax.tick_params(labelsize=20)
cbr.ax.set_xticklabels(['{:.0f}'.format(0),'{:.2f}'.format(0.01),'{:.1f}'.format(0.1),'{:.0f}'.format(1),'{:.1f}'.format(4.5),'{:.0f}'.format(8),'{:.1f}'.format(11.5),'{:.0f}'.format(15),'{:.1f}'.format(18.5),'{:.0f}'.format(22)])
cs = axs[0,0].contour(x,y, sigma_0,levels=levels, colors='gray', zorder=1)
cl=axs[0,0].clabel(cs,fontsize=20,colors='black',inline=True,inline_spacing=-15)

# customizing axes
axs[0,0]=custome(axs[0,0],initial_section1)
plt.figtext(0.5,0.646, "Transport percentage contribution",ha='center',va='center',size=27)

#==================================================================
# Northern Tyrrhenian entry section
#==================================================================
bounds = [0,0.01,0.02,0.03,0.04,0.05]
norm = BoundaryNorm(bounds, all_cmp.N)
mycmp=cm.ScalarMappable(norm=norm, cmap=all_cmp)

mesh2=axs[1,0].pcolormesh(x,y,final_mask2,cmap=all_cmp,shading='nearest',norm=norm) #vmin=0.01,vmax=0.05
mesh =axs[1,0].pcolor(x,y,tyrr[1],cmap='gray',alpha=0.15,facecolor = "none",edgecolor = 'k')
cbr=plt.colorbar(mesh2,ax=axs[1,0],ticks=[0,0.01,0.02,0.03,0.04,0.05],location='bottom')
cbr.ax.tick_params(labelsize=20)
cs = axs[1,0].contour(x,y, sigma_0,levels=levels, colors='gray', zorder=1)
cl=axs[1,0].clabel(cs,fontsize=20,colors='black',inline=True,inline_spacing=-15)
cbr.ax.set_xticklabels(['{:.0f}'.format(0),'{:.2f}'.format(0.01),'{:.2f}'.format(0.02),'{:.2f}'.format(0.03),'{:.2f}'.format(0.04),'{:.2f}'.format(0.05)])

# customizing axes
axs[1,0]=custome(axs[1,0],initial_section2)

#===================================================================
# final section entering at the Northern Tyrrhenian section
#===================================================================
bounds = [0,0.01, 0.05, 0.09, 0.13, 0.17, 0.21, 0.25]
norm = BoundaryNorm(bounds, all_cmp.N)
mycmp=cm.ScalarMappable(norm=norm, cmap=all_cmp)

mesh1=axs[1,1].pcolormesh(x,y,init_mask2,cmap=all_cmp,shading='nearest',norm=norm) #vmin=0.01,vmax=0.25
mesh =axs[1,1].pcolor(x,y,tyrr[0],cmap="gray",alpha=0.15,facecolor = "none",edgecolor = 'k')
cs =axs[1,1].contour(x,y, sigma_0,levels=levels, colors='gray', zorder=1)
cl=axs[1,1].clabel(cs,fontsize=20,colors='black',inline=True,inline_spacing=-15)  
cbr=plt.colorbar(mesh1,ax=axs[1,1],ticks=[round(0,0),round(0.01,2),round(0.05,2),round(0.09,2),0.13,round(0.17,2),0.21,0.25],location='bottom')
cbr.ax.tick_params(labelsize=20)
cbr.ax.set_xticklabels(['{:.0f}'.format(0),'{:.2f}'.format(0.01),'{:.2f}'.format(0.05),'{:.2f}'.format(0.09),'{:.2f}'.format(0.13),'{:.2f}'.format(0.17),'{:.2f}'.format(0.21),'{:.2f}'.format(0.25)])

# customizing axes
axs[1,1]=custome(axs[1,1],final_section)
plt.figtext(0.5,0.385, "Transport percentage contribution",ha='center',va='center',size=27)

#==================================================================
# Sicily Strait entry section
#==================================================================
bounds = [0,0.01,0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5 ]
norm = BoundaryNorm(bounds, all_cmp.N)
mycmp=cm.ScalarMappable(norm=norm, cmap=all_cmp)

mesh2=axs[2,0].pcolormesh(x,y,final_mask3,cmap=all_cmp,shading='nearest',norm=norm) #vmin=0.1,vmax=1.5
mesh =axs[2,0].pcolor(x,y,sicily[1],cmap='gray',alpha=0.15,facecolor = "none",edgecolor = 'k')
cbr=plt.colorbar(mesh2,ax=axs[2,0],ticks=[round(0,0),round(0.01,2),round(0.1,1),round(0.3,1),round(0.5,1),round(0.7,1),round(0.9,1),round(1.1,1),round(1.3,1),round(1.5,1)],location='bottom')
cbr.ax.tick_params(labelsize=20)
cs = axs[2,0].contour(x,y, sigma_0,levels=levels, colors='gray', zorder=1)
manual = [(0,0),(37.76,14),(0,0),(0,0),(38.6, 14), (38.65, 13.2), (38.94, 13.3)]
cl=axs[2,0].clabel(cs,fontsize=20,colors='black',inline=True,inline_spacing=-15, manual=manual)     
cbr.ax.set_xticklabels(['{:.0f}'.format(0),'{:.2f}'.format(0.01),'{:.1f}'.format(0.1),'{:.1f}'.format(0.3),'{:.1f}'.format(0.5),'{:.1f}'.format(0.7),'{:.1f}'.format(0.9),'{:.1f}'.format(1.1),'{:.1f}'.format(1.3),'{:.1f}'.format(1.5)])

# customizing axes
axs[2,0]=custome(axs[2,0],initial_section3)

#===================================================================
# final section entering at the Sicily Strait section
#===================================================================
bounds = [0,0.01, 0.1, 1, 1.5,2, 2.5, 3, 3.5, 4, 4.3 ]
norm = BoundaryNorm(bounds, all_cmp.N)
mycmp=cm.ScalarMappable(norm=norm, cmap=all_cmp)

mesh1=axs[2,1].pcolormesh(x,y,init_mask3,cmap=all_cmp,shading='nearest',norm=norm) # vmin=0.1,vmax=5
mesh =axs[2,1].pcolor(x,y,sicily[0],cmap="gray",alpha=0.15,facecolor = "none",edgecolor = 'k')
cs =axs[2,1].contour(x,y, sigma_0,levels=levels, colors='gray', zorder=1)
cl=axs[2,1].clabel(cs,fontsize=20,colors='black',inline=True, inline_spacing=-15)  
cbr=plt.colorbar(mesh1,ax=axs[2,1],ticks=[round(0,0),round(0.01,2),round(0.1,1),round(1,0),round(1.5,1),round(2,0),round(2.5,1),round(3,0),round(3.5,1),round(4,0),round(4.3,1)],location='bottom')
cbr.ax.tick_params(labelsize=20)
cbr.ax.set_xticklabels(['{:.0f}'.format(0),'{:.2f}'.format(0.01),'{:.1f}'.format(0.1),'{:.0f}'.format(1),'{:.1f}'.format(1.5),'{:.0f}'.format(2),'{:.1f}'.format(2.5),'{:.0f}'.format(3),'{:.1f}'.format(3.5),'{:.0f}'.format(4),'{:.1f}'.format(4.3)])

# customize axes
axs[2,1]=custome(axs[2,1],final_section)

# plotting the weighted mean
axs[0, 1].scatter(s1i, t1i, s=250, c='green', marker=(5, 1))
axs[0, 0].scatter(s1f, t1f, s=250, c='green', marker=(5, 1))
axs[1, 1].scatter(s2i, t2i, s=250, c='green', marker=(5, 1))
axs[1, 0].scatter(s2f, t2f, s=250, c='green', marker=(5, 1))
axs[2, 1].scatter(s3i, t3i, s=250, c='green', marker=(5, 1))
axs[2, 0].scatter(s3f, t3f, s=250, c='green', marker=(5, 1))

plt.subplots_adjust(hspace=0.1)
plt.figtext(0.5,0.125, "Transport percentage contribution",ha='center',va='center',size=27)

plt.savefig("TS.jpg",dpi=500, bbox_inches='tight')
