import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import math
import matplotlib.colors
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib import cm
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap


def psi_dens(path):
  
  # extrapolate variables
  
  ds=xr.open_dataset(path)
  xy_zonal=ds.xy_zonal.isel(nb_sect=0)
  xy_mer=ds.xy_mer.isel(nb_sect=0)
  xy_uh=ds.xy_uh.isel(nb_sect=0)
  xy_ruh=ds.xy_ruh.isel(nb_sect=0)
  xy_vh=ds.xy_vh.isel(nb_sect=0)
  xy_rvh=ds.xy_rvh.isel(nb_sect=0)

  imt_reg=ds.attrs["imt_reg"]
  jmt_reg=ds.attrs["jmt_reg"]
  imt_reg_start=ds.attrs["imt_reg_start"]-1
  jmt_reg_start=ds.attrs["jmt_reg_start"]-1
  imt_reg_end=ds.attrs["imt_reg_end"]
  jmt_reg_end=ds.attrs["jmt_reg_end"]
  kmt_reg=ds.attrs["kmt_reg"]
  kmt_reg_start=ds.attrs["kmt_reg_start"]-1
  kmt_reg_end=ds.attrs["kmt_reg_end"]
  
  bs=xr.open_dataset('mesh_mask_nemo.nc')
  tmask_reg=bs.tmask[0,0,jmt_reg_start:jmt_reg_end,imt_reg_start:imt_reg_end]
  
  # lat and lon of F points on C-grid
  xp=bs.glamf[0,:,:]
  yp=bs.gphif[0,:,:]
  
  # lat and lon of T points on C-grid
  xt=bs.glamt[0,:,:]
  yt=bs.gphit[0,:,:]
  
  tmask2=tmask_reg
  mtseg=np.zeros((np.shape(tmask_reg)))
  umask_reg=tmask_reg*0.
  vmask_reg=tmask_reg*0.
  
  # construction of umask and vmask
  for j in range(0,jmt_reg):
      for i in range(0,imt_reg-1):
          if (tmask_reg[j,i] + tmask_reg[j,i+1]) >= 1 :
              umask_reg[j,i]=1
                
  for j in range(0,jmt_reg-1):
      for i in range(0,imt_reg):
          if (tmask_reg[j,i] + tmask_reg[j+1,i]) >= 1 :
              vmask_reg[j,i]=1
  # define regional domain 
  with open(r"sections.txt", 'r') as fp:
      nb_sec=len(fp.readlines())
    
  segind=np.zeros(nb_sec)
  i1=np.zeros(nb_sec,int)
  i2=np.zeros(nb_sec,int)
  j1=np.zeros(nb_sec,int)
  j2=np.zeros(nb_sec,int)
  k1=np.zeros(nb_sec,int)
  k2=np.zeros(nb_sec,int)
  i1_reg=np.zeros(nb_sec,int)
  i2_reg=np.zeros(nb_sec,int)
  j1_reg=np.zeros(nb_sec,int)
  j2_reg=np.zeros(nb_sec,int)

  ip=0
  with open(r"sections.txt", 'r') as fp:
      for line in fp:
          data = line.split()
          segind[ip]=int(data[0])
          i1[ip]=abs(int(data[1]))-1
          i2[ip]=abs(int(data[2]))-1
          j1[ip]=abs(int(data[3]))-1
          j2[ip]=abs(int(data[4]))-1
          k1[ip]=abs(int(data[5]))-1
          k2[ip]=abs(int(data[6]))-1
        
          i1_reg[ip]=int(i1[ip]-imt_reg_start)
          i2_reg[ip]=int(i2[ip]-imt_reg_start)

          j1_reg[ip]=int(j1[ip]-jmt_reg_start)
          j2_reg[ip]=int(j2[ip]-jmt_reg_start) 
          ip+=1

  for ip in range(0,nb_sec):
      if i1_reg[ip] == i2_reg[ip]:
          for j in range(j1_reg[ip],j2_reg[ip]+1):
            if tmask_reg[j,i1_reg[ip]] == 1:
                tmask2[j,i1_reg[ip]] = 0
                mtseg[j,i1_reg[ip]] = 1
      


  for ip in range(0,nb_sec):
      if j1_reg[ip] == j2_reg[ip]:
          for i in range(i1_reg[ip],i2_reg[ip]+1):
              if tmask_reg[j1_reg[ip],i] == 1.:
                  tmask2[j1_reg[ip],i] = 0.
                  mtseg[j1_reg[ip],i] = 1.

  # pmask for Psi
  pmask_reg=np.zeros((jmt_reg,imt_reg))+1.
  
  for j in range(0,jmt_reg-1):
      for i in range(0,imt_reg -1):
          if  xy_zonal[j,i]==0.  and \
              xy_zonal[j,i+1]==0. and \
              xy_mer[j,i]==0.     and \
              xy_mer[j+1,i]==0. :
              pmask_reg[j,i]=math.nan

  # Compute divergence
  div=np.zeros((jmt_reg,imt_reg))
  for j in range(1,jmt_reg):
      for i in range(1,imt_reg):
          div[j,i]=xy_zonal[j,i]-xy_zonal[j,i-1]+ \
              xy_mer[j,i]-xy_mer[j-1,i]
        
  divmax=np.amax(np.absolute(div))

  ipb=np.zeros((jmt_reg,imt_reg),int)
  
  for j in range(0,jmt_reg-1):
      for i in range(0,imt_reg-1):
          if (tmask_reg[j,i]==tmask_reg[j+1,i+1]) and \
          (tmask_reg[j+1,i]==tmask_reg[j,i+1]):
              if (tmask_reg[j,i]+tmask_reg[j,i+1])==1. :
                  ipb[j,i]=1
                

  ip0=0
  jp0=0
  for j in range(0,jmt_reg-1):
      for i in range(0,imt_reg-1):
          if abs(div[j,i]) <= divmax and \
              abs(div[j,i+1]) <= divmax and \
              abs(div[j+1,i]) <= divmax and \
              abs(div[j+1,i+1]) <= divmax :
              if (xy_zonal[j,i] * xy_zonal[j+1,i])!=0. or \
              (xy_mer[j,i] * xy_mer[j,i+1])!=0.:
                  ip0=i
                  jp0=j
                  break
      if (ip0 != 0) and (jp0 != 0):
          break



  mp=np.zeros((jmt_reg,imt_reg),int)
  mp[jp0,ip0]=1

  mpold=0
  totmp=1
  while totmp>mpold:
      mpold=totmp
      for j in range(0,jmt_reg-1):
          for i in range(0,imt_reg-1):
              if ipb[j,i] == 0 and mp[j,i] ==1 :
                  if mp[j,i+1]==0 and (tmask2[j,i+1]== 1 \
                  or tmask2[j+1,i+1] == 1 ):
                      mp[j,i+1]=1
                  if mp[j+1,i] == 0 and (tmask2[j+1,i]==1 or \
                                      tmask2[j+1,i+1]==1):
                      mp[j+1,i]=1
                    
      for j in range(1,jmt_reg-1):
          for i in range(1,imt_reg-1):
              if ipb[j,i] == 0 and mp[j,i] ==1 :
                  if mp[j,i-1]==0 and (tmask2[j,i]== 1 \
                  or tmask2[j+1,i] == 1 ):
                      mp[j,i-1]=1
                  if mp[j-1,i] == 0 and (tmask2[j,i]==1 or \
                                      tmask2[j,i+1]==1):
                      mp[j-1,i]=1
                    
      for j in range(jmt_reg-2,-1,-1):
          for i in range(0,imt_reg-1):
              if ipb[j,i] == 0 and mp[j,i] ==1 :
                  if mp[j,i+1]==0 and (tmask2[j,i+1]== 1 \
                 or tmask2[j+1,i+1] == 1 ):
                      mp[j,i+1]=1
                  if mp[j+1,i] == 0 and (tmask2[j+1,i]==1 or \
                                      tmask2[j+1,i+1]==1):
                      mp[j+1,i]=1   
                    
      for j in range(jmt_reg-2,0,-1):
          for i in range(1,imt_reg-1):
              if ipb[j,i] == 0 and mp[j,i] ==1 :
                  if mp[j,i-1]==0 and (tmask2[j,i]== 1 \
                 or tmask2[j+1,i] == 1 ):
                      mp[j,i-1]=1
                  if mp[j-1,i] == 0 and (tmask2[j,i]==1 or \
                                      tmask2[j,i+1]==1):
                      mp[j-1,i]=1
                    
      for j in range(0,jmt_reg-1):
          for i in range(imt_reg-2,-1,-1):
              if ipb[j,i] == 0 and mp[j,i] ==1 :
                  if mp[j,i+1]==0 and (tmask2[j,i+1]== 1 \
                  or tmask2[j+1,i+1] == 1 ):
                      mp[j,i+1]=1
                  if mp[j+1,i] == 0 and (tmask2[j+1,i]==1 or \
                                      tmask2[j+1,i+1]==1):
                      mp[j+1,i]=1
                    
      for j in range(1,jmt_reg-1):
          for i in range(imt_reg-2,0,-1):
              if ipb[j,i] == 0 and mp[j,i] ==1 :
                  if mp[j,i-1]==0 and (tmask2[j,i]== 1 \
                  or tmask2[j+1,i] == 1 ):
                      mp[j,i-1]=1
                  if mp[j-1,i] == 0 and (tmask2[j,i]==1 or \
                                      tmask2[j,i+1]==1):
                      mp[j-1,i]=1
                    
      for j in range(jmt_reg-2,-1,-1):
          for i in range(imt_reg-2,-1,-1):
              if ipb[j,i] == 0 and mp[j,i] ==1 :
                  if mp[j,i+1]==0 and (tmask2[j,i+1]== 1 \
                  or tmask2[j+1,i+1] == 1 ):
                      mp[j,i+1]=1
                  if mp[j+1,i] == 0 and (tmask2[j+1,i]==1 or \
                                      tmask2[j+1,i+1]==1):
                      mp[j+1,i]=1    
                    
      for j in range(jmt_reg-2,0,-1):
          for i in range(imt_reg-2,0,-1):
              if ipb[j,i] == 0 and mp[j,i] ==1 :
                  if mp[j,i-1]==0 and (tmask2[j,i]== 1 \
                  or tmask2[j+1,i] == 1 ):
                      mp[j,i-1]=1
                  if mp[j-1,i] == 0 and (tmask2[j,i]==1 or \
                                      tmask2[j,i+1]==1):
                      mp[j-1,i]=1
                    
      totmp=sum(sum(mp))

  iref=np.zeros((jmt_reg,imt_reg),int)
  
  # Initialization on the Spain coast near the Gibraltar initial section
  iref_psi=5     
  jref_psi=21

  iref[jref_psi,iref_psi]=1
  irefold=0
  totiref=1

  while totiref != irefold:
      print(irefold,totiref)
      irefold=totiref
    
      for j in range(0,jmt_reg-1):
          for i in range(0,imt_reg-1):
              if iref[j,i] == 1 :
                  if iref[j,i+1] == 0 and tmask2[j+1,i+1]==0 and tmask2[j,i+1]==0:
                      iref[j,i+1]=1
    
                  if iref[j+1,i] == 0 and tmask2[j+1,i+1]==0 and tmask2[j+1,i]==0:
                      iref[j+1,i]=1
      for j in range(1,jmt_reg-1):
          for i in range(1,imt_reg-1):
              if iref[j,i] == 1:
                  if iref[j,i-1] == 0 and tmask2[j+1,i]==0 and tmask2[j,i]==0:
                      iref[j,i-1]=1

                  if iref[j-1,i] == 0 and tmask2[j,i+1]==0 and tmask2[j,i]==0:
                      iref[j-1,i]=1
                    
      for j in range(jmt_reg-2,-1,-1):
          for i in range(0,imt_reg-1):
              if iref[j,i] == 1:
                  if iref[j,i+1] == 0 and tmask2[j+1,i+1]==0 and tmask2[j,i+1]==0:
                      iref[j,i+1]=1

                  if iref[j+1,i] == 0 and tmask2[j+1,i+1]==0 and tmask2[j+1,i]==0:
                      iref[j+1,i]=1
                    
      for j in range(jmt_reg-2,0,-1):
          for i in range(1,imt_reg-1):
              if iref[j,i] == 1:
                  if iref[j,i-1] == 0 and tmask2[j+1,i]==0 and tmask2[j,i]==0:
                      iref[j,i-1]=1
    
                  if iref[j-1,i] == 0 and tmask2[j,i+1]==0 and tmask2[j,i]==0:
                      iref[j-1,i]=1
            
      for j in range(0,jmt_reg-1):
          for i in range(imt_reg-2,-1,-1):
              if iref[j,i] == 1:
                  if iref[j,i+1] == 0 and tmask2[j+1,i+1]==0 and tmask2[j,i+1]==0:
                      iref[j,i+1]=1

                  if iref[j+1,i] == 0 and tmask2[j+1,i+1]==0 and tmask2[j+1,i]==0:
                      iref[j+1,i]=1
                    
      for j in range(1,jmt_reg-1):
          for i in range(imt_reg-2,0,-1):
              if iref[j,i] == 1:
                  if iref[j,i-1] == 0 and tmask2[j+1,i]==0 and tmask2[j,i]==0:
                      iref[j,i-1]=1
    
                  if iref[j-1,i] == 0 and tmask2[j,i+1]==0 and tmask2[j,i]==0:
                      iref[j-1,i]=1
            
      for j in range(jmt_reg-2,-1,-1):
          for i in range(imt_reg-2,-1,-1):
              if iref[j,i] == 1:
                  if iref[j,i+1] == 0 and tmask2[j+1,i+1]==0 and tmask2[j,i+1]==0:
                      iref[j,i+1]=1
    
                  if iref[j+1,i] == 0 and tmask2[j+1,i+1]==0 and tmask2[j+1,i]==0:
                      iref[j+1,i]=1
            
      for j in range(jmt_reg-2,0,-1):
          for i in range(imt_reg-2,0,-1):
              if iref[j,i] == 1:
                  if iref[j,i-1] == 0 and tmask2[j+1,i]==0 and tmask2[j,i]==0:
                      iref[j,i-1]=1
    
                  if iref[j-1,i] == 0 and tmask2[j,i+1]==0 and tmask2[j,i]==0:
                      iref[j-1,i]=1
            
      totiref=sum(sum(iref))
                    
        
  jp1=0
  ip1=0

  for j in range(0,jmt_reg):
      for i in range(0,imt_reg):
          if iref[j,i]==1 and mp[j,i]==1 and ip1==0:
              ip1=i
              jp1=j
        
  ip1=2
  jp1=20


  psi=np.zeros((jmt_reg,imt_reg))-1e12
  psi[jp1,ip1]=0.
  ipsi=np.zeros((jmt_reg,imt_reg),int)
  ipsi[jp1,ip1]=1
  totipsi=1
  while totipsi<mpold:
      print(totipsi,mpold)
      for j in range(0,jmt_reg):
          for i in range(0,imt_reg-1):
              if ipsi[j,i] == 1 and ipb[j,i] == 0 :
                  if ipsi[j,i+1] == 0 and mp[j,i+1] == 1 and vmask_reg[j,i+1] == 1 :
                      psi[j,i+1]=psi[j,i]+xy_mer[j,i+1]
                      ipsi[j,i+1]=1
                  
      for j in range(0,jmt_reg):
          for i in range(1,imt_reg):
              if ipsi[j,i] == 1 and ipb[j,i] == 0:
                  if ipsi[j,i-1] == 0 and mp[j,i-1] == 1 and vmask_reg[j,i] == 1:
                      psi[j,i-1]=psi[j,i]-xy_mer[j,i]
                      ipsi[j,i-1]=1
                    
   
      for j in range(0,jmt_reg-1):
          for i in range(0,imt_reg):
              if ipsi[j,i] == 1 and ipb[j,i] == 0 :
                  if ipsi[j+1,i] == 0 and mp[j+1,i] == 1 and umask_reg[j+1,i] == 1 :
                      psi[j+1,i]=psi[j,i]-xy_zonal[j+1,i]
                      ipsi[j+1,i]=1
                  

      for j in range(1,jmt_reg):
          for i in range(0,imt_reg):
              if ipsi[j,i] == 1 and ipb[j,i] == 0:
                  if ipsi[j-1,i] == 0 and mp[j-1,i] == 1 and umask_reg[j,i] == 1:
                      psi[j-1,i]=psi[j,i]+xy_zonal[j,i]
                      ipsi[j-1,i]=1
                    
 
      for j in range(jmt_reg-1,-1,-1):
          for i in range(0,imt_reg-1):
              if ipsi[j,i] == 1 and ipb[j,i] == 0:
                  if ipsi[j,i+1] == 0 and mp[j,i+1] == 1 and vmask_reg[j,i+1] == 1:
                      psi[j,i+1]=psi[j,i]+xy_mer[j,i+1]
                      ipsi[j,i+1]=1
                   

      for j in range(jmt_reg-1,-1,-1):
          for i in range(1,imt_reg):
              if ipsi[j,i] == 1 and ipb[j,i] == 0:
                  if ipsi[j,i-1] == 0 and mp[j,i-1] == 1 and vmask_reg[j,i] == 1:
                      psi[j,i-1]=psi[j,i]-xy_mer[j,i]
                      ipsi[j,i-1]=1
                    

      for j in range(jmt_reg-2,-1,-1):
          for i in range(0,imt_reg):
              if ipsi[j,i] == 1 and ipb[j,i] == 0:
                  if ipsi[j+1,i] == 0 and mp[j+1,i] == 1 and umask_reg[j+1,i] == 1:
                      psi[j+1,i]=psi[j,i]-xy_zonal[j+1,i];
                      ipsi[j+1,i]=1
                    
  
      for j in range(jmt_reg-1,0,-1):
          for i in range(0,imt_reg):
              if ipsi[j,i] == 1 and ipb[j,i] == 0:
                  if ipsi[j-1,i] == 0 and mp[j-1,i] == 1 and umask_reg[j,i] == 1 :
                      psi[j-1,i]=psi[j,i]+xy_zonal[j,i]
                      ipsi[j-1,i]=1
                    
  
      for j in range(jmt_reg-1,-1,-1):
          for i in range(imt_reg-2,-1,-1):
              if ipsi[j,i] == 1 and ipb[j,i] == 0:
                  if ipsi[j,i+1] == 0 and mp[j,i+1] == 1 and vmask_reg[j,i+1] == 1:
                      psi[j,i+1]=psi[j,i]+xy_mer[j,i+1]
                      ipsi[j,i+1]=1
                   

      for j in range(jmt_reg-1,-1,-1):
          for i in range(imt_reg-1,0,-1):
              if ipsi[j,i] == 1 and ipb[j,i] == 0:
                  if ipsi[j,i-1] == 0 and mp[j,i-1] == 1 and vmask_reg[j,i] == 1:
                      psi[j,i-1]=psi[j,i]-xy_mer[j,i]
                      ipsi[j,i-1]=1
                    
 
      for j in range(jmt_reg-2,-1,-1):
          for i in range(imt_reg-1,-1,-1):
              if ipsi[j,i] == 1 and ipb[j,i] == 0:
                  if ipsi[j+1,i] == 0 and mp[j+1,i] == 1 and umask_reg[j+1,i] == 1:
                      psi[j+1,i]=psi[j,i]-xy_zonal[j+1,i]
                      ipsi[j+1,i]=1
                    
 
      for j in range(jmt_reg-1,0,-1):
          for i in range(imt_reg-1,-1,-1):
              if ipsi[j,i] == 1 and ipb[j,i] == 0:
                  if ipsi[j-1,i] == 0 and mp[j-1,i] == 1 and umask_reg[j,i] == 1:
                      psi[j-1,i]=psi[j,i]+xy_zonal[j,i]
                      ipsi[j-1,i]=1
                    
   
      for j in range(0,jmt_reg):
          for i in range(imt_reg-2,-1,-1):
              if ipsi[j,i] == 1 and ipb[j,i] == 0:
                  if ipsi[j,i+1] == 0 and mp[j,i+1] == 1 and vmask_reg[j,i+1] == 1:
                      psi[j,i+1]=psi[j,i]+xy_mer[j,i+1]
                      ipsi[j,i+1]=1
                    
   
      for j in range(0,jmt_reg):
          for i in range(imt_reg-1,0,-1):
              if ipsi[j,i] == 1 and ipb[j,i] == 0:
                  if ipsi[j,i-1] == 0 and mp[j,i-1] == 1 and vmask_reg[j,i] == 1:
                      psi[j,i-1]=psi[j,i]-xy_mer[j,i]
                      ipsi[j,i-1]=1
                    
 
      for j in range(0,jmt_reg-1):
          for i in range(imt_reg-1,-1,-1):
              if ipsi[j,i] == 1 and ipb[j,i] == 0:
                  if ipsi[j+1,i] == 0 and mp[j+1,i] == 1 and umask_reg[j+1,i] == 1:
                      psi[j+1,i]=psi[j,i]-xy_zonal[j+1,i]
                      ipsi[j+1,i]=1
                    
    
      for j in range(1,jmt_reg):
          for i in range(imt_reg-1,-1,-1):
              if ipsi[j,i] == 1 and ipb[j,i] == 0:
                  if ipsi[j-1,i] == 0 and mp[j-1,i] == 1 and umask_reg[j,i] == 1:
                      psi[j-1,i]=psi[j,i]+xy_zonal[j,i]
                      ipsi[j-1,i]=1
                    
    
    
                    
      totipsi=sum(sum(ipsi))

  psi[psi<-1e10]=np.nan
  Psi=psi/1e6
  Psi=Psi*pmask_reg
  
  # Forward in time
  Psi=-1*Psi
  
  pmax=np.nanmax(Psi)
  pmin=np.nanmin(Psi)

  xt_reg=xt[jmt_reg_start:jmt_reg_end,imt_reg_start:imt_reg_end]
  yt_reg=yt[jmt_reg_start:jmt_reg_end,imt_reg_start:imt_reg_end]

  xp_reg=xp[jmt_reg_start:jmt_reg_end,imt_reg_start:imt_reg_end]
  yp_reg=yp[jmt_reg_start:jmt_reg_end,imt_reg_start:imt_reg_end]

  dens_v=np.zeros(np.shape(xy_rvh))
  dens_u=np.zeros(np.shape(xy_ruh))
  dens_v[:,0:imt_reg]=xy_rvh[:,0:imt_reg]/xy_vh[:,0:imt_reg]
  dens_u[0:jmt_reg,:]=xy_ruh[0:jmt_reg,:]/xy_uh[0:jmt_reg,:]

  dens=(dens_u+dens_v)/2

  dens2=dens*pmask_reg
  
  return Psi, dens2

def axes(ax):
  ax.set_extent([-8, 17, 33, 44])
  ax.coastlines('10m')
  ax.add_feature(cartopy.feature.NaturalEarthFeature("physical","land",'10m', \
                                                edgecolor='none', facecolor="lightgrey"), zorder=0)
  bc=ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,zorder=0)
  bc.top_labels=False
  bc.right_labels=False
  bc.xformatter = LONGITUDE_FORMATTER
  bc.yformatter = LATITUDE_FORMATTER
  bc.xlabel_style = {'size': 24}
  bc.ylabel_style = {'size': 24}
  return ax


# start plotting

fig, (ax1, ax2, ax3) = plt.subplots(3,sharex=True,figsize=(18,22),subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw = {'hspace':0.01})


path='Tyrrhenian/ariane_statistics_quantitative.nc'


#plt.figure(figsize=(13,7))
#ax3=plt.axes(projection=ccrs.PlateCarree())

# Tyrrhenian
levels=np.arange(-0.31,0.11,0.06)

for i in range (0,levels.size,1):
    if levels[i]>0:
        bs=axe(ax2).contour(xp_reg,yp_reg,Psi_dens(path)[0],[levels[i]],colors='black',linestyles='dashed',transform=ccrs.PlateCarree(),linewidths=1.5)
        axe(ax2).clabel(bs, inline=True,fmt='% 2.2f', fontsize=13)
    if levels[i]<0:
        bs=axe(ax2).contour(xp_reg,yp_reg,Psi_dens(path)[0],[levels[i]],colors='black',linestyles='solid',transform=ccrs.PlateCarree(),linewidths=1.5)
        axe(ax2).clabel(bs, inline=True,fmt='% 2.2f', fontsize=13)
        



# plot the sections

for ip in range(0,nb_sec-1):
    axe(ax2).plot([xt_reg[j1_reg[ip],i1_reg[ip]],xt_reg[j2_reg[ip],i2_reg[ip]]],[yt_reg[j1_reg[ip],i1_reg[ip]]\
                                                                              ,yt_reg[j2_reg[ip],i2_reg[ip]]],\
         transform=ccrs.PlateCarree(),color='black', linestyle='-')


# Definying colormap
colors=['#c6e1f1','#c2c5f3','#7ee98e','#f69253','#dd4c94',"#0d9edf"]
cmap2=ListedColormap(colors)
ticks=[27.0,27.5,28.6,28.8,29.0,29.1,29.2]
norm = BoundaryNorm(ticks, cmap2.N)
x,y=np.meshgrid(xt_reg[0,:],yt_reg[:,0])
vmax=29.2
vmin=27.0
ff=axe(ax2).pcolormesh(x,y,Psi_dens(path)[1],shading='nearest',cmap=cmap2,norm=norm,vmin=round(vmin,1),vmax=round(vmax,1),transform=ccrs.PlateCarree(),zorder=1,alpha=0.8)


axe(ax2).text(-7.5, 41.1, 'Northern Tyrrhenian \nentry 1% = 0.01 Sv',fontweight='bold', fontsize = 27,alpha=1)


axes=np.array([ax1,ax2,ax3])
cbar=fig.colorbar(ff,ax=axes.ravel().tolist(),ticks=ticks,norm=norm,boundaries=ticks,spacing='proportional', shrink=0.9, location='right')
cbar.set_label(r'$\sigma_0$(kg/m$^3$)',size=28)
cbar.ax.set_yticklabels(ticks,fontsize=26) 

plt.savefig('Psi.jpg',dpi=450,bbox_inches='tight')


