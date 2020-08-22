import numpy as np
import pylab as plt
import sys

#
# Disable warning when running code
#
import warnings
warnings.filterwarnings('ignore')


#
# Crude Gaussian plume model
#
def newmodel(source,wind,y,z,release_height,usetau):

    # Using notes from U. Wash. Using Pasquill Stability Category D (neutral)
    x = (ii+1)*usewind
    sigma_z = 32.093*(x/1000.)**0.81066
    theta   = 0.017453293*(8.3330 - 0.72382*np.log(x/1000.))
    sigma_y = 465.11628*(x/1000.)*np.tan(theta)
    
    C        = source/(2*np.pi*wind*sigma_y*sigma_z) 
    Zdiffuse = np.exp(-(z-release_height)**2/(2.*sigma_z**2)) + \
               np.exp(-(z+release_height)**2/(2.*sigma_z**2))
    Ydiffuse = np.exp(-y**2/(2*sigma_y**2))
    UseTime = x/wind
    chem = np.exp(-UseTime/usetau)

    return C*Zdiffuse*Ydiffuse,chem,C*Zdiffuse*Ydiffuse*chem

def ppt_ugm3(MW):
    return (2.69e25/6.023e23)*1e-12*MW/1e-6


#
# HFC-41
#
#usetau   = 365*24*3600  # 2.8 years
#MW       = 34.03
#LoD      = ppt_ugm3(MW) # ug/m3
#compound = 'HFC-41'
#
#
# HFC-161
#
usetau   = 73*24*3600 # 73 days ~ 0.2 years
MW       = 48.06
LoD      = ppt_ugm3(MW) # ug/m3
compound = 'HFC-161'
#
##
## HFC-152
##
#usetau   = 0.5*365*24*3600 # 0.5 years
#MW       = 66.05
#LoD      = ppt_ugm3(MW) # ug/m3
#compound = 'HFC-152'
#
##
## Propyl fluoride
##
#usetau   = 2*24*3600 # 1-2 days
#MW       = 62.09
#LoD      = ppt_ugm3(MW) # ug/m3
#compound = 'C3H7F'
#
##
## Isopropyl fluoride
##
#usetau   = 2*24*3600 # 1-2 days
#MW       = 62.09
#LoD      = ppt_ugm3(MW) # ug/m3
#compound = 'C3H7F'



tsteps    = 200*24*3600 # second
usewind   = 5 # m/s
source    = 9e6#2.5e7 # ug/sec
y = 0; z = 0
release_height = 0 # m

zsteps = np.arange(51)*5
ind    = np.where(zsteps==100)

Call  = []; Cchem = []; Cmix = []


dtsteps = 1*3600 # in jumps of a day

timeuse = []

for ii in np.arange(0,tsteps,dtsteps): 

    tmp   = []
    timeuse.append(ii)
    
    for zz in np.arange(len(zsteps)):

        usez = zsteps[zz]
    
        Cmix_out,Cchem_out,Call_out = newmodel(source,usewind,y,usez,release_height,usetau)

        tmp.append(Call_out)
        
    Call.append(tmp)#; Cchem.append(Cchem_out); Cmix.append(Cmix_out)


    
X,Y = np.meshgrid(np.arange(0,tsteps,dtsteps)*usewind/1000,zsteps)
Call = np.transpose(Call)


print(np.min(np.log(Call)),np.max(np.log(Call)))

minlvl = -1e-10#np.min(np.log(Call))
maxlvl = 20#np.max(np.log(Call))
print(minlvl,maxlvl)
lvls = np.arange(minlvl, maxlvl, (maxlvl - minlvl)/50)

#im = plt.contourf(X,Y,np.log(Call),levels=lvls,cmap=plt.get_cmap('jet'),extend='both')
#cb = plt.colorbar(im,label='ln(pollutant concentration) [$\mu$g/m$^3$]')
#plt.xlabel('Distance from source [km]')
#plt.ylabel('Sampling height above surface [m]')
#
#plt.savefig('contourplot_gauss.png')


plt.figure(20)

xdistance = np.arange(0,tsteps,dtsteps)*usewind/1000.

plt.plot(xdistance,\
         #np.log(Call[0,:]),c='black',label='Total: Mixing + Chemistry; Sampling @ surface')
         Call[0,:],c='black',label='Mixing & Chemistry; surface sampling')
#plt.plot(np.arange(0,tsteps,dtsteps)*usewind/1000.,\
#         Call[np.squeeze(ind),:], c='red',label='Total: Mixing + Chemistry; Sampling @ 100 m')

# Assuming here that we are using existing GC-MS instrumentation at the AGAGE sites which concentrates
# a 2-liter  volume air sample. We can detect 1 ppt with good precision, so that equates to 0.00139 ug/m3
# for HFC-41.  

ind = np.squeeze(np.where(Call[0,:] <= LoD))

print(xdistance[ind[0]])

plt.plot([0,6000],[LoD,LoD],'r--', label='LoD '+compound)
plt.title('Release height:'+str(release_height)+' km; Emission = '+str(source/1e6)+' g/s')#; LoD distance: '+str(xdistance[ind[0]])+' km')

plt.plot([xdistance[ind[0]],xdistance[ind[0]]],[0,0.25],'r-.',label='Distance where X ~ LoD '+compound)

plt.ylim([0,0.25])
plt.xlim([0,2000])

#plt.plot(np.arange(tsteps)*usewind/1000.,Cmix,c='red',label='Mixing')
plt.xlabel('Distance from source [km]'); plt.ylabel('Pollutant concentration [$\mu$g/m$^3$]')
plt.legend()

plt.savefig('ts_gauss.png')

plt.show()
