import numpy as np
import pylab as plt
import sys

def readfile(filename,var1,var2,var3,var4,var5,var6,var7):

    f = open(filename,'r')
    header = f.readline()
    headercols = header.rsplit(';')
    #print(headercols)

    allcols = []

    while 1:
        line = f.readline()
        if not line: break
        cols = np.array(line.rsplit(';'))
        ind = np.where(cols == '')
        cols[ind] = -999
        cols = list(cols)
        allcols.append(cols)
    f.close()

    var1ind = np.squeeze(np.where(np.array(headercols) == var1))
    var2ind = np.squeeze(np.where(np.array(headercols) == var2))
    var3ind = np.squeeze(np.where(np.array(headercols) == var3))
    var4ind = np.squeeze(np.where(np.array(headercols) == var4))
    var5ind = np.squeeze(np.where(np.array(headercols) == var5))
    var6ind = np.squeeze(np.where(np.array(headercols) == var6))
    var7ind = np.squeeze(np.where(np.array(headercols) == var7))

    nlines = np.shape(allcols)[0]
    allcols = np.array(allcols)
    
    var1out = []
    var2out = []
    var3out = []
    var4out = []
    var5out = []
    var6out = []
    var7out = []
    sectors = []

    for ii in np.arange(nlines):
        var1out.append(np.float(allcols[ii,var1ind]))
        var2out.append(np.float(allcols[ii,var2ind]))
        var3out.append(allcols[ii,var3ind])
        var4out.append(np.float(allcols[ii,var4ind]))
        var5out.append(np.float(allcols[ii,var5ind]))
        var6out.append(allcols[ii,var6ind])
        var7out.append(np.float(allcols[ii,var7ind]))
        sectors.append(allcols[ii,4])
        
    return var1out, var2out, var3out, var4out, \
        var5out, var6out, var7out, sectors


def sanitycheck(sectors,activity,ef,emission):
    nlines = len(activity)
    print('*SANITY CHECK: ACTIVITY*EF = EMISSION*')
    for ii in np.arange(nlines):
        if ((ef[ii] != -999) and (emission[ii] != -999)):
            print(sectors[ii],activity[ii]*ef[ii]/1e3,emission[ii])
    else:
        print(sectors[ii],-999, -999)
    print('*--*')


def getrandomnumber(sigma,errmodel,nensemble):

    #
    # Sigma = 1-sigma uncertainty of variable being studies
    # Assume mu = 0 for these calculations
    #
    if errmodel == 'norm':
        values = np.random.normal(0,sigma,nensemble)
    elif errmodel == 'logn':
        values = np.random.lognormal(0,sigma,nensemble)
        
    return values
        

def runMCmodel(activity,activity_eps,\
               ef,ef_err,ef_errmodel,nensemble):

    # random numbers (normal or log-normal with a sigma determined by st of err)
    ef_eps       = getrandomnumber(ef_err,ef_errmodel,nensemble)

    emissionensemble = []

    for ii in np.arange(nensemble):
        emissionensemble.append((activity+activity_eps[ii])*(ef+ef_eps[ii]))

    return emissionensemble


    
        
    

filename = 'Uncertainty_result_DEU_2017.csv'

#
# CO2
#
activity, activity_err, activity_errmodel, \
    efco2, efco2_err, efco2_errmodel, \
    emissco2, sectors = readfile(filename,'AD','AD_stdev','AD_distr',\
                                 'EFCO2','EFCO2_stdev','EFCO2_distr',\
                                 'emission_CO2')
#
# CO
#
activity, activity_err, activity_err_model, \
    efco, efco_err, efco_errmodel, \
    emissco, sectors = readfile(filename,'AD','AD_stdev','AD_distr',\
                                'EFCO','EFCO_stdev','EFCO_distr',\
                                'emission_CO')
#
# NOX
#
activity, activity_err, activity_errmodel, \
    efnox, efnox_err, efnox_errmodel, \
    emissnox, sectors = readfile(filename,'AD','AD_stdev','AD_distr',\
                                'EFNOX','EFNOX_stdev','EFNOX_distr',\
                                 'emission_NOX')

#
# Sanity check: 
#
#sanitycheck(sectors,activity,efnox,emissnox)


#
# Monte Carlo calculation, approach based on Palmer et al 2006
#
nensemble = 10000

nsectors  = len(activity)

CO2emissionSectorMC = []
COemissionSectorMC  = []

# to study alt ensemble correlation approach
CO2emissionSectorMCALT = []
COemissionSectorMCALT  = []

import scipy.stats as stats

for ii in np.arange(nsectors):

    if efco[ii] != -999:

        activity_eps = getrandomnumber(activity_err[ii],activity_errmodel[ii],nensemble)

        CO2emissEnsemble = runMCmodel(activity[ii],activity_eps,
                                      efco2[ii],efco2_err[ii],efco2_errmodel[ii],nensemble)
        CO2emissionSectorMC += CO2emissEnsemble
        CO2emissionSectorMCALT.append(CO2emissEnsemble)        

        COemissEnsemble = runMCmodel(activity[ii],activity_eps,
                                     efco[ii],efco_err[ii],efco_errmodel[ii],nensemble)
        COemissionSectorMC += COemissEnsemble
        COemissionSectorMCALT.append(COemissEnsemble)

        # Correlations between individual sectors
        print(np.corrcoef(CO2emissEnsemble,COemissEnsemble)[0,1],sectors[ii])
        #print(stats.pearsonr(CO2emissEnsemble,COemissEnsemble),sectors[ii])
    
#plt.hist(CO2emissEnsemble,alpha=0.5)
#plt.hist(COemissEnsemble,alpha=0.5)


# Correlation from the sum of all sectors
sumCO2 = np.sum(CO2emissionSectorMCALT,axis=0)
sumCO  = np.sum(COemissionSectorMCALT,axis=0)
print('Correlation from sum of all sectors for all ensembles: ', np.corrcoef(sumCO2,sumCO)[0,1])
#plt.plot(sumCO2,sumCO,'ro')

# Correlation of all ensembles from all sectors
#print('Correlation from all ensembles from all sectors: ', np.corrcoef(CO2emissionSectorMC,COemissionSectorMC)[0,1])
#plt.subplot(212)
#plt.plot(CO2emissionSectorMC,COemissionSectorMC,'ro')

plt.show()
