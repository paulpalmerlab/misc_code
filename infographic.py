# Google LLC "Google COVID-19 Community Mobility Reports".
# https://www.google.com/covid19/mobility/ Accessed: <date>.

import pandas as pd
import numpy  as np
import pylab  as plt
import sys
from datetime import datetime
from time import mktime


import matplotlib
# Say, "the default sans-serif font is COMIC SANS"
matplotlib.rcParams['font.sans-serif'] = "Trebuchet MS"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"
# Define font size
matplotlib.rcParams.update({'font.size': 15})




fig, ax = plt.subplots(figsize=(15,5))




#
# Plot Google mobility data
#

UKdata = pd.read_excel('Global_Mobility_Report-220620.xls',\
                       sheet_name="GB_report")


UKdata.date = pd.to_datetime(UKdata.date,format='%d%m%Y')


#plt.plot(UKdata.date,UKdata.retail_and_recreation_percent_change_from_baseline,\
#         label='Retail and recreation')
#plt.plot(UKdata.date,UKdata.grocery_and_pharmacy_percent_change_from_baseline,\
#         label='Grocery and pharmacy')
#plt.plot(UKdata.date,UKdata.parks_percent_change_from_baseline,\
#         label='Parks')
#plt.plot(UKdata.date,UKdata.transit_stations_percent_change_from_baseline,\
#         label='Transit')



#
# Convert pandas dataframe to numpy
#
#UKdates   = UKdata.date.to_numpy()
#UKtransit = UKdata.transit_stations_percent_change_from_baseline.to_numpy()
#UKdates = np.array(UKdates)

startind = 15
ind = 108 #np.where(UKdates == '2020-06-01T00:00:00.000000000')

p1 = ax.bar(UKdata.date[startind:ind],UKdata.transit_stations_percent_change_from_baseline[startind:ind],\
       label='Transit',color='red',alpha=0.2)

#plt.plot(UKdata.date,UKdata.workplaces_percent_change_from_baseline,\
#         label='Workplaces')
#plt.plot(UKdata.date,UKdata.residential_percent_change_from_baseline,\
#         label='Residential')
plt.ylim([-100,100])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.xlabel('Date')
plt.ylabel('% change in transport sector from\n Jan 3–Feb 6 2020 baseline')

ax.yaxis.label.set_color('red')

plt.text(UKdata.date[60],-95, 'https://www.google.com/covid19/mobility/ Accessed: 22 June 2020',\
         fontsize=8,color='red')






#
# Plot daily UK deaths linked with Covid-19
# https://www.gov.uk/guidance/coronavirus-covid-19-information-for-the-public
#

CV19Death = pd.read_excel('2020-06-22_COVID-19_UK_deaths_time_series.xls',\
                       sheet_name="2020-06-22_COVID-19_UK_deaths_t")

death_dates = CV19Death['Publicly confirmed as deceased as of 5pm this day']
deaths_5pm  = CV19Death['UK Daily count of deaths in all settings (revised)']

deaths_5pm = np.divide(deaths_5pm,15)

ax2 = ax.twinx()

#
# Convert pandas dataframe to numpy
#
death_dates = death_dates.to_numpy()
deaths_5pm  = deaths_5pm.to_numpy()

#ind = np.squeeze(np.where(death_dates == '2020-06-01'))
ind = 89
print(death_dates[ind])

ax2.bar(death_dates[0:ind],deaths_5pm[0:ind],\
       color='green',alpha=0.25)

ax2.set_ylabel('UK Daily count of deaths in all settings')
ax2.yaxis.label.set_color('green')

ax2.spines['top'].set_visible(False)

plt.ylim([-100,100])

#ind = np.where(UKdata.date == '2020-06-01 00:00:00')
#UKdates = UKdata.date.to_numpy()
#plt.xlim([UKdates[0],UKdates[ind]])

plt.text(UKdata.date[60],95, \
         'https://coronavirus.data.gov.uk/about#total-and-daily-uk-cases Accessed: 22 June 2020',\
         fontsize=8,color='green')

#
# Cosmetic change - remove ticks above zero
#
yticks = ax.yaxis.get_major_ticks()
for ii in np.arange(5,9): yticks[ii].label1.set_visible(False)

yticks = ax2.yaxis.get_major_ticks()
for ii in np.arange(0,4): yticks[ii].label2.set_visible(False)

labels = ['','','','','0','375','750','1125','1500']
ax2.set_yticklabels(labels)








#
# CV19 timeline over the UK
#
cv19epochs = ['First reported UK death from Covid-19;\nmove from containment phase to delay phase',\
              'First reported death\nin Scotland',\
              'FCO advises against\nall non-essential international\ntravel',\
              'Start of Coronovirus job\nretention scheme',\
              'Start of UK lockdown',\
              'Wider range of businesses closed by law',\
              'Largest daily count of\nCOVID-19 deaths (1172)',\
              'First phase of lockdown\neasing in England',\
              'Lowest daily increase in deaths since\nintroducton of lockdown restrictions',\
              'First phase of\nlockdown easing\nin Scotland',\
              'Second phase of lockdown\neasing in England']#,\
              #'Social bubbles in England NI',\
              #'Opening of non-essential\nshops in England and NI']

cv19dates = ["2020-03-05",\
             "2020-03-13",\
             "2020-03-17",\
             "2020-03-20",\
             "2020-03-23",\
             "2020-03-26",\
             "2020-04-20",\
             "2020-05-13",\
             "2020-05-17",\
             "2020-05-29",\
             "2020-06-01"]#,\
             #"2020-06-13",\
             #"2020-06-15"]

halign = ['left',\
          'left',\
          'left',\
          'right',\
          'center',\
          'left',\
          'center',\
          'right',\
          'right',\
          'right',\
          'right']#,\
          #'center',\
          #'center']

cv19dates = pd.DataFrame(cv19dates)
cv19dates = cv19dates.astype('datetime64[ns]')

stemlengths = np.array([60,40,20,-30,-80,-10,20,20,40,20,-80])#,50,-30])

distancefromstem = np.array([5,5,5,-20,-12,-10,5,5,5,5,-20])#,5,-20])

plt.stem(cv19dates,stemlengths,'k',use_line_collection=True,basefmt=' ',markerfmt='ko')

for x, y, z, hh, zz in zip(np.array(cv19dates),np.array(stemlengths),np.array(cv19epochs), np.array(halign), distancefromstem):
    fontsize = 8
    plt.annotate(z, xy=(x,y), xytext=(0,zz), textcoords='offset points',ha=hh,fontsize=fontsize)

#print(death_dates[47]) 
#print(deaths_5pm[47])
    
#plt.plot([death_dates[46],death_dates[108]],[deaths_5pm[46],deaths_5pm[46]],'g--',alpha=0.2)



fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('infographic.png',dpi=500)
plt.show()
