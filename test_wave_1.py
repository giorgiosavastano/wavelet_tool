import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

periods_massimi = []
for i in xrange( len( lista ) ):
    inde = np.argmax( maxs_general[i] )
    periods_massimi.append( periods_general[i][inde] )

periods_massimi_array = np.asarray( periods_massimi )    
maschera = (periods_massimi_array > 400) & (periods_massimi_array < 1800)
percentuale1 = len( periods_massimi_array[maschera] ) / 262.0 * 100.0

plt.figure( 'tst2hist' )
plt.hist( periods_massimi )
plt.ylabel('number of values')
plt.xlabel('Periods [s]')
plt.title('Histogram of Wave Periods')
plt.grid()
plt.show()
  


period_4 = []
period_7 = []
period_8 = []
period_10 = []
period_20 = []

for i in xrange( len(lista) ):
    if '4' in lista[i]:
        period_4.append( periods_massimi[i] )
    elif '7' in lista[i]:
        period_7.append( periods_massimi[i] )    
    elif '8' in lista[i]:
        period_8.append( periods_massimi[i] )
    elif '10' in lista[i]:
        period_10.append( periods_massimi[i] )
    elif '20' in lista[i]:
        period_20.append( periods_massimi[i] )