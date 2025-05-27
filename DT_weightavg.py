import numpy
import matplotlib.pyplot as plt
import math

a = numpy.loadtxt('cat_tonga_nz.txt')
t0 = 1326250800
nhr = 12
tvec = (a[:,0]-t0)/3600+3
t = a[:,0]
tstep = 30
tmax = t0+nhr*3600

niter = nhr*3600/tstep

TEC = a[:,1]
SIDist = a[:,3]
el = a[:,2]
ew = numpy.sin(el*math.pi/180.0)


outfile = 'tonga_weightavg_timeseries_nz.txt'
ffo = open(outfile,'w')


tind = numpy.linspace(t0,t0+1440*tstep,num=1441)
dind = numpy.linspace(0,6000,num=121)


ATEC = numpy.zeros((len(tind),len(dind)))

for i  in range(0,len(tind)):
	for j in range(0,len(dind)):
		a1 = numpy.where((abs(t-tind[i]) < 300) & (abs(SIDist-dind[j]) < 200))[0]
		if (len(a1) > 3):
			wt = numpy.exp(-numpy.power((t[a1]-tind[i]),2)/2/30/30)
			wd = numpy.exp(-numpy.power((SIDist[a1]-dind[j]),2)/2/50/50)
			w = wt*wd
			ATEC[i,j] = numpy.sum(numpy.multiply(TEC[a1],w))/numpy.sum(w)

			tout = "{0:.1f}".format(float(tind[i]))
			dout = "{0:.1f}".format(float(dind[j]))
			tecout = "{0:.4e}".format(float(ATEC[i,j]))
			ffo.write(tout+' '+dout+' '+tecout+'\n')
		else:
			tout = "{0:.1f}".format(float(tind[i]))
			dout = "{0:.1f}".format(float(dind[j]))
			tecout = "{0:.4e}".format(float(ATEC[i,j]))
			ffo.write(tout+' '+dout+' '+tecout+'\n')

ffo.close()
#print(numpy.median(numpy.diff(dind)))
#plt.plot(meanTEC)
#plt.savefig('testmedian.png')
#plt.close()
