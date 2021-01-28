
from matplotlib.pylab import*
from numpy import*
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def get_pop(filename):
    
    data = genfromtxt(filename)
    tax = unique(data[:,2])
    wax = unique(data[:,3])
    Nt = len(tax)
    Nw = len(wax)
    pop = data[:,4].reshape([Nt,Nw])
    
    return tax,wax,pop

def get_rates(pop):   

    def expo_fit(x,y,p0):
        
        from scipy.optimize import curve_fit   
        def fitfunc(x,a,b,c):
            return a*exp(-x*b)+c  

        if any(p0)==None:
            p0 = [y[0],0.005,y[-1]]        
        (popt, pcov) = curve_fit(fitfunc,x,y,p0=p0)    
        p0=popt
        return popt[1],sqrt(pcov[1,1]),p0
    rate = []
    err = []
    en =[]
    p0 = None
    for iwax,xwax in enumerate(wax):
        if xwax>0.0 and xwax<=0.5 and iwax%10==0:
            z = pop[:,iwax]
            z1 = amax(z)
            imin = argmin(abs(z-z1))
            imin += 275
            imax = imin + 40
            icen = argmin(abs(z-amax(z)))
            xtax = copy(tax[imin:imax])
            zz = copy(z[imin:imax])
            r,yerr,p0=expo_fit(xtax,zz,p0)
            rate.append(r)
            err.append(yerr)
            en.append(xwax)
   
    rate = array(rate)
    err = array(err)
    en = array(en)
    return en,rate,err  

tax,wax,pop = get_pop('f.dat')
en,rate,rateerr = get_rates(pop)

fig,ax = subplots(1,1,figsize=(8,7))
rscale = 100 #To scale up the decay rate unit as it is arbitrary
ax.errorbar(en,rate*rscale,yerr=rateerr,fmt='o-',color = 'b',markersize=8,markeredgewidth=0,linewidth=2,label=r'$T = 100\,\mathrm{K}$')
ax.legend(loc=0,fontsize=20)
ax.tick_params(axis='both',labelsize=20)
ax.set_ylabel(r'$\mathrm{Decay\:Rates[arb.u]}$',fontsize=28)
ax.set_xlabel(r'$E-E_\mathrm{F}\mathrm{[eV]}$',fontsize=28)
# ax1.set_xlim(0,0.5)
# ax1.set_ylim(0,0.09*rscale)

savefig('decay_rates.png',bbox_inches='tight')



