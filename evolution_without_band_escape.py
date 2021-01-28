#!/usr/bin/env python

from matplotlib.pylab import*
from numpy import*
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import matplotlib.gridspec as gridspec
norm = matplotlib.colors.Normalize(vmin=0, vmax=390, clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.plasma)

def fermi(x,T):    
    return 1.0/(exp(x/T) + 1.0)

def get_pop(filename):
    
    data = genfromtxt(filename)
    tax = unique(data[:,2])
    wax = unique(data[:,3])
    Nt = len(tax)
    Nw = len(wax)
    pop = data[:,4].reshape([Nt,Nw])    
    return tax,wax,pop

def get_electron_temp(pop):       

    def fermi_fit(x,y,p0):
        
        from scipy.optimize import curve_fit   
        def fitfunc(x,T,mu,a):
            return  a*fermi(x-mu,T)
        if all(p0)==None:
            p0 = [0.1,0.1,1.1]
        itmin = argmin(abs(x-0))
        itmax = argmin(abs(x-1.0))
        xx = copy(x[itmin:itmax])
        yy = y[itmin:itmax]
        (popt, pcov) = curve_fit(fitfunc,xx,yy,p0=p0)    
        p0=popt
        return popt[0],sqrt(pcov[0,0]),p0
    
    temp = []
    temperr = []
    time =[]
    mu = []
    p0=None
    for it,tt in enumerate(tax):
        if it>0 and it%2==0:
            z = pop[it,:]
            T,Terr,p0 = fermi_fit(wax,z,p0)
            temp.append(T)
            temperr.append(Terr)
            time.append(tt)
            mu.append(p0[1])            
   
    temp = array(temp)
    Terr = array(Terr)
    time = array(time)
    mu = array(mu)
    return time,temp,Terr,mu  

tax,wax,f = get_pop('f.dat')
time,temp,temperr,mu = get_electron_temp(f)

# Plot the transient electron distribution
fig,ax = subplots(1,1,figsize=(8,6))
ax.set_ylabel(r'$f(E-E_F)$',fontsize=24)

def do_plots(ax,pop,title):

    ax.tick_params(axis='both',labelsize=16)
    ax.set_ylim(0,1.1)

    ax.set_xlim(0,0.499)

    ax.set_xlabel(r'$E-E_F\mathrm{\,[eV]}$',fontsize=24)
    for it,tt in enumerate(tax):
        if it%10 ==0 and tt<499:
            ax.plot(wax,pop[it,:],'-',linewidth=2,color = mapper.to_rgba(tt))
do_plots(ax,f,'')
ax.tick_params(axis='both',labelsize=20)
ax.annotate(r'$\mathrm{Time}$',xy=(0.995,-0.08),xycoords='axes fraction',color='k',fontsize=22) 
ax.plot(wax,fermi(wax-mu[-1],0.0086),'k--',linewidth=1,label=r'$n_\mathrm{F}(T_\mathrm{lattice})$')
ax.legend(fontsize=28,frameon=False)
gs2 = gridspec.GridSpec(1, 1)
gs2.update(left=0.92, right=0.94, hspace=0.05)
ax3 = plt.subplot(gs2[0,0])
cb1 = matplotlib.colorbar.ColorbarBase(ax3, cmap=cm.inferno,norm=norm,orientation='vertical')  
cb1.set_ticks([0,100,200,300,390])
cb1.set_ticklabels(["0","2","4","6","8"])
cb1.ax.tick_params(labelsize=18) 

savefig('transient_distrbution.png',bbox_inches='tight')


# Plot the transient temperatures
fig,ax = subplots(1,1,figsize=(8,6))
conv = 11604 # 1eV = 11604 Kelvin
Tph = 100 # Lattice temperature
tscale = 50 # To scale down the time unit scale as it is arbitrary
ax.errorbar(time/tscale,temp*conv,yerr=temperr,fmt='r-',color = 'r',markeredgewidth=0,linewidth=3,label=r'$\mathrm{Electron}$')
ax.set_xlabel(r'$\mathrm{Time\:[arb.u]}$',fontsize=28)
ax.set_ylabel(r'$\mathrm{Temperature\,[K]}$',fontsize=28)
ax.tick_params(axis='both',labelsize=20)
ax.axhline(y=Tph,color='k',linestyle='--',linewidth=3,label=r'$\mathrm{Lattice}$')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1],loc=0,fontsize=22)
ax.set_xlim(0,400/tscale)
ax.set_ylim(0,temp[0]*conv+200)
savefig('transient_temp.png',bbox_inches='tight')



