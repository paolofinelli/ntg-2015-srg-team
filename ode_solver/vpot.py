import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import subprocess

# set origin
# ----------

origin = 'lower'

# Define functions for reading data:
# =================================

#(1) real original
# ----------------

def get_demo_image():

#   from matplotlib.cbook import get_sample_data 
#   f = get_sample_data("axes_grid/bivariate_normal.npy", asfileobj=False)
    f=np.loadtxt("vkpot1.out")
   
    # z is a numpy array of 15x15
    return f, (0.,9.,0.,9.)


#(2) imaginary original
# ---------------------

def get_demo_image1():

#   from matplotlib.cbook import get_sample_data
   
#   f = get_sample_data("axes_grid/bivariate_normal.npy", asfileobj=False)
    f=np.loadtxt("vkpot100.out")
    
    # z is a numpy array of 15x15
    return f, (0.,9.,0.,9.)


#(3) real separable
# -----------------
  
def get_demo_image2():

#   from matplotlib.cbook import get_sample_data    
#   f = get_sample_data("axes_grid/bivariate_normal.npy", asfileobj=False)
    f=np.loadtxt("vkpot400.out")
    # z is a numpy array of 15x15
    return f, (0.,9.,0.,9.)

#(4) imaginary separable
# ----------------------

def get_demo_image3():

#   from matplotlib.cbook import get_sample_data
    
#   f = get_sample_data("axes_grid/bivariate_normal.npy", asfileobj=False)
    f=np.loadtxt("vkpot1000.out")
    # z is a numpy array of 15x15
    return f, (0.,9.,0.,9.)
 
#(5) momenta
# ----------
 
def get_demo_image4():

#   from matplotlib.cbook import get_sample_data
    
#   f = get_sample_data("axes_grid/bivariate_normal.npy", asfileobj=False)
    f=np.loadtxt("k1.out")
    # z is a numpy array of 15x15
    return f,(0., 9.) 
  

#(6) on-shel momenta
#-------------------
    
def get_demo_image5():

#   from matplotlib.cbook import get_sample_data

#   f = get_sample_data("axes_grid/bivariate_normal.npy", asfileobj=False)
    f=np.loadtxt("k0.out")
    # z is a numpy array of 15x15
    return f,(0., 9.)  

#------------------------------------------------------------
#------------------------------------------------------------

# Making the figure
# ================= 

extends = ["max", "min", "min", "min"]
headers= [ "(a) $V_{s=0}(k',k)$", "(b) $V_{s=0.001}(k'k)$",
"(c) $V_{s=0.004}(k',k)$", "(d) $V_{s=0.01}(k',k)$"]
cmap = plt.cm.get_cmap("gnuplot")
#cmap.set_under("magenta")
#cmap.set_over("yellow")
anx=['$s=0$','$s=0.001$','$s=0.004$','$s=0.01$'] 
# contour levels
# --------------

lv1 = [-.06, -.04, -.02,0.,.02,0.04,.06,.08]    # fig. (1,1)
lv2=  [ -0.04,-.03, -0.02,-.01,0.,.01,0.02,.03,0.04]   # fig. (1,2)
lv3= lv2#[-.12, -.08, -0.04, 0.,0.04,.08]   # fig. (2,1)
lv4= [-0.04, -.02,0.,.02,0.04,.06,0.08]   # fig. (2,2)

cbl1= [-1.2, -1.0,-0.8, -.6,-.4, -.2,1e-210,0.2]
#cbl1=[-1.2, -.8, -0.4, 0.,0.2,0.4]
cbl2=cbl1#cbl2= [-8., -6., -4., -2.,-1.,0.2,0.3]#
cbl3=cbl1 #[-.6, -.4, -0.2, 0.,0.2,0.4]
cbl4=cbl1#[-.3, -.2, -0.1, 0.,0.1,0.2,0.]

#ls=['solid','dashed','dashdot','dotted']
ls=['solid','dashed','dashdot','dashed','dotted']
lws=[1.5]

levs= [lv1,lv2,lv3,lv4]
cbls= [cbl1,cbl2,cbl3,cbl4]

#cbls=levs
levs=cbls

cls= ['orange','black','red','green','blue','black','green']


# x & ylabels
# -----------
 
xs= ["","","$k$","$k$"]
ys= ["$k$'","","$k$'",""]

# ticks
# -----

#xt=[0.,0.,1.,1.5,2.,2.5,3.,3.5,4.]

yt=[0.,2.5,5.,7.5,10.,12.5,15,17.5,20.]
#yt= [0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.]
xt=yt

# tick labels
# -----------

xt1= ['','','','']
xt2= ['','','','']
xt3= ['','',5.,'',10.,'',15.,'',20.]
xt4= [0.,'',5.,'',10.,'',15.,'',20.]  

yt1= [0.,'',5.,'',10.,'',15.,'',20.] 
yt2= ['','','','']
yt3= [0.,'',5.,'',10.,'',15.,'',20.]
yt4= ['','','','']

xtls= [xt1,xt2,xt3,xt4]
ytls= [yt1,yt2,yt3,yt4]

# axes limits
# -----------
dlim= 20.0

# start (2,2) plot
# ----------------

f, axs = plt.subplots(2,2)

# read data
# ---------

Z0, extent = get_demo_image()     # orig. real
Z1, extent1 = get_demo_image2()   # sep.  real
Z2, extent2 = get_demo_image1()   # orig. imag.
Z3, extent3 = get_demo_image3()   # sep.  imag. 

extents= [extent,extent2,extent1,extent3]

ZS=[Z0,Z2,Z1,Z3]

vmax, vmin = np.max(ZS), np.min(ZS)
import matplotlib.colors
norm = matplotlib.colors.Normalize(vmax=vmax, vmin=vmin)



# read k & k' values
# ------------------

X, xte= get_demo_image4()
Y, xte= get_demo_image4()

# read k_0 (on shell mom.)
# ------------------------
K0, xte= get_demo_image5()


#set colorplot limits
#-------------------- 

vmns= [-1.4, -1.4,-1.4,-1.4]
vmxs=[0.005,0.005,0.005,0.005]

#loop over figure panels
#-----------------------

for an,ax,xtl,ytl,xl,yl,vmin,vmax,levels,cbl,header, Z in zip(anx,axs.ravel(),xtls,ytls,xs,ys,vmns,vmxs,levs,cbls, headers, ZS):    

#   contour plot
#   ------------

    cs = ax.contour(X, Y, Z,colors=cls,origin=origin,levels=levels,linestyles=ls,interpolation='nearest',linewidths=lws)

#  line width
#  ----------

    lws=2.

#   panel title
#   -----------

    ax.set_title(header,size=30.5,position=(.5,1.05),color='black')

#   axis limits
#   ------------

    ax.set_xlim([0,dlim])
    ax.set_ylim([0,dlim])
      
#   color map
#   ---------

    im=ax.pcolor(X,Y,Z,cmap=cmap,vmin=vmin,vmax=vmax)

#    plt.clabel(cs, fmt = '%3.2f', colors = 'w', fontsize=14)

#   colorboxes
#   ----------


    cb=f.colorbar(cs, ax=ax,shrink=1.5,aspect=5,format='%4.2f')

    cb2=f.colorbar(im, ax=ax,ticks=cbl, shrink=1.,format='%4.2f')
    cb.add_lines(cs)

    
    cb.ax.tick_params(labelsize=30.5)
    cb2.ax.tick_params(labelsize=30.5)

    l,b,w,h = plt.gca().get_position().bounds
    ll,bb,ww,hh = cb2.ax.get_position().bounds


    cb.ax.set_position([ll+0.08,bb, ww,hh])
    cb2.ax.set_position([ll+0.005,bb, ww,hh])

    plt.xlabel(xl,fontsize=30.5)
    plt.ylabel(yl,fontsize=30.5)
    plt.xticks(xt)
    plt.yticks(yt)
    ax.set_xticklabels(xtl,fontsize=30.5)
    ax.set_yticklabels(ytl,fontsize=30.5)
  

    ax.xaxis.set_tick_params(width=3)
    ax.yaxis.set_tick_params(width=3)


    for axis in ['top','bottom','left','right']:
       ax.spines[axis].set_linewidth(3)

#    p1=plt.plot(X, K0,color='black',linewidth=2.,linestyle=ls[3],label='half-shell points')
#    p2=plt.plot(K0, Y,color='black',linewidth=2.,linestyle=ls[3])    


    cb.ax.get_children()[4].set_linewidths(3.0)
    cb.ax.get_children()[4].set_linestyles(ls)

    cb2.set_clim(vmin,vmax)

    ax.set_aspect(1,adjustable='box')
#    ax.annotate(an,xy=(15.5,15.5), xytext=(5, 8),color='blue',
#    fontsize=30.5)

#    ax.annotate('', xy=(3.5,3.5), xytext=(3.5, 5),color='white'
#    fontsize=20.5)
#    ax.set_position([0.1,0.1, .2, .2])  

legend = ax.legend(loc='upper left', shadow=True,fontsize=30.5)
#f.savefig('tsep_ecm5_os.pdf', format='PDF')
#f.savefig('tsep_ecm5_os.png', format='PNG')
#f.savefig('tsep_ecm5_os.eps', format='EPS')
plt.show()
