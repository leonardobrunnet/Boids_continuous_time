#from pympler.tracker import SummaryTracker
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import math as math
import random as rand
import os
import sys

def print_name_location_output():
    ParStr='-V0_'+str(V0)+'-Dv_'+str(Dv)+'-N_'+str(N)+'-noise_'+str(noise)+'-Fadh_'+str(Fadh)
    path='data/SPP'+ParStr
    mkdir='mkdir -p '+path
    os.system(mkdir)
    fi=path+'/'+'parameter.in-'+ParStr
    fpar=open(fi,'w')
    ParDeltaPhi=path+'/'+'delta_phi.txt'
    fdelta_phi=open(ParDeltaPhi,'w')
    fdelta_phi.write("#time delta phi\n")
    fpar.write('path %s \n' % path)
    fpar.write('V0=%f \n' % V0)
    fpar.write('v1=%f \n' % v1)
    fpar.write('v2=%f \n' % v2)
    fpar.write('mu=%f \n' % mu)
    fpar.write('Frep=%f \n' % Frep)
    fpar.write('Fadh=%f \n' % Fadh)
    fpar.write('Req=%f \n' % Req)
    fpar.write('R0=%f \n' % R0)
    fpar.write('R=%f \n' % R)
    fpar.write('noise=%f \n' % noise)
    fpar.write('N=%f \n' % N)
    fpar.write('a=%f \n' % a)
    fpar.write('tf=%f \n' % tf)
    fpar.write('dt=%f \n' % dt)
    fpar.write('Lx=%d Ly=%d \n' % (L[0],L[1]))
    fpar.write('lbox=%d \n' % lbox)
    return path,ParStr,fdelta_phi


#tracker=SummaryTracker()
a=sys.argv[0:4]
#global variables and initialization
global N,ll,L,lbox,nb,dt,nb2,t,noise


#particle parameters
#V0=0.2      
#Dv=0.999*v0
#v1=v0-Dv/2
#v2=v0+Dv/2
v1=float(a[1])
v2=float(a[2])
Fadh=float(a[3])

Dv=v1-v2
V0=v1+Dv/2

mu=1
Frep=15.
#Fadh=.75     #    Fadh=0.75 #original model
Req=5./6.
R0=1.        #    interaction radius
noise=0.2

N=10
a=int(N/2)   # fraction of particles with speed v1
ll=60
L=np.array([ll,ll])
lbox=1       # box length depends on the interaction radius
nb=np.int64(L/lbox)
nb2=nb[1]*nb[0]
tf= 10010.0
dt=0.01
passos=int(tf/dt)
tau=1.0
exit_fig=10000
rand.seed()

c=[30,30] #center of the init circle
R=Req*math.sqrt(1.*N)/2 #radius of the init circle

#name and location of the output files 
path,ParStr,fdelta_phi=print_name_location_output()

#initialize N particles
V=np.zeros(N)
x=[]
y=[]
    
for k in range(N):
#option:  start with cells in a circle 
#    u = rand.random()
#    radius = R * math.sqrt(u)
#    theta = rand.uniform(0., 2 * math.pi)
#    x.append(ll/2 + radius * math.cos(theta))
#    y.append(ll/2 + radius * math.sin(theta))

    if k <= int(N/2):# 2 sub populations with different motility (red particles are slower)
        V[k]=v1
    else:
        V[k]=v2
 
#tracker.print_diff()


# Scatter image plots start here
def images(figindex,sizes,ParStr):

        plt.figure(figsize=(8,8))
        plt.axis([0,L[0]+Delta,0,L[1]+Delta])
        plt.axes().set_aspect(1.0)
        plt.title(ParStr)
        x,y,nx,ny=[],[],[],[]
        list(map(lambda i:x.append(i.r[0]), part))
        list(map(lambda i:y.append(i.r[1]), part))
        list(map(lambda i:nx.append(i.n[0]), part))
        list(map(lambda i:ny.append(i.n[1]), part))
        X=x[0:int(N/2)]
        XX=x[int(N/2):N]
        Y=y[0:int(N/2)]
        YY=y[int(N/2):N]
        plt.scatter(X,Y,s=sizes,c='r',alpha=0.3)  
        plt.scatter(XX,YY,s=sizes,c='b',alpha=0.3)
        name=path+'/'+str(figindex)+".png"
        fig = plt.gcf()
        plt.rc("savefig",dpi=200)
        fig.savefig(name,bbox_inches='tight')
        plt.close()
        figindex+=1
        return x,y,figindex

    
#Particle class definition
class particle:

    def __init__(self, x, y, vx, vy, ident,identclust,v0):
        
        self.v0=v0
        self.r = np.array([x,y])
        self.v= np.array([vx*self.v0,vy*self.v0])
        self.theta = np.arctan2(vy,vx)
        self.n = np.array([np.cos(self.theta),np.sin(self.theta)])
        self.ident = ident
        self.identclust= identclust
        self.Mybox = int(self.r[0]/lbox)+nb[0]*int(self.r[1]/lbox)
        self.Force = np.array([0.,0.])
        self.delta = 0.
        self.all_neigh=[]
        self.neigh_counter=0
        
    def mybox(self): #Each particle calculates the box it is in
        j=int(self.r[0]/lbox)+nb[0]*int(self.r[1]/lbox)
        return int(j)

    def changebox(self):
        newbox=self.mybox()
        if(newbox!=self.Mybox): #verify particle box change 
            box[self.Mybox].mylist.remove(self.ident) #this is safe since particles ident is unique
            box[newbox].mylist.append(self.ident)
            self.Mybox=newbox
        return self.Mybox

    def mov(self): #Particle movement
        
        self.v=self.v0*self.n+mu*self.Force
        dr=self.v*dt
        theta_noise=math.pi*(rand.random()-0.5)*2. # noise in the movment direction
        dr_noise=np.array([np.cos(theta_noise),np.sin(theta_noise)])*noise # noise in the movment itself
        self.r+=dr+dr_noise*np.sqrt(dt)
        self.autovelchange()
        self.theta+=self.dtheta*dt/tau 
        self.n[0]=np.cos(self.theta)
        self.n[1]=np.sin(self.theta)
        self.contour()

        return self.r,self.v,self.n

    
    def contour(self):
        self.r=(self.r+L)%L
        return self.r

    def autovelchange(self):
        vnorm=np.linalg.norm(self.v)
        u=self.v/vnorm
        self.dtheta=np.arcsin(self.n[0]*u[1]-self.n[1]*u[0])
        return self.dtheta

    def forces_between_particles(self):
        def force(self,dr):
            normdr=np.linalg.norm(dr)
            if(normdr<Req):
                f=Frep*(1/Req-1/normdr)*dr
            else:
                if(normdr<R0):
                    f=Fadh*dr*(1-Req/normdr)/(R0-Req)
                else:
                    f=0.0
            return f
            
        for i in box[self.Mybox].mylist:
            if(self.ident!=i):
                dr=part[i].r-self.r
                self.Force+=force(self,dr)
        for i in box[self.Mybox].neighboxlist:
                dr=part[i].r-self.r
                self.Force+=force(self,dr)
                part[i].Force-=force(self,dr)

    def zeroforce(self):
        self.Force=np.array([0.,0.])

    def zero_neigh_count(self):
        self.neigh_count=0
        
    def list_all_neigh(self):
        self.all_neigh=[]
        self.all_neigh.extend(box[self.Mybox].mylist)
        self.all_neigh.extend(box[self.Mybox].neighboxlist)
        return self.all_neigh

    def delta_neigh(self,x,y):
        if len(self.all_neigh) > 1 :
            for i in self.all_neigh:
                if i != self.ident:
                    dr=(self.r-part[i].r)
                    dr2=dif_in_box_square(dr)
                    dr_old=np.array([x[self.ident]-x[i],y[self.ident]-y[i]])
                    dr2_old=dif_in_box_square(dr_old)
                    this_delta=1-dr2_old/dr2
                    self.delta+=this_delta
                    part[i].delta+=this_delta
                    self.neigh_counter+=1
                    part[i].neigh_counter+=1
                    
    def my_delta(self):
        if self.neigh_counter > 0 :
            self.delta = self.delta/self.neigh_counter


    def randpos(self):
        x=ll*rand.random()
        y=ll*rand.random()
        self.r = np.array([x,y])
            
def dif_in_box_square(dr):
    dr[dr>ll/2]=float(ll)-dr[dr>ll/2]
    dr[dr<-ll/2]=float(ll)+dr[dr<-ll/2]
    dr2=np.linalg.norm(dr)**2
    return dr2
        

#Box class definition
class boite:
    def __init__(self,index):
        self.index = index
        self.mylist = []
        self.neighboxlist = []

    def neighbor_list(self):
        self.neighboxlist=[]
        zz1=(np.int64(self.index/nb[0])*nb[0]+(self.index+1)%nb[0])%nb2
        zz2=(np.int64((self.index+nb[0])/nb[0])*nb[0]+(self.index+nb[0]-1)%nb[0])%nb2
        zz3=(np.int64((self.index+nb[0])/nb[0])*nb[0]+(self.index+nb[0])%nb[0])%nb2
        zz4=(np.int64((self.index+nb[0])/nb[0])*nb[0]+(self.index+nb[0]+1)%nb[0])%nb2
        self.neighboxlist.extend(box[int(zz1)].mylist)
        self.neighboxlist.extend(box[int(zz2)].mylist)
        self.neighboxlist.extend(box[int(zz3)].mylist)
        self.neighboxlist.extend(box[int(zz4)].mylist)


    def neighbor_list_round(self): #this list includes also the left and upper boxes, it is used for he clustering algorithm
        self.neighboxlist_round=self.neighboxlist[:]
        zz1=(np.int64(self.index/nb[0])*nb[0]+(self.index-1)%nb[0])%nb2
        zz2=((np.int64(((self.index-nb[0]+nb2)%nb2)/nb[0])*nb[0]+(self.index+nb[0]-1)%nb[0])+nb2)%nb2
        zz3=((np.int64(((self.index-nb[0]+nb2)%nb2)/nb[0])*nb[0]+(self.index+nb[0])%nb[0])+nb2)%nb2
        zz4=((np.int64(((self.index-nb[0]+nb2)%nb2)/nb[0])*nb[0]+(self.index+nb[0]+1)%nb[0])+nb2)%nb2
        self.neighboxlist_round.extend(box[int(zz1)].mylist)
        self.neighboxlist_round.extend(box[int(zz2)].mylist)
        self.neighboxlist_round.extend(box[int(zz3)].mylist)
        self.neighboxlist_round.extend(box[int(zz4)].mylist)


#Main program                                  

# part=list(particle(res[i][0],res[i][1],rand.random()-0.5,rand.random()-0.5, i,i,V[i]) for i in range(N))

part=list(particle(L[0]*rand.random(),L[1]*rand.random(),rand.random()-0.5,rand.random()-0.5,i, i,V[i]) for i in range(N))

#part=list(particle(x[i],y[i],rand.random()-0.5,rand.random()-0.5, i,i,V[i]) for i in range(N))

#Image parameters*******

figindex=0
sizes=30.   # particles size on the image plots
Delta=0.

x,y,figindex=images(figindex,sizes,ParStr) #image of Init. Cond.
#************************

# boxes initialisation

nb2=np.int64(nb2)
box=list(boite(i) for i in range(nb2))
t=0

#tracker.print_diff()

# Construct list of particles in each box

for i in range(N):
    part[i].Mybox=part[i].mybox()
    box[np.int64(part[i].Mybox)].mylist.append(i)


#tracker.print_diff()

    # Construct list of particles in neighboring boxes

list(map(lambda i:i.neighbor_list(), box))
list(map(lambda i:i.list_all_neigh(), part)) #this is for the calculus of delta

#tracker.print_diff()


#System evolution

intt=0
traj1,traj11=list(),list()
traj2,traj22=list(),list()
traj3,traj33=list(),list()
traj4,traj44=list(),list()
traj5,traj55=list(),list()
traj6,traj66=list(),list()
traj7,traj77=list(),list()
traj8,traj88=list(),list()
#tracker.print_diff()


while(t<passos*dt):

    #Calculate the forces
    list(map(lambda i:i.forces_between_particles(), part))
    #Move all particles
    list(map(lambda i:i.mov(), part))
    t+=dt #update time
    #Find the newboxes 
    list(map(lambda i:i.changebox(), part))
    #Construct the list of particles in neighboring boxes
    list(map(lambda i:i.neighbor_list(), box))
    #Reset forces to zero
    list(map(lambda i:i.zeroforce(), part))

    #    tracker.print_diff()

    #Put particles in random positons to test delta
#    list(map(lambda i:i.randpos(), part))

    # Make a scatter graph
    intt+=1
    delta_average=0.0
    if intt%(exit_fig/10)== 0:
        print(t)
    if(intt%exit_fig==0):

        #Calculate delta ( Guillaume Gregoire, Hugues Chate, Yuhai
        #Tu. Moving and staying together without a leader. Physica D:
        #Nonlinear Phenomena, Elsevier, 2003, 181,
        #p. 157-170. <hal-00001037>

        if figindex > 0:
            list(map(lambda i:i.zero_neigh_count(), part))
            list(map(lambda i:i.delta_neigh(x,y), part))
            list(map(lambda i:i.my_delta(), part))
            delta_average=sum(map(lambda i:i.delta, part))/float(N)

        #Order parameter for collective movement

        phix=sum(map(lambda i:i.v[0], part))/float(N)
        phiy=sum(map(lambda i:i.v[1], part))/float(N)
        phi=math.sqrt(phix**2+phiy**2)

        #Printing order parameters
        
        fdelta_phi.write('%.2f %.4f %.4f\n' % (t,delta_average,phi))
        print('\n%.2f %.4f %.4f\n' % (t,delta_average,phi))
        
        #Instantaneous particle positions

        x,y,figindex=images(figindex,sizes,ParStr)

#        tracker.print_diff()


        #List of particles in boxes neighboring box i (right,down left, down center and down right) 

        list(map(lambda i:i.neighbor_list, box))

        #List of neighbors of particle i
        
        list(map(lambda i:i.list_all_neigh(), part))

        # record the positions of some individual cells over time

        # traj1.append(x[1]) 
        # traj11.append(y[1])             
        # traj2.append(x[a+1])
        # traj22.append(y[a+1])
        # traj3.append(x[3])
        # traj33.append(y[3])
        # traj4.append(x[a+3])
        # traj44.append(y[a+3])
        # traj5.append(x[10])
        # traj55.append(y[10])
        # traj6.append(x[a+10])
        # traj66.append(y[a+10])
        # traj7.append(x[49])
        # traj77.append(y[49])
        # traj8.append(x[a+49])
        # traj88.append(y[a+49])




# Plot the trajectories of some individual cells
plt.figure(figsize=(8,8))
plt.axis([0,L[0]+Delta,0,L[1]+Delta])
plt.axes().set_aspect(1.0)
plt.plot(traj1,traj11,'o')
plt.plot(traj2,traj22,'o')
plt.plot(traj3,traj33,'o')
plt.plot(traj4,traj44,'o')
plt.plot(traj5,traj55,'o')
plt.plot(traj6,traj66,'o')
plt.plot(traj7,traj77,'o')
plt.plot(traj8,traj88,'o')
name=path+"/trajectories"
fig = plt.gcf()
plt.rc("savefig",dpi=200)
fig.savefig(name,bbox_inches='tight')
fdelta_phi.close()

#********************Clustering***********************************


#2 particles are in the same cluster if the distance between them is < R0
for n in range(10): # 10 iterations to let the algorithm converge
    for p in part:
        for i in box[p.Mybox].mylist:
            if(p.identclust!=part[i].identclust):
                dr=part[i].r-p.r
                normdr=np.linalg.norm(dr)
                if normdr< R0:
                    part[i].identclust=p.identclust
        for i in box[p.Mybox].neighboxlist_round:
            if(p.identclust!=part[i].identclust):
                dr=part[i].r-p.r
                normdr=np.linalg.norm(dr)
                if normdr< R0:
                    part[i].identclust=p.identclust


# some clusters may not have converged (e.g. when 2 patches of cells are connected by only one cell)
# if the distance between 2 particles is < R0 but they have different cluster IDs, the cluster IDs of the 2 particles and
# of all the other particles sharing this ID are set equal to the cluster IDs of the first of the 2 particles
for p in part:
    for i in box[p.Mybox].mylist:
        dr=part[i].r-p.r
        normdr=np.linalg.norm(dr)
        if(p.identclust!=part[i].identclust) & (normdr< R0):
            for j in range(N):
                if part[j].identclust==part[i].identclust:
                    part[j].identclust=p.identclust
    for i in box[p.Mybox].neighboxlist_round:
        dr=part[i].r-p.r
        normdr=np.linalg.norm(dr)
        if(p.identclust!=part[i].identclust) & (normdr< R0):
            for j in range(N):
                if part[j].identclust==part[i].identclust:
                    part[j].identclust=p.identclust


clusterID =np.zeros((3,N))
for j in range(N):
    clusterID[0,j]=part[j].identclust # a list of the cluster IDs to which the particle belongs
    if part[j].v0==v1:# the subpopulation to which the particles belong (V1 or V2)
        clusterID[1,j]=1
    else:
        clusterID[2,j]=1

clusterSize=np.zeros((3,N))
for j in range(N):
    index= np.int64(clusterID[0,j])
    clusterSize[0,index]+=1 # the number of particles belonging to each cluster (=cluster size)
    clusterSize[1,index]+=clusterID[1,j] # the composition of each cluster (nb of V1 or V2 particles)
    clusterSize[2,index]+=clusterID[2,j]

compoclust = clusterSize[:,clusterSize[0,:]!=0]

np.savetxt('path/Clustcompo.csv',compoclust)
#l0: Clust size
#l1: Nb of " V1 particles" within the cluster



# ******************the order parameters*****************


#--------------------------------------------------------
#- for the gas-liquid interface

n=max(clusterSize[0,:])# n = size of the biggest cluster
print('n/N=',str(n/N))

#--------------------------------------------------------        
#- for the liqui-solid interface
print('delta=', str(delta_average))


#--------------------------------------------------------
#- for collective movement

print( 'phi=',str(phi))

#--------------------------------------------------------


#update the file for the phase diagram

path='/home/forget/data'
os.chdir(path)
import csv
row = [Fadh , V0, n/N , delta_average, phi ]

with open('final_output.csv', 'a') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerow(row)

csvFile.close()
#--------------------------------------------------------
