import numpy as np
import matplotlib.pyplot as m
import math as math
import random as r

#Particle class definition
class particle:
    noise=0.6
    v0 = 1.0
    mu=1.0
    Frep=30.0
    Fadh=0.75
    Req=5./6.
    R0=1.0
    def __init__(self, x, y, vx, vy, ident):
        self.r = np.array([x,y])
        self.v =  np.array([vx*self.v0,vy*self.v0])
        self.theta = np.arctan2(vy,vx)
        self.n = np.array([np.cos(self.theta),np.sin(self.theta)])
        self.ident = ident
        self.Mybox = int(self.r[0]/lbox)+n*int(self.r[1]/lbox)
        self.Force =np.array([0.,0.])
        
    
    def mybox(self): #Each particle calculates the box it is in
        j=int(self.r[0]/lbox)+n*int(self.r[1]/lbox)
        return j

    def changebox(self):
        newbox=self.mybox()
        if(newbox!=self.Mybox): #verify particle box change 
            box[self.Mybox].mylist.remove(self.ident) #this is safe since particles ident is unique
            box[newbox].mylist.append(self.ident)
            self.Mybox=newbox
        return self.Mybox

    def mov(self): #Particle moviment
        self.v=self.v0*self.n+self.mu*self.Force
        self.r+=self.v*dt
        self.autovelchange(self.v,self.n)
        self.theta+=self.dtheta*dt/tau+self.noise*(r.random()-0.5)*np.sqrt(dt)
        self.contour(self.r)
        self.n[0]=np.cos(self.theta)
        self.n[1]=np.sin(self.theta)
        return self.r,self.v,self.n

    
    def contour(self,r):
        self.r=(self.r+L)%L
        return self.r

    def autovelchange(self,v,n):
        vnorm=np.linalg.norm(v)
        u=self.v/vnorm
        self.dtheta=np.arcsin(self.n[0]*u[1]-self.n[1]*u[0])
        return self.dtheta

    def forces_between_particles(self):
        def force(self,dr):
            if(np.linalg.norm(dr)<self.Req):
                f=self.Frep*(np.linalg.norm(dr)/self.Req-1)*dr/np.linalg.norm(dr)
            else:
                if(np.linalg.norm(dr)<self.R0):
                    f=self.Fadh*(dr-self.Req)/(self.R0-self.Req)
                else:
                    f=dr*0
            return f
        for i in box[self.Mybox].mylist:
            if(self.ident!=i):
                dr=part[i].r-self.r
                self.Force+=force(self,dr)
        for i in box[self.Mybox].neighboxlist:
                dr=part[i].r-self.r
                f=force(self,dr)
                self.Force+=force(self,dr)
                part[i].Force-=force(self,dr)

    def zeroforce(self):
        self.Force=np.array([0.,0.])
        

#Box class definition
class boite:
    def __init__(self,index):
        self.index = index
        self.mylist = []
        self.neighboxlist = []

    def neighbor_list(self):
        self.neighboxlist=[]
        i1=(self.index/n*n+(self.index%n+1)%n)%n2  #right box
        i2=((self.index/n+1)*n+(self.index%n-1+n)%n)%n2  #lower left box
        i3=((self.index/n+1)*n+(self.index%n))%n2  #lower center box
        i4=((self.index/n+1)*n+(self.index%n+1)%n)%n2  #lower right box
        for i in [i1, i2,i3,i4]:
            self.neighboxlist.extend(box[i].mylist)

#Main program                                  
#global variables
global N,L,lbox,n,dt,n2
N=500
L=60
lbox=1
n=L/lbox
n2=n*n                  
dt=0.05
exit_fig=30
tau=1.0
passos=25000

#initialize N particles

part=list(particle(L*r.random(),L*r.random(), r.random()-0.5,r.random()-0.5, i) for i in range(N))
box=list(boite(i) for i in range(n2))
t=0

# Construct list of particles in each box

for i in range(N):
    part[i].Mybox=part[i].mybox()
    box[part[i].Mybox].mylist.append(i)

# Construct list of particles in neighboring boxes

map(lambda i:i.neighbor_list(), box)

#System evolution
figindex=0
intt=0
while(t<passos*dt):

#Calculate the forces
    map(lambda i:i.forces_between_particles(), part)
#Move all particles
    map(lambda i:i.mov(), part)
    t+=dt #update time
#Find the newboxes 
    map(lambda i:i.changebox(), part)
#Construct the list of particles in neighboring boxes
    map(lambda i:i.neighbor_list(), box)
#Reset forces to zero
    map(lambda i:i.zeroforce(), part)

# Make a scatter graph
    intt+=1
    if(intt%exit_fig==0):
        print(t)
        m.axis([0,L,0,L])
        x,y,nx,ny=[],[],[],[]
        map(lambda i:x.append(i.r[0]), part)
        map(lambda i:y.append(i.r[1]), part)
        map(lambda i:nx.append(i.n[0]), part)
        map(lambda i:ny.append(i.n[1]), part)
        m.quiver(x,y,nx,ny)
        name=str(figindex)+".png"
        m.savefig(name, bbox_inches='tight')
        figindex+=1
        m.clf()
