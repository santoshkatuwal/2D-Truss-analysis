import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
 
#tutorial video source: https://www.youtube.com/watch?v=Y-ILnLMZYMw
#%% Truss iproperties input
E=float(input("Enter value of E: "))
A=float(input("Enter value of A: "))
#E=1e4 #Ksi
#A=0.111 #in^2

path=os.getcwd()
data=pd.read_csv(path+'\\2D.csv')
NODES=data.Nodes.dropna()
X=data.x.dropna()
Y=data.y.dropna()

N_i=data.node_i.dropna()
N_j=data.node_j.dropna()
loadX=data.load_x.dropna()
loadY=data.load_y.dropna()
dofX=data.dof_x.dropna()
dofY=data.dof_y.dropna()


l=len(NODES)
elem=len(N_i)
nodes=[]    #defining matrix of nodes
bars=[]     #defining matrix for elements
P=[]        #defining matrix for loads
DOFCON=[]   #defining matrix for DOF

for i in range(l):
    nodes.append([X.iloc[i],Y.iloc[i]]) #Reading node geometry from csv file
for i in range(elem):
    bars.append([N_i.iloc[i],N_j.iloc[i]])  #Reading element geometry from excel file
for i in range(l):
    P.append([loadX.iloc[i],loadY.iloc[i]]) #Reading load value from csv file
for i in range(l):
    DOFCON.append([dofX.iloc[i],dofY.iloc[i]]) #Reading degree of freedoms from csv file

nodes=np.array(nodes) #converting list to ndarray
bars=np.array(bars) #converting list values to ndarray
P=np.array(P)
DOFCON=np.array(DOFCON)


sup=len(dofX)+len(dofY)-np.count_nonzero(dofX)-np.count_nonzero(dofY) 
    #counting numbers of dof where movement is not possible

Ur=[0]*sup #making displacement matrix where movement is not possible

#%% Truss analysis
def TrussAnalysis():
    NN=len(nodes) #number of nodes
    NE=len(bars) #number of elements
    DOF=2
    NDOF= DOF*NN #Total numbers of DOF = 2* number of nodes for 2D truss
    
    #Truss analysis
    d=nodes[bars[:,1],:]-nodes[bars[:,0],:] #creating matrix of [(xj-xi) (yj-yi)] for each bars
    L=np.sqrt((d**2).sum(axis=1))
    angle=d.T/L
    a=np.concatenate((-angle.T,angle.T),axis=1)
    K=np.zeros([NDOF,NDOF])
    for k in range(NE):
        aux=2*bars[k,:]
        index=np.r_[aux[0]:aux[0]+2,aux[1]:aux[1]+2]
        ES=np.dot(a[k][np.newaxis].T*E*A,a[k][np.newaxis])/L[k] #local stiffness matrix
        K[np.ix_(index,index)]=K[np.ix_(index,index)]+ES
        
    freeDOF=DOFCON.flatten().nonzero()[0] #Eleminating DOF at support
    supportDOF=(DOFCON.flatten()==0).nonzero()[0]
    Kff=K[np.ix_(freeDOF,freeDOF)]
    Kfr=K[np.ix_(freeDOF,supportDOF)]
    Krf=Kfr.T
    Krr=K[np.ix_(supportDOF,supportDOF)]
    Pf=P.flatten()[freeDOF]
    Uf=np.linalg.solve(Kff,Pf)
    U=DOFCON.astype(float).flatten()
    U[freeDOF]=Uf
    U[supportDOF]=Ur
    U=U.reshape(NN,DOF)
    u=np.concatenate ((U[bars[:,0]],U[bars[:,1]]),axis=1)
    N=E*A/L[:]*(a[:]*u).sum(axis=1)
    R=(Krf[:]*Uf).sum(axis=1)+(Krr[:]*Ur).sum(axis=1)
    R=R.reshape(int(sup/2),DOF)
    return np.array(N), np.array(R), U

#%%Plot
def Plot(nodes,c,lt,lw,lg):
    for i in range(len(bars)):
        xi, xf = nodes[bars[i,0],0],nodes[bars[i,1],0]
        yi, yf = nodes[bars[i,0],1],nodes[bars[i,1],1]
        line, = plt.plot([xi,xf],[yi,yf],color=c,linestyle=lt,lw=lw)
    line.set_label(lg)
    plt.legend(prop={'size':8})


#%% Result
N, R, U = TrussAnalysis() 
print('Axial Forces (+ve: tension, -ve: Compression)')
print(N[np.newaxis].T)
print('Reaction Forces (+ve: Upward, -ve: Downward)')
print(R)
print('Deformation at nodes')
print(U)

Plot(nodes,'gray','--',1,'Undeformed')
scale=1
Dnodes=U*scale+nodes
Plot(Dnodes,'k','-',3,'Deformed')

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    





