

import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
import scipy.special as spl
import scipy.integrate as integrate
import pandas as pd
import numba_scipy as scp
from numba import jit,njit
import numba as nb

from NumbaQuadpack import quadpack_sig, dqags





#tBlist=np.array((8)) #enter the list of bath bandwidths that you want to simulate for
tB=8
pie=np.pi
# Time parameters
h= 0.1 # enter the timestep to be used for simulation
Time_max = 100
Net_time=np.int64(np.floor(Time_max/h)) #
testrange =Net_time


Lambda= 1 #system bath coupling
t_electron=1 #electron bandwidth= = σ
a1=1  #lattice constant
Temp_electron=0.1# electron temperature
Temp_bath= 0.8 #the bath temperature
mu_electron = 1  # chemical potential of the electron
mu_bath = 1 #chem potential of the bath



#volume parameters
sitenum = 19 #the no. of sites in the lattice, same as number of momentum modes to be simulated - enter an ODD number
a2=2*pie*(1/(sitenum*a1)) #reciprocal space lattice constant
V_ph = np.arange(-0.5*sitenum*a2,0.5*a2*(sitenum+1),a2)#collect(-0.5*sitenum*a2:a2:0.5*a2*(sitenum+1))


# electron parameters
@njit
def fermi(e,T,mu):
    return 1/(np.exp((e-mu)*(1/T))+1)

@njit
def energy_electron(k):
    return t_electron*(1-np.cos(V_ph[k]*a1))

@njit
def Gretarded(k,t1,t2):
    return (t1>=t2)*complex(0,-1)*np.exp(complex(0,-1)*energy_electron(k)*(t1-t2)*h)

@njit
def Gkeldysh(k,t1,t2,Telectron,mu):
    return -complex(0,1)*np.tanh((energy_electron(k)-mu)/(2*Telectron))*np.exp(-complex(0,1)*energy_electron(k)*(t1-t2)*h)

#Gkeldysh =lambda k,t1,t2,Telectron,mu: -complex(0,1)*np.tanh((energy_electron(k)-mu)/(2*Telectron))*np.exp(-complex(0,1)*energy_electron(k)*(t1-t2)*h)



#bath parameters
@njit
def J(w,tB):
    return (2/tB)*np.sqrt( 1- (w/(2*tB))**2 )
#J = lambda w,tB: (2/tB)*np.sqrt( 1- (w/(2*tB))**2 )




@njit
def SigmaR(t1,t2,tB):
    if t1>t2:
        sum1= -complex(0,1)*(Lambda**2)*(1)*(1/(tB))*( spl.j1(2*tB*abs(t1-t2)*h) / (abs(t1-t2)*h ) )
        return sum1
    else:
        return 0
@njit
def sigmak(w,tB):
    return -complex(0,1)*(Lambda**2)*J(w,tB)*np.tanh((w-mu_bath)/(2*Temp_bath))
#=lambda w,tB: -complex(0,1)*(Lambda**2)*J(w,tB)*np.tanh((w-mu_bath)/(2*Temp_bath))

@njit
def SigmaK(t1,t2,tB):
    dw=(1/10000)*4*tB
    steps = np.arange(-2*tB,2*tB,dw) #collect(-2*tB:dω:2*tB)
    result=0
    for w in steps:
        result = result + dw*sigmak(w,tB)*np.exp(-complex(0,1)*w*(t1-t2)*h)

    return result/(2*pie)



### spawning Gr & Gk matrices

@njit
def matinit():
    A=[]
    for i in range(len(V_ph)+2):
        A.append( np.zeros((Net_time+5,Net_time+5),dtype=np.complex128))#Array{ComplexF64,2}(undef,Net_time+5,Net_time+5)
        #Gkmatrix.append( np.zeros((Net_time+5,Net_time+5),dtype=np.complex128))
    return A


## Box Initialization ###
boxinitindex=1

@njit
def boxinit_gr(A):

    for k in range(len(V_ph)):
        for i in range(Net_time):
            A[k][i,i] = -complex(0,1) #exactly true           ## Gr(t,t)≂̸0

                ###actual Box Initialization ###

    #GF Initialization

    for k in range(len(V_ph)):
        for i in range(boxinitindex+1):
            for j in range(boxinitindex+1):
                A[k][i,j] = Gretarded(k,i,j)



@njit
def boxinit_gk(A):
    for k in range(len(V_ph)):
        for i in range(boxinitindex+1):
            for j in range(boxinitindex+1):
                A[k][i,j] = Gkeldysh(k,i,j,Temp_electron,mu_electron)


#matrix definitions
@njit
def gimme_a_matrix():
    return  np.zeros(shape=(Net_time+5,Net_time+5),dtype=np.complex128)



Grmatrix =np.array(matinit())      #  #numba doesn't know what Gr array contains, I use a tuple to get around this
Gkmatrix = np.array(matinit()) #Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

SigmaRmatrix = np.array(gimme_a_matrix())
SigmaKmatrix = np.array(gimme_a_matrix())

boxinit_gr(Grmatrix)
boxinit_gk(Gkmatrix)





Grmatrix[5].shape
Gkmatrix[5].shape
#### Convolution definitions
@njit
def Sgma_conv_Gr(k,t1,t2):
    sum=0
    if t1>t2:
        sum= SigmaRmatrix[t1,t2]*Grmatrix[k][t2,t2]*(h/2)
        for i in range(t2+1,t1):
            sum = sum + SigmaRmatrix[t1,i]*Grmatrix[k][i,t2]*h
        return sum

    else:
        return 0


Sgma_conv_Gr(2,10,1)

@njit
def Sgma_conv_GK(k,t1,t2):
    if t1>1:
        sum=0
        sum = sum + SigmaRmatrix[t1,1]*Gkmatrix[k][1,t2]*(h/2)

        for i in range(2,t1):               #2:t1-1
            sum = sum+SigmaRmatrix[t1,i]*Gkmatrix[k][i,t2]*h

        return sum

    else:
        return 0



@njit
def Sgma_conv_GA(k,t1,t2):
    if t2>1:
        sum=0
        sum=sum+SigmaKmatrix[t1,1]* np.conjugate(Grmatrix[k][t2,1]) * h*(1/2)  #starting 1/2
        sum=sum+ SigmaKmatrix[t1,t2]* np.conjugate(Grmatrix[k][t2,t2]) * h*(1/2) #ending 1/2

        for i in range(2,t2):# i=2:t2-1:
            sum=sum+ SigmaKmatrix[t1,i]* np.conjugate(Grmatrix[k][t2,i]) * h   #middle ones, they get h & not h/2 due to double addition

        return sum

    else:
        return 0


### Keeping stuff here that helps in plotting the results

energyrange=np.zeros(len(V_ph))
for k in range(len(V_ph)):
    energyrange[k]=energy_electron(k)

@njit
def quadtest(x,en,tB):
    return (1/pie)*fermi(x,Temp_bath,mu_bath)*(Lambda**2)*(J(x,tB))*( 1/( (x-en)**2 + (Lambda**2*J(x,tB))**2 ) )

Net_time
testrange


def newres(e,hop):  # computes n(E) i.e. equilibrium occupation at energy E using bath spectral function
    return integrate.quad(quadtest, -2*hop,2*hop,args=(e,hop))


#
# @njit
# def even_newer_res(e,hop):
#     dw=(1/10000)*4*tB
#     steps = np.arange(-2*tB,2*tB,dw) #collect(-2*tB:dω:2*tB)
#     sum=0
#     for x in steps:
#         sum=sum+quadtest(x,e,hop)*dw
#     return sum
#
#
# even_newer_res(1,3)
#
#
# newres(1,3)
#
# testrange

##########################################################################################################
# Code to update Sigma R and Sigma K matrices


@njit
def sigmaK_init(A):
    # Sigma K matrix setup

    for i in range(testrange):
        A[i,0] = SigmaK(i,0,tB)
        A[0,i] = -np.conj(A[i,0])

    #
    # for i in range(testrange):
    #     A[0,i] = SigmaK(0,i,tB)

    for j in range(testrange):
        for i in range(j,testrange):
            A[i,j] = A[i-j,0]

    for j in range(testrange):
        for i in range(j,testrange):
            A[j,i] = A[0,i-j]


@njit
def sigmaR_init(B):
    # Sigma R matrix setup
    for i in range(testrange):
        B[i,0] =SigmaR(i,0,tB)

    for j in range(testrange):
        for i in range(j,testrange):
            B[i,j] = B[i-j,0]

sigmaR_init(SigmaRmatrix)


############################################################################################################


def simulation():
    sigmaK_init(SigmaKmatrix)
    sigmaR_init(SigmaRmatrix)

    ###################################### Evolution equations start here #############################

    ####### Gr evolution #########

    for i in range(boxinitindex,testrange+1):     ### The diagonal value #should probably start from 2
        # Update GR, GK edges
        for k in range(len(V_ph)):
            for j in range(i+1):
                bessellimit = -complex(0,1)*(Lambda**2)*(1)#*(1)*(1/(tB))*2*tB*(1/2)# *(besselj1(2*tB*abs(t1-t2)*h)/(abs(t1-t2)*h))
                endpoint = (h/2)* Gretarded(k,i+1,i+1)*bessellimit*(h/2)
                endpnt = 1/(1-endpoint)
                Grmatrix[k][i+1,j] = ( complex(0,1)*Gretarded(k,i+1,i)*Grmatrix[k][i,j]+ (h/2)* Gretarded(k,i+1,i)*(Sgma_conv_Gr(k,i,j) + (h/2)*bessellimit*Grmatrix[k][i,j] ) + (h/2)* Gretarded(k,i+1,i+1)*Sgma_conv_Gr(k,i+1,j) )*endpnt





    ######### GK evolution ##############

    for i in range(boxinitindex,testrange+1):     ### The diagonal value #should probably start from 2

        # Update GR, GK edges
        for k in range(len(V_ph)):
            for j in range(i+1):
                bessellimit = -complex(0,1)*(Lambda**2)*(1) #*(1)*(1/(tB))*2*tB*(1/2)# *(besselj1(2*tB*abs(t1-t2)*h)/(abs(t1-t2)*h))
                endpoint = (h/2)* Gretarded(k,i+1,i+1)*bessellimit*(h/2)
                endpnt=1/(1-endpoint)
                Gkmatrix[k][i+1,j] = (complex(0,1)*Gretarded(k,i+1,i)*Gkmatrix[k][i,j]+ (h/2) * Gretarded(k,i+1,i) * ( Sgma_conv_GK(k,i,j)+ h/2*bessellimit*Gkmatrix[k][i,j] + Sgma_conv_GA(k,i,j) )+ (h/2) * Gretarded(k,i+1,i+1) * ( Sgma_conv_GK(k,i+1,j) + Sgma_conv_GA(k,i+1,j) ) )*endpnt
                Gkmatrix[k][j,i+1] = - np.conj(Gkmatrix[k][i+1,j]) # iGᴷ is hermitian  ⟹ iGᴷ(1,2) = conj((iGᴷ(2,1)) ⟹ Gᴷ(1,2) = - conj(Gᴷ(2,1))


        ## Diagonal terms update ##
        #Update GK(t+ϵ,t+ϵ) i.e GK(i+1,i+1) here  - needs Σₑᴿ on the i+1 block edges  i.e.
        for k in range(len(V_ph)):
            bessellimit = -complex(0,1)*(Lambda**2)*(1)#*(1)*(1/(tB))*2*tB*(1/2)# *(besselj1(2*tB*abs(t1-t2)*h)/(abs(t1-t2)*h))
            endpoint = (h/2)* Gretarded(k,i+1,i+1)*bessellimit*(h/2)
            endpnt=1/(1-endpoint)
            Gkmatrix[k][i+1,i+1] = (complex(0,1)*Gretarded(k,i+1,i)*Gkmatrix[k][i,i+1]+ (h/2)*Gretarded(k,i+1,i)*(Sgma_conv_GK(k,i,i+1)+h/2*bessellimit*Gkmatrix[k][i,i+1] + Sgma_conv_GA(k,i,i+1))+ (h/2) * Gretarded(k,i+1,i+1) * ( Sgma_conv_GK(k,i+1,i+1) + Sgma_conv_GA(k,i+1,i+1) ) )*endpnt

    ########################### Evolution equations end here #############################
    #
    # occupations= np.zeros((len(V_ph),testrange),dtype=np.complex128) #Array{ComplexF64,2}(undef,length(V_ph),testrange)
    # for i in range(testrange):
    #     for j in range(len(V_ph)):
    #         occupations[j,i]= (np.imag(Gkmatrix[j][i,i])+1)*0.5
    #
    # newoccu=np.zeros(len(V_ph),dtype=np.int64)
    #
    # # thermal value storage - from bath spectral function
    #
    # for k in range(len(V_ph)):
    #    val = newres(energy_electron(k),tB)[0]
    #    newoccu[k] = val



    ######################
    #Code ends here
    #now printing values onto .csv file
    #######################

    #outputfile=open("Bath bandwidth=%g,Telectron = %g,\n Tbath =%g, mu_bath=%g,mu_electron=%g, time= %g * %g.csv" %(tB,Temp_electron,Temp_bath,mu_bath,mu_electron,h,testrange),"w")

    #outputfile.write("Energy, Occupation \n")
    #outputfile.write("Value from formula")

    #outputfile.write("\n")

#     for k in range(len(V_ph)):
#         row_string="{},{}".format(energyrange[k],newoccu[k])
#         outputfile.write(row_string)
#         outputfile.write("\n")

#     #### CAUTION: CHANGE '/Users/gurukalyanjayasingh/Desktop/' TO respective path on your system ####

#     #pd.DataFrame(occupations).to_csv("/Users/debikalyanjayasingh/Desktop/Temp_guru/Bath_bandwidth=%g,Telectron = %g,\n Tbath =%g, mu_bath=%g,mu_electron=%g, time= %g * %g.csv"%(tB,Temp_electron,Temp_bath,mu_bath,mu_electron,h,testrange))



#     for k in range(len(V_ph)):
#         row_string="{},{}".format(energyrange[k],np.real(occupations[k,testrange-1]))
#         outputfile.write(row_string)
#         outputfile.write("\n")

    #outputfile.write("total time %g" %(Net_time))

    #outputfile.close()



#the big for loop for tB ends here


# In[30]:




jitted_simulation = jit()(simulation)


# In[32]:


get_ipython().run_line_magic('time', 'simulation()')





occupations= np.zeros((len(V_ph),testrange),dtype=np.complex128) #Array{ComplexF64,2}(undef,length(V_ph),testrange)
for i in range(testrange):
    for j in range(len(V_ph)):
        occupations[j,i]= (np.imag(Gkmatrix[j][i,i])+1)*0.5

newoccu=np.zeros(len(V_ph),dtype=np.int64)

# thermal value storage - from bath spectral function

for k in range(len(V_ph)):
   val = newres(energy_electron(k),tB)[0]
   newoccu[k] = val




2+3

get_ipython().run_line_magic('time', 'jitted_simulation()')




def numm(x):
    a = np.arange(1,x)
    val=0
    for i in a:
        val= val+i**2
    return val


# In[5]:



np.arange(1,10)


# In[15]:


get_ipython().run_line_magic('time', 'numm(100000)')


# In[6]:


from numba import jit,njit


# In[28]:


jitted_func = njit()(numm)


# In[44]:


get_ipython().run_line_magic('time', 'jitted_func(10000000)')


# In[53]:


def testing(x):
    val=0
    for i in range(x):
        val=val+ energy_electron(i%len(V_ph))
    return val


# In[173]:


get_ipython().run_line_magic('time', 'testing(10000000)')


# In[59]:


jitted_fur=  njit()(testing)


# In[179]:


get_ipython().run_line_magic('time', 'jitted_fur(10000000)')


# In[180]:


def guru(z):
    return np.tanh(z+1)


# In[181]:


jit_guru = njit()(guru)


# In[182]:


def sequence(x):
    val=0
    for i in range(x):
        val+=guru(i)


# In[187]:


def sequence2(x):
    val=0
    for i in range(x):
        val+=jit_guru(i)


# In[185]:


jit_seq = njit()(sequence)


# In[191]:


jit_seq2 = njit(sequence2)


# In[245]:


get_ipython().run_line_magic('time', 'sequence(10000000)')


# In[246]:


get_ipython().run_line_magic('time', 'sequence2(10000000)')


# In[251]:


get_ipython().run_line_magic('time', 'jit_seq2(10000000)')


# ## Ok! So njit for loops must have njit functions in them. The only choice is then tonake everything njit.

# In[ ]:


@njit
def action():
    car=0
    for i in range(10):
        car = car+i
    return car

action()
