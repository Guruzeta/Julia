import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
import scipy.special as spl
import scipy.integrate as integrate
import pandas as pd

pie=np.pi

tBlist=np.array([8]) #enter the list of bath bandwidths that you want to simulate for

h= 0.1 # enter the timestep to be used for simulation
Time_max = 10
Net_time=int(np.floor(Time_max/h)) #
Lambda= 1 #system bath coupling
t_electron=1 #electron bandwidth= = σ
#tB =6 #bath bandwidth
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
fermi = lambda e,T,mu: 1/(np.exp((e-mu)*(1/T))+1)
energy_electron = lambda k : t_electron*(1-np.cos(V_ph[k]*a1))
Gretarded = lambda k,t1,t2: (t1>=t2)*complex(0,-1)*np.exp(complex(0,-1)*energy_electron(k)*(t1-t2)*h)
Gkeldysh =lambda k,t1,t2,Telectron,mu: -complex(0,1)*np.tanh((energy_electron(k)-mu)/(2*Telectron))*np.exp(-complex(0,1)*energy_electron(k)*(t1-t2)*h)



#bath parameters
J = lambda w,tB: (2/tB)*np.sqrt( 1- (w/(2*tB))**2 )

def SigmaR(t1,t2,tB):
    if t1>t2:
        sum1= -complex(0,1)*(Lambda**2)*(1)*(1/(tB))*( spl.jv(1,2*tB*abs(t1-t2)*h) / (abs(t1-t2)*h ) )
        return sum1
    else:
        return 0

sigmak=lambda w,tB: -complex(0,1)*(Lambda**2)*J(w,tB)*np.tanh((w-mu_bath)/(2*Temp_bath))
def SigmaK(t1,t2,tB):
    dw=(1/10000)*4*tB
    steps = np.arange(-2*tB,2*tB,dw) #collect(-2*tB:dω:2*tB)
    result=0
    for w in steps:
        result = result + dw*sigmak(w,tB)*np.exp(-complex(0,1)*w*(t1-t2)*h)

    return result/(2*pie)




#matrix definitions

Grmatrix =np.zeros(len(V_ph)+2,dtype='O')      # Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

Gkmatrix =np.zeros(len(V_ph)+2,dtype='O')  #Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

SigmaRmatrix = np.zeros(shape=(Net_time+5,Net_time+5),dtype=complex)#Array{ComplexF64,2}(undef,Net_time+5,Net_time+5)

SigmaKmatrix = np.zeros(shape=(Net_time+5,Net_time+5),dtype=complex)#Array{ComplexF64,2}(undef,Net_time+5,Net_time+5)







### spawning Gr & Gk matrices
def matinit():
    for i in range(len(V_ph)+2):
        Grmatrix[i] = np.zeros(shape=(Net_time+5,Net_time+5),dtype=complex)#Array{ComplexF64,2}(undef,Net_time+5,Net_time+5)
        Gkmatrix[i] = np.zeros(shape=(Net_time+5,Net_time+5),dtype=complex)


## Box Initialization ###
boxinitindex=1

def boxinit():

    for k in range(len(V_ph)):
        for i in range(Net_time):
            Grmatrix[k][i,i] = -complex(0,1) #exactly true           ## Gr(t,t)≂̸0

                ###actual Box Initialization ###

    #GF Initialization

    for k in range(len(V_ph)):
        for i in range(boxinitindex+1):
            for j in range(boxinitindex+1):
                Grmatrix[k][i,j] = Gretarded(k,i,j)
                Gkmatrix[k][i,j] =  Gkeldysh(k,i,j,Temp_electron,mu_electron)


#### Convolution definitions

def Sgma_conv_Gr(k,t1,t2):
    if t1>t2:
        sum=0
        sum=sum + SigmaRmatrix[t1,t2]*Grmatrix[k][t2,t2]*(h/2)

        for i in range(t2+1,t1):
            sum = sum + SigmaRmatrix[t1,i]*Grmatrix[k][i,t2]*h
        return sum

    else:
        return 0

def Sgma_conv_GK(k,t1,t2):
    if t1>1:
        sum=0
        sum = sum + SigmaRmatrix[t1,1]*Gkmatrix[k][1,t2]*(h/2)

        for i in range(2,t1):               #2:t1-1
            sum = sum+SigmaRmatrix[t1,i]*Gkmatrix[k][i,t2]*h

        return sum

    else:
        return 0


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

energyrange=[]
for k in range(len(V_ph)):
    energyrange.append(energy_electron(k))

quadtest =lambda x,en,tB: (1/pie)*fermi(x,Temp_bath,mu_bath)*(Lambda**2)*(J(x,tB))*( 1/( (x-en)**2 + (Lambda**2*J(x,tB))**2 ) )

def newres(e,hop):  # computes n(E) i.e. equilibrium occupation at energy E using bath spectral function
    return integrate.quad(quadtest, -2*hop,2*hop,args=(e,hop))





for tB in tBlist:


    ## Convolution function definition

    # def Sgma_conv_Gr(k,t1,t2):
    #     if t1>t2:
    #         sum=0
    #         sum=sum + SigmaRmatrix[t1,t2]*Grmatrix[k][t2,t2]*(h/2)
    #
    #         for i in range(t2+1,t1):
    #             sum = sum + SigmaRmatrix[t1,i]*Grmatrix[k][i,t2]*h
    #         return sum
    #
    #     else:
    #         return 0
    #
    # def Sgma_conv_GK(k,t1,t2):
    #     if t1>1:
    #         sum=0
    #         sum = sum + SigmaRmatrix[t1,1]*Gkmatrix[k][1,t2]*(h/2)
    #
    #         for i in range(2,t1):               #2:t1-1
    #             sum = sum+SigmaRmatrix[t1,i]*Gkmatrix[k][i,t2]*h
    #
    #         return sum
    #
    #     else:
    #         return 0
    #
    #
    # def Sgma_conv_GA(k,t1,t2):
    #     if t2>1:
    #         sum=0
    #         sum=sum+SigmaKmatrix[t1,1]* np.conjugate(Grmatrix[k][t2,1]) * h*(1/2)  #starting 1/2
    #         sum=sum+ SigmaKmatrix[t1,t2]* np.conjugate(Grmatrix[k][t2,t2]) * h*(1/2) #ending 1/2
    #
    #         for i in range(2,t2):# i=2:t2-1:
    #             sum=sum+ SigmaKmatrix[t1,i]* np.conjugate(Grmatrix[k][t2,i]) * h   #middle ones, they get h & not h/2 due to double addition
    #
    #         return sum
    #
    #     else:
    #         return 0

    #### Bath definitions


    #J = lambda w: (2/tB)*np.sqrt( 1- (w/(2*tB))**2 )

    # def SigmaR(t1,t2,tB):
    #     if t1>t2:
    #         sum1= -complex(0,1)*(Lambda**2)*(1)*(1/(tB))*( spl.jv(1,2*tB*abs(t1-t2)*h) / (abs(t1-t2)*h ) )
    #         return sum1
    #     else:
    #         return 0

    #sigmak=lambda w: -complex(0,1)*(Lambda**2)*J(w,tB)*np.tanh((w-mu_bath)/(2*Temp_bath))

    # def SigmaK(t1,t2,tB):
    #     dw=(1/10000)*4*tB
    #     steps = np.arange(-2*tB,2*tB,dw) #collect(-2*tB:dω:2*tB)
    #     result=0
    #     for w in steps:
    #         result = result + dw*sigmak(w,tB)*np.exp(-complex(0,1)*w*(t1-t2)*h)
    #
    #     return result/(2*pie)





    # def matinit():
    #     for i in range(len(V_ph)+2):
    #         Grmatrix[i] = np.zeros(shape=(Net_time+5,Net_time+5),dtype=complex)#Array{ComplexF64,2}(undef,Net_time+5,Net_time+5)
    #         Gkmatrix[i] = np.zeros(shape=(Net_time+5,Net_time+5),dtype=complex)




    # ######### Box Initialization ##########
    # boxinitindex=1

    # def boxinit():
    #
    #     for k in range(len(V_ph)):
    #         for i in range(Net_time):
    #             Grmatrix[k][i,i] = -complex(0,1) #exactly true           ## Gr(t,t)≂̸0
    #
    # ######## Box Initialization ############
    #
    #     #GF Initialization
    #
    #     for k in range(len(V_ph)):
    #         for i in range(boxinitindex+1):
    #             for j in range(boxinitindex+1):
    #                 Grmatrix[k][i,j] = Gretarded(k,i,j)
    #                 Gkmatrix[k][i,j] =  Gkeldysh(k,i,j,Temp_electron,mu_electron)




    #totaltime = Net_time
    testrange =Net_time
    matinit() #spawn Gr matrix, Gk matrix
    boxinit()

    ### Code to update the Sigma R, Sigma K matrix


    # Sigma R matrix setup
    for i in range(testrange+1):
        SigmaRmatrix[i,0] =SigmaR(i,0,tB)


    for j in range(testrange+1):
        for i in range(j,testrange+1):
            SigmaRmatrix[i,j] = SigmaRmatrix[i-j,0]


    # Sigma K matrix setup

    for i in range(testrange+1):
        SigmaKmatrix[i,0] = SigmaK(i,0,tB)


    for i in range(testrange+1):
        SigmaKmatrix[0,i] = SigmaK(0,i,tB)


    for j in range(testrange+1):
        for i in range(j,testrange+1):
            SigmaKmatrix[i,j] = SigmaKmatrix[i-j,0]


    for j in range(testrange+1):
        for i in range(j,testrange+1):
            SigmaKmatrix[j,i] = SigmaKmatrix[0,i-j]

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



    # energyrange=[]
    # for k in range(len(V_ph)):
    #     energyrange.append(energy_electron(k))



    occupations= np.zeros((len(V_ph),testrange),dtype=complex) #Array{ComplexF64,2}(undef,length(V_ph),testrange)
    for i in range(testrange):
        for j in range(len(V_ph)):
            occupations[j,i]= (np.imag(Gkmatrix[j][i,i])+1)*0.5



    # quadtest =lambda x,en,tB: (1/pie)*fermi(x,Temp_bath,mu_bath)*(Lambda**2)*(J(x,tB))*( 1/( (x-en)**2 + (Lambda**2*J(x,tB))**2 ) )
    #
    # def newres(e,hop):  # computes n(E) i.e. equilibrium occupation at energy E using bath spectral function
    #     return integrate.quad(quadtest, -2*hop,2*hop,args=(e,hop))


    newoccu=[]

    # thermal value storage - from bath spectral function

    for k in range(len(V_ph)):
       val = newres(energy_electron(k),tB)[0]
       newoccu.append(val)



    ######################
    #Code ends here
    #now printing values onto .csv file
    #######################

    outputfile=open("Bath bandwidth=%g,Telectron = %g,\n Tbath =%g, mu_bath=%g,mu_electron=%g, time= %g * %g.csv" %(tB,Temp_electron,Temp_bath,mu_bath,mu_electron,h,testrange),"w")

    outputfile.write("Energy, Occupation \n")
    outputfile.write("Value from formula")

    outputfile.write("\n")

    for k in range(len(V_ph)):
        row_string="{},{}".format(energyrange[k],newoccu[k])
        outputfile.write(row_string)
        outputfile.write("\n")

    #### CAUTION: CHANGE '/Users/gurukalyanjayasingh/Desktop/' TO respective path on your system ####

    pd.DataFrame(occupations).to_csv("/Users/debikalyanjayasingh/Desktop/Temp_guru/Bath_bandwidth=%g,Telectron = %g,\n Tbath =%g, mu_bath=%g,mu_electron=%g, time= %g * %g.csv"%(tB,Temp_electron,Temp_bath,mu_bath,mu_electron,h,testrange))



    for k in range(len(V_ph)):
        row_string="{},{}".format(energyrange[k],np.real(occupations[k,testrange-1]))
        outputfile.write(row_string)
        outputfile.write("\n")

    outputfile.write("total time %g" %(Net_time))

    outputfile.close()



    plt.scatter(energyrange,newoccu,label="using form of bath J(w) directly", marker="*")
    plt.scatter(energyrange,np.real(occupations[:,testrange-1]),label="From code")
    plt.xlabel("Energies")
    plt.legend()
    plt.ylabel("Occupations nk")
    plt.title("Bath bandwidth=%g,Temp_electron = %g,\n Tbath =%g, mu_bath=%g,mu_electron=%g, Stepsize= %g\n totaltime = %g" %(tB,Temp_electron,Temp_bath,mu_bath,mu_electron,h, h*Net_time))
    

#the big for loop for tB ends here
