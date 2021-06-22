using  Debugger


######### Model/simulation parameters
t𝑝=2 #phonon bandwidth
tₑ=1 #electron bandwidth= = σ
a1=1  #lattice constant
λ𝑐= 1 #phonon-electron coupling
Tₑ=0.75# electron temperature
T𝑝=1.75# phonon temperature
μ = -1  # chemical potential of the electron


#time-simulation parameters
h= 0.08 #the time spacing
Time_max = 50 #the net time
N𝑡= Int64(Time_max/h) #


#phonon volume parameters
sitenum = 4 #gives the no. of sites in the lattice
a2=2*π*(1/(sitenum*a1)) #reciprocal space lattice constant
V_ph = collect(-0.5*sitenum*a2:a2:0.5*a2*(sitenum+1))
#filter!(e->e!=0,V_ph) # not taking k=0 mode currently


#%%

#Phonon definitions

### Disperion relation
function ω𝑝(k)
    return t𝑝*abs(sin(V_ph[k]*a1*0.5))+0.2
end


### Definition of Bare D_0, Dzerobar, and D_zero_K
function D₀ᴿ(k,t1,t2)
    if t1>t2
        return (-(1)*sin(ω𝑝(k)*(t1-t2)*h))/(ω𝑝(k))     # the equal to case shall give 0
    else
        return 0
    end
end





function D̄₀ᴿ(k,t1,t2)
    if t1>=t2                            ### What does D̄ do at equal times? produce 1? What if it rigorously doesn't hold?
        return (-1*cos(ω𝑝(k)*(t1-t2)*h))
    else
        return 0
    end                     #remember this is only true if t1>t2
end


function D₀ᴷ(k,t,t1,Tphonon)
    a= (-im)*(cos(ω𝑝(k)*(t-t1)*h) * coth(ω𝑝(k)*0.5/(Tphonon)) )* (1/ω𝑝(k))
    return a
end

D₀ᴷ(1,2,1,T𝑝)

function D̄₀ᴷ(k,t,t1,Tphonon)
    return im*sin(ω𝑝(k)*(t-t1)*h)*coth(ω𝑝(k)*0.5/(Tphonon))
end



#%%

# Electron Definitions

function ϵₑ(k)
    return tₑ*(1-cos(V_ph[k]*a1))+0.2
end


function G₀ᴿ(k,t1,t2)
    if t1>=t2
        return -im*exp(-im*ϵₑ(k)*(t1-t2)*h)
    else
        return 0
    end
end

    #prints 0 for t1<t2

function G₀ᴷ(k,t1,t2,Telectron,μ)

    return -im*tanh((ϵₑ(k)-μ)/(2*Telectron))*exp(-im*ϵₑ(k)*(t1-t2)*h)
end



G₀ᴷ(1,2,3,4,5)


#%%
#Matrix definitions: making Array of arrays : each inner array is 2dim with currently undefine size, the outer array is 1d and holds
#total k points+ 10 elements


    Dᴿmatrix = Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

    D̄ᴿmatrix = Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

    Dᴷmatrix = Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

    D̄ᴷmatrix = Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

    Gᴿmatrix = Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

    Gᴷmatrix = Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

    Σ𝑝ᴿ= Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

    Σ𝑝ᴷ = Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

    Σₑᴿ = Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

    Σₑᴷ = Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

matinit = function ()
    for i=1:length(V_ph)+2
        Dᴿmatrix[i] = Array{ComplexF64,2}(undef,N𝑡+5,N𝑡+5)
        D̄ᴿmatrix[i] = Array{ComplexF64,2}(undef,N𝑡+5,N𝑡+5)
        Dᴷmatrix[i] = Array{ComplexF64,2}(undef,N𝑡+5,N𝑡+5)
        D̄ᴷmatrix[i] = Array{ComplexF64,2}(undef,N𝑡+5,N𝑡+5)
        Gᴿmatrix[i] = Array{ComplexF64,2}(undef,N𝑡+5,N𝑡+5)
        Gᴷmatrix[i] = Array{ComplexF64,2}(undef,N𝑡+5,N𝑡+5)
        Σ𝑝ᴿ[i] = Array{ComplexF64,2}(undef,N𝑡+5,N𝑡+5)
        Σ𝑝ᴷ[i] = Array{ComplexF64,2}(undef,N𝑡+5,N𝑡+5)
        Σₑᴿ[i] = Array{ComplexF64,2}(undef,N𝑡+5,N𝑡+5)
        Σₑᴷ[i] = Array{ComplexF64,2}(undef,N𝑡+5,N𝑡+5)
    end
end
matinit()




##%

### DEFINITIONS OF CONVOLUTION FUNCTIONS

# Definitions of convolutions for phonons

function F(k,t₁,t₂)
    if t₁>t₂
        return sum(t->Σ𝑝ᴿ[k][t₁,t]*Dᴿmatrix[k][t,t₂]*h, collect(t₂:t₁))
    elseif t₁==t₂
        return 0
    else
        return "You're convoluting in the opposite direction. Possible error at RR conv"
    end
end


F(3,9,10)

function RK(k,t1,t2) #∫₀ᵗ Σ𝑝ᴿ⋅Dᴷ
    if t1>1
        return sum(t->Σ𝑝ᴿ[k][t1,t]*Dᴷmatrix[k][t,t2]*h, collect(1:t1))
    else
        return 0
    end
end

RK(2,2,3)

function KA(k,t1,t2) #∫₀⋅Dᴿ
    if t2>1
        return sum(t->Σ𝑝ᴷ[k][t1,t]*Dᴿmatrix[k][t2,t] * h, collect(1:t2)) ## Dᴬ is just transpose of Dᴿ
    else
        return 0
    end
end
KA(2,3,2)

###### Equivalent definitions of convolutions for electrons ###########

function Fₑ(k,t₁,t₂)
    if t₁>t₂
        return sum(t->Σₑᴿ[k][t₁,t]*Gᴿmatrix[k][t,t₂]*h, collect(t₂:t₁))
    elseif t₁==t₂
        return 0
    else
        return "You're convoluting in the opposite direction. Possible error at RR/electron conv"
    end
end

function RKₑ(k,t1,t2) #∫₀ᵗ Σₑᴿ⋅Dᴷ
    if t1>1
        return sum(t->Σₑᴿ[k][t1,t]*Gᴷmatrix[k][t,t2]*h, collect(1:t1))
    else
        return 0
    end
end


RK(2,2,3)

function KAₑ(k,t1,t2) #∫₀⋅Dᴿ
    if t2>1
        return sum(t->Σₑᴷ[k][t1,t]*conj(Gᴿmatrix[k][t2,t] )*h, collect(1:t2))
    else
        return 0
    end
end


function sumBZ1d(k,p) # returns the index of element in the BZ array that reflects the physical sum of two input indices, taking the periodic behaviour into account
    if V_ph[k]*V_ph[p] ==(π/a1)^2
        return p

    elseif -1*π/a1<=V_ph[k]+V_ph[p]<=π/a1
        res = V_ph[k]+V_ph[p]
        return findmin(abs.(V_ph.-res))[2]

    elseif V_ph[k]+V_ph[p]>π/a1
        res = V_ph[k]+V_ph[p]-2*π/a1
        return findmin(abs.(V_ph.-res))[2]

    else V_ph[k]+V_ph[p]<-1*π/a1
        res = V_ph[k]+V_ph[p]+2*π/a1
        return findmin(abs.(V_ph.-res))[2]
    end
end


function negative(k) # returns array index of -k vector
    middle = (length(V_ph)+1)*0.5
    l=length(V_ph)
    return Int((k<middle)*(l-k+1) + (k==middle)*middle+ (k>middle)*(middle-(k-middle)))
end

negative(4)


# lesson - .$ where $ is a  binary operation is the julia equivalent of handling doing array +-* with a scalar on each element
sumBZ1d(1,6) #probably need to define zero mode




#%%

### INITIALIZATIONS


######### diagonal Initialization ##########
boxinitindex=20

boxinit=function()

    for k =1:length(V_ph)
        for i=1:N𝑡
            Dᴿmatrix[k][i,i]=0 #exactly 0
            D̄ᴿmatrix[k][i,i]=-1 #exactly 1
            Gᴿmatrix[k][i,i] = -im #exactly true           ## Gr(t,t)≂̸0
            ##Dᴷmatrix[k][i,i] =1                   ### Why am I initializing the DK?

            ##Gᴷmatrix[k][i,i] = G₀ᴷ(k,t,t,Tₑ,μ)      ### Why am I even initializing this? Aren't we supposed to get this from the code?

            ##Dᴿmatrix[k][i+1,i] = D₀ᴿ(k,i+1,i) #only 2,1 or all i+1,i???
            ##D̄ᴿmatrix[k][i+1,i] = D̄₀ᴿ(k,+i,i)
            ##Σ𝑝ᴿ[k][i+1,i] =
        end
        #println(k)
    end



    ######## Box Initialization ############

    #GF Initialization

    for k=1:length(V_ph)
        for i=1:boxinitindex
            for j=1:boxinitindex
                Dᴿmatrix[k][i,j] = D₀ᴿ(k,i,j)
                D̄ᴿmatrix[k][i,j] = D̄₀ᴿ(k,i,j)
                Dᴷmatrix[k][i,j] = D₀ᴷ(k,i,j,T𝑝)
                D̄ᴷmatrix[k][i,j] = D̄₀ᴷ(k,i,j,T𝑝)
                Gᴿmatrix[k][i,j] = G₀ᴿ(k,i,j)
                Gᴷmatrix[k][i,j] =  G₀ᴷ(k,i,j,Tₑ,μ)

            end
        end
    end
end







        # Dᴿmatrix[k][2,1] = D₀ᴿ(k,2,1) #lower traingular
        # D̄ᴿmatrix[k][2,1] = D̄₀ᴿ(k,2,1) #lower triangular
        #
        # Dᴷmatrix[k][1,1] = D₀ᴷ(k,1,1,T𝑝)
        # Dᴷmatrix[k][1,2] = D₀ᴷ(k,1,2,T𝑝)
        # Dᴷmatrix[k][2,1] = D₀ᴷ(k,2,1,T𝑝)
        # Dᴷmatrix[k][2,2] = D₀ᴷ(k,2,2,T𝑝)
        #
        # D̄ᴷmatrix[k][1,1] = D̄₀ᴷ(k,1,1,T𝑝)
        # D̄ᴷmatrix[k][1,2] = D̄₀ᴷ(k,1,2,T𝑝)
        # D̄ᴷmatrix[k][2,1] = D̄₀ᴷ(k,2,1,T𝑝)
        # D̄ᴷmatrix[k][2,2] = D̄₀ᴷ(k,2,2,T𝑝)
        #
        # Gᴿmatrix[k][2,1] = G₀ᴿ(k,2,1)
        #
        # Gᴷmatrix[k][1,1] =  G₀ᴷ(k,1,1,Tₑ,μ)
        # Gᴷmatrix[k][1,2] =  G₀ᴷ(k,1,2,Tₑ,μ)
        # Gᴷmatrix[k][2,1] =  G₀ᴷ(k,2,1,Tₑ,μ)
        # Gᴷmatrix[k][2,2] =  G₀ᴷ(k,2,2,Tₑ,μ)

boxinit()


#Self energy Initialization


#%%
Dᴿmatrix
matinit()
boxinit()
#Actual for loop
testrange=100
for i=boxinitindex:test       ### The diagonal value #should probably start from 2

    #Update DR
    for k=1:length(V_ph)
        for j=1:i

            if j<i
                D̄ᴿmatrix[k][i,j] = ω𝑝(k)^2 * D₀ᴿ(k,i,i-1) * Dᴿmatrix[k][i-1,j] - D̄₀ᴿ(k,i,i-1) * D̄ᴿmatrix[k][i-1,j] + (h/2)*( D̄₀ᴿ(k,i,i)* F(k,i,j) + D̄₀ᴿ(k,i,i-1) * F(k,i-1,j) )
            end
            @bp
            Dᴿmatrix[k][i+1,j] = D̄₀ᴿ(k,i+1,i) * Dᴿmatrix[k][i,j] + D₀ᴿ(k,i+1,i) * D̄ᴿmatrix[k][i,j] + (h/2)*D₀ᴿ(k,i+1,i)*F(k,i,j)
        end
    end

     #Update DK
     for k = 1:length(V_ph)
         for j=1:i
             D̄ᴷmatrix[k][i,j] = ω𝑝(k)^2 * D₀ᴿ(k,i,i-1) * Dᴷmatrix[k][i-1,j] - D̄₀ᴿ(k,i,i-1) * D̄ᴷmatrix[k][i-1,j] + (h/2)*(  D̄₀ᴿ(k,i,i)* RK(k,i,j) + D̄₀ᴿ(k,i,i-1)* RK(k,i-1,j) + D̄₀ᴿ(k,i,i)* KA(k,i,j) + D̄₀ᴿ(k,i,i-1)* KA(k,i-1,j) )
             Dᴷmatrix[k][i+1,j] = D̄₀ᴿ(k,i+1,i) * Dᴷmatrix[k][i,j] + D₀ᴿ(k,i+1,i) * D̄ᴷmatrix[k][i,j] + (h/2)*( D₀ᴿ(k,i+1,i)* RK(k,i,j) + D₀ᴿ(k,i+1,i)* KA(k,i,j) )
             Dᴷmatrix[k][j,i+1] = -conj(Dᴷmatrix[k][i+1,j])
             D̄ᴷmatrix[k][j,i] = +conj(D̄ᴷmatrix[k][i,j])#what abt i,i entry? If Dk is imaginary, then it will just flip sign here.....? This term is to take care of that...Not sure
         end
     end

    # Update GR, GK
    for k = 1 : length(V_ph)
        for j=1:i
            Gᴿmatrix[k][i+1,j] = im*G₀ᴿ(k,i+1,i)*Gᴿmatrix[k][i,j] + (h/2)* G₀ᴿ(k,i+1,i)*(Fₑ(k,i,j))
            Gᴷmatrix[k][i+1,j] = im*G₀ᴿ(k,i+1,i)*Gᴷmatrix[k][i,j]+ (h/2)*G₀ᴿ(k,i+1,i)* (RKₑ(k,i,j) + KAₑ(k,i,j))
            Gᴷmatrix[k][j,i+1] = - conj(Gᴷmatrix[k][i+1,j]) # iGᴷ is hermitian  ⟹ iGᴷ(1,2) = conj((iGᴷ(2,1)) ⟹ Gᴷ(1,2) = - conj(Gᴷ(2,1))
        end
    end


    # Extract Phonon self energy Σ𝑝ᴷ,Σ𝑝ᴿ in the n+1,n+1 block (use in Dr calculation in next loop)
    # Σ𝑝ᴿ update
    for j=1:i
        for k=1:length(V_ph)
            sum=0
            for p=1:length(V_ph)
                sum = sum -im*(λ𝑐^2/2)*(Gᴿmatrix[sumBZ1d(k,p)][i+1,j]*Gᴷmatrix[p][j,i+1] + Gᴷmatrix[sumBZ1d(k,p)][i+1,j]* conj(Gᴿmatrix[p][i+1,j]) )
            end
            Σ𝑝ᴿ[k][i+1,j]=sum
            sum=nothing
        end
    end

    # Σ𝑝ᴷ Update
    for j=1:i
        for k=1:length(V_ph)
            sum1=0
            sum2=0
            for p=1:length(V_ph)
                ### Code for lower edge
                sum1 = sum1 - im*(λ𝑐^2/2)*( Gᴿmatrix[sumBZ1d(k,p)][i+1,j]*conj(Gᴿmatrix[p][i+1,j]) - Gᴷmatrix[sumBZ1d(k,p)][i+1,j]*conj(Gᴷmatrix[p][i+1,j])  )

                ### Code for upper edge
                sum2 = sum2- im*(λ𝑐^2/2)*( conj(Gᴿmatrix[sumBZ1d(k,p)][i+1,j])*Gᴿmatrix[p][i+1,j] - Gᴷmatrix[sumBZ1d(k,p)][j,i+1]*conj(Gᴷmatrix[p][j,i+1]) )
            end
            Σ𝑝ᴷ[k][i+1,j]= sum1
            Σ𝑝ᴷ[k][j,i+1]= sum2
            sum1=nothing
            sum2=nothing
        end
    end

    #Now extract self energies Σₑᴿ,Σₑᴷ in the n+1,n+1 block (shall be used for calculation of GR in n+2,n+2 loop i.e. next big loop's GR,GK update since you're dropping the self consistent term)

    ## Σₑᴿ Update
    for j=1:i
        for k=1:length(V_ph)
            sum=0
            for p=1:length(V_ph)
                sum = sum + im*( λ𝑐^2/2 )*( Gᴿmatrix[sumBZ1d(k,p)][i+1,j] * Dᴷmatrix[negative(p)][i+1,j] + Gᴷmatrix[sumBZ1d(k,p)][i+1,j] * Dᴿmatrix[negative(p)][i+1,j] )
            end
            Σₑᴿ[k][i+1,j]= sum
            sum=nothing
        end
    end

    ## Σₑᴷ Update
    for j=1:i
        for k=1:length(V_ph)
            sum1=0
            sum2=0
            for p=1:length(V_ph)
                ## code for upper edges
                sum1 = sum1+im*(λ𝑐^2/2)*( Gᴿmatrix[sumBZ1d(k,p)][i+1,j] * Dᴿmatrix[negative(p)][i+1,j] + Gᴷmatrix[sumBZ1d(k,p)][i+1,j] * Dᴷmatrix[negative(p)][i+1,j]  )

                ## code for lower edges
                sum2 = sum2+ im*(λ𝑐^2/2)*( conj( Gᴿmatrix[sumBZ1d(k,p)][i+1,j] )* Dᴿmatrix[p][i+1,j] - conj( Gᴷmatrix[sumBZ1d(k,p)][i+1,j] ) * Dᴷmatrix[p][i+1,j] )
            end
            Σₑᴷ[k][i+1,j]=sum1
            Σₑᴷ[k][j,i+1]=sum2
            sum1=nothing
            sum2=nothing
        end
    end

    ############## Diagonal terms update #############


    #Update GK(t+ϵ,t+ϵ) i.e GK(i+1,i+1) here  - needs Σₑᴿ on the i+1 block edges  i.e.
    for k=1:length(V_ph)
        Gᴷmatrix[k][i+1,i+1] = im*G₀ᴿ(k,i+1,i)*Gᴷmatrix[k][i,i+1]+ (h/2)*G₀ᴿ(k,i+1,i)* (RKₑ(k,i,i+1) + KAₑ(k,i,i+1))
    end
    #initialie Gr(i+1,i+1) = -im from before i.e. outside this loop

    #updating Σ𝑝ᴿ[i+1,i+1],Σ𝑝ᴷ[i+1,i+1]
    for k=1:length(V_ph)
        sum=0
        for p=1:length(V_ph)
            sum = sum -im*(λ𝑐^2/2)*( Gᴿmatrix[sumBZ1d(k,p)][i+1,i+1]*Gᴷmatrix[p][i+1,i+1] + Gᴷmatrix[sumBZ1d(k,p)][i+1,i+1]* conj(Gᴿmatrix[p][i+1,i+1]) )
        end
        Σ𝑝ᴿ[k][i+1,i+1]=sum
        sum=nothing
    end
    #- this will be used to update the DR[i+1,i+1]

    #updating Σ𝑝ᴷ(i+1,i+1)
    for k=1:length(V_ph)
        sum=0
        for p=1:length(V_ph)
            sum = sum - im*(λ𝑐^2/2)*( Gᴿmatrix[sumBZ1d(k,p)][i+1,i+1]*conj(Gᴿmatrix[p][i+1,i+1]) - Gᴷmatrix[sumBZ1d(k,p)][i+1,i+1]*conj(Gᴷmatrix[p][i+1,i+1])  )
        end
        Σ𝑝ᴷ[k][i+1,i+1] = sum
        sum=nothing
    end


    #Update DK(t+ϵ,t+ϵ) here, D̄(i,i) block is calculated already
    for k=1:length(V_ph)
        D̄ᴷmatrix[k][i+1,i] = ω𝑝(k)^2 * D₀ᴿ(k,i+1,i) * Dᴷmatrix[k][i,i] - D̄₀ᴿ(k,i+1,i) * D̄ᴷmatrix[k][i,i] + (h/2)*(  D̄₀ᴿ(k,i+1,i+1)* RK(k,i+1,i) + D̄₀ᴿ(k,i+1,i)* RK(k,i,i) + D̄₀ᴿ(k,i+1,i+1)* KA(k,i+1,i) + D̄₀ᴿ(k,i+1,i)* KA(k,i,i) )
        D̄ᴷmatrix[k][i,i+1] = conj( D̄ᴷmatrix[k][i+1,i] )
        Dᴷmatrix[k][i+1,i+1] = D̄₀ᴿ(k,i+1,i) * Dᴷmatrix[k][i,i+1] + D₀ᴿ(k,i+1,i) * D̄ᴷmatrix[k][i,i+1] + (h/2)*( D₀ᴿ(k,i+1,i)* RK(k,i,i+1) + D₀ᴿ(k,i+1,i)* KA(k,i,i+1) )
    end

    #updating sigmElR, sigmaElK here (given that now you've all Gfuncs at the i+1,i+1)
    for k=1:length(V_ph)
        sum=0
        for p=1:length(V_ph)
            sum = sum + im*( λ𝑐^2/2 )*( Gᴿmatrix[sumBZ1d(k,p)][i+1,i+1] * Dᴷmatrix[negative(p)][i+1,i+1] + Gᴷmatrix[sumBZ1d(k,p)][i+1,i+1] * Dᴿmatrix[negative(p)][i+1,i+1] )
        end
        Σₑᴿ[k][i+1,i+1]= sum
        sum=nothing
    end

    for k=1:length(V_ph)
        sum=0
        for p=1:length(V_ph)
            sum = sum+im*(λ𝑐^2/2)*( 0 + Gᴷmatrix[sumBZ1d(k,p)][i+1,i+1] * Dᴷmatrix[negative(p)][i+1,i+1]  )
        end
        Σₑᴷ[k][i+1,i+1]=sum
        sum=nothing
    end

end


#%%

#Plotting

set
using Plots
Tₑ
T𝑝
testk=1
testrange =80
#a,b = Array{ComplexF64}(undef,testrange),Array{ComplexF64}(undef,testrange)
b2=[]
a2=[]
for m=1:length(V_ph)
    a= Array{Float64}(undef,testrange)
    b = Array{Float64}(undef,testrange)
    for i=1:testrange
        #a[i] = Dᴿmatrix[testk][i,1]
        #b[i] = Gᴿmatrix[testk][i,1]
        a[i] =  (real(im*Dᴷmatrix[m][i,i])*ω𝑝(m)-1)*0.5            #real(im*Dᴷmatrix[testk][i,i])*ω𝑝(testk)
        b[i] = (imag(Gᴷmatrix[m][i,i])+1 )*0.5
    end
    push!(b2,b)
    push!(a2,a)
end


#real(b2)
plot(real(b2),lw= 3,title = "Tₑ = $(Tₑ), Tphonon = $(T𝑝)")



savefig(figpath*"testing_order_of_filling1")
plot(real(a2),legend = false)

t= collect(h:h:testrange*h)

plot(t,real(b),label="nₑ(k,t,t) electrons")
plot!(t,real(a),label="n𝑝(k,t,t) phonons", title= "For Tₑ=$(Tₑ), T_phonon = $(T𝑝), # of modes=$(sitenum)")

#plot!(t,imag(b2),label="nₑ(k,t,t) electrons")


#plot(t,real(a),label="phonons" )



figpath = "/Users/gurukalyanjayasingh/Desktop/Julia/Plots/"

savefig(figpath*"fig4.png")

for m=1:length(b)
    println(1+imag(b[m]))
end



#Disperion plots
disp_phonon=[]
disp_electron= []
for m=1:length(V_ph)
    push!(disp_phonon,ω𝑝(m))
    push!(disp_electron,ϵₑ(m))
end

plot(disp_phonon)
plot(disp_electron)
savefig(figpath*"dispersion2.png")
#%% Code ends

Tₑ
T𝑝

arr1 =[]

for m=1:length(V_ph)
    push!(arr1, im*G₀ᴷ(m,1,1,Tₑ,μ))
end
arr1

plot(arr1)

#%%
Debugging

matinit()
boxinit()
f4 = function()
    for i=2:12  ### The diagonal value #should probably start from 2



        #Update DR
        for k=1:length(V_ph)
            for j=1:i

                if j<i
                    D̄ᴿmatrix[k][i,j] = ω𝑝(k)^2 * D₀ᴿ(k,i,i-1) * Dᴿmatrix[k][i-1,j] - D̄₀ᴿ(k,i,i-1) * D̄ᴿmatrix[k][i-1,j] + (h/2)*( D̄₀ᴿ(k,i,i)* F(k,i,j) + D̄₀ᴿ(k,i,i-1) * F(k,i-1,j) )
                end
                Dᴿmatrix[k][i+1,j] = D̄₀ᴿ(k,i+1,i) * Dᴿmatrix[k][i,j] + D₀ᴿ(k,i+1,i) * D̄ᴿmatrix[k][i,j] + (h/2)*D₀ᴿ(k,i+1,i)*F(k,i,j)
            end
        end

         #Update DK
         for k = 1:length(V_ph)
             for j=1:i
                 D̄ᴷmatrix[k][i,j] = ω𝑝(k)^2 * D₀ᴿ(k,i,i-1) * Dᴷmatrix[k][i-1,j] - D̄₀ᴿ(k,i,i-1) * D̄ᴷmatrix[k][i-1,j] + (h/2)*(  D̄₀ᴿ(k,i,i)* RK(k,i,j) + D̄₀ᴿ(k,i,i-1)* RK(k,i-1,j) + D̄₀ᴿ(k,i,i)* KA(k,i,j) + D̄₀ᴿ(k,i,i-1)* KA(k,i-1,j) )
                 Dᴷmatrix[k][i+1,j] = D̄₀ᴿ(k,i+1,i) * Dᴷmatrix[k][i,j] + D₀ᴿ(k,i+1,i) * D̄ᴷmatrix[k][i,j] + (h/2)*( D₀ᴿ(k,i+1,i)* RK(k,i,j) + D₀ᴿ(k,i+1,i)* KA(k,i,j) )
                 Dᴷmatrix[k][j,i+1] = -conj(Dᴷmatrix[k][i+1,j])
                 D̄ᴷmatrix[k][j,i] = +conj(D̄ᴷmatrix[k][i,j])#what abt i,i entry? If Dk is imaginary, then it will just flip sign here.....? This term is to take care of that...Not sure
             end
         end

        # Update GR, GK
        for k = 1 : length(V_ph)
            for j=1:i
                Gᴿmatrix[k][i+1,j] = im*G₀ᴿ(k,i+1,i)*Gᴿmatrix[k][i,j] + (h/2)* G₀ᴿ(k,i+1,i)*(Fₑ(k,i,j))
                Gᴷmatrix[k][i+1,j] = im*G₀ᴿ(k,i+1,i)*Gᴷmatrix[k][i,j]+ (h/2)*G₀ᴿ(k,i+1,i)* (RKₑ(k,i,j) + KAₑ(k,i,j))
                Gᴷmatrix[k][j,i+1] = - conj(Gᴷmatrix[k][i+1,j]) # iGᴷ is hermitian  ⟹ iGᴷ(1,2) = conj((iGᴷ(2,1)) ⟹ Gᴷ(1,2) = - conj(Gᴷ(2,1))
            end
        end


        # Extract Phonon self energy Σ𝑝ᴷ,Σ𝑝ᴿ in the n+1,n+1 block (use in Dr calculation in next loop)
        # Σ𝑝ᴿ update
        for j=1:i
            for k=1:length(V_ph)
                sum=0
                for p=1:length(V_ph)
                    sum = sum -im*(λ𝑐^2/2)*(Gᴿmatrix[sumBZ1d(k,p)][i+1,j]*Gᴷmatrix[p][j,i+1] + Gᴷmatrix[sumBZ1d(k,p)][i+1,j]* conj(Gᴿmatrix[p][i+1,j]) )
                end
                Σ𝑝ᴿ[k][i+1,j]=sum
                sum=nothing
            end
        end

        # Σ𝑝ᴷ Update
        for j=1:i
            for k=1:length(V_ph)
                sum1=0
                sum2=0
                for p=1:length(V_ph)
                    ### Code for lower edge
                    sum1 = sum1 - im*(λ𝑐^2/2)*( Gᴿmatrix[sumBZ1d(k,p)][i+1,j]*conj(Gᴿmatrix[p][i+1,j]) - Gᴷmatrix[sumBZ1d(k,p)][i+1,j]*conj(Gᴷmatrix[p][i+1,j])  )

                    ### Code for upper edge
                    sum2 = sum2- im*(λ𝑐^2/2)*( conj(Gᴿmatrix[sumBZ1d(k,p)][i+1,j])*Gᴿmatrix[p][i+1,j] - Gᴷmatrix[sumBZ1d(k,p)][j,i+1]*conj(Gᴷmatrix[p][j,i+1]) )
                end
                Σ𝑝ᴷ[k][i+1,j]= sum1
                Σ𝑝ᴷ[k][j,i+1]= sum2
                sum1=nothing
                sum2=nothing
            end
        end

        #Now extract self energies Σₑᴿ,Σₑᴷ in the n+1,n+1 block (shall be used for calculation of GR in n+2,n+2 loop i.e. next big loop's GR,GK update since you're dropping the self consistent term)

        ## Σₑᴿ Update
        for j=1:i
            for k=1:length(V_ph)
                sum=0
                for p=1:length(V_ph)
                    sum = sum + im*( λ𝑐^2/2 )*( Gᴿmatrix[sumBZ1d(k,p)][i+1,j] * Dᴷmatrix[negative(p)][i+1,j] + Gᴷmatrix[sumBZ1d(k,p)][i+1,j] * Dᴿmatrix[negative(p)][i+1,j] )
                end
                Σₑᴿ[k][i+1,j]= sum
                sum=nothing
            end
        end

        ## Σₑᴷ Update
        for j=1:i
            for k=1:length(V_ph)
                sum1=0
                sum2=0
                for p=1:length(V_ph)
                    ## code for upper edges
                    sum1 = sum1+im*(λ𝑐^2/2)*( Gᴿmatrix[sumBZ1d(k,p)][i+1,j] * Dᴿmatrix[negative(p)][i+1,j] + Gᴷmatrix[sumBZ1d(k,p)][i+1,j] * Dᴷmatrix[negative(p)][i+1,j]  )

                    ## code for lower edges
                    sum2 = sum2+ im*(λ𝑐^2/2)*( conj( Gᴿmatrix[sumBZ1d(k,p)][i+1,j] )* Dᴿmatrix[p][i+1,j] - conj( Gᴷmatrix[sumBZ1d(k,p)][i+1,j] ) * Dᴷmatrix[p][i+1,j] )
                end
                Σₑᴷ[k][i+1,j]=sum1
                Σₑᴷ[k][j,i+1]=sum2
                sum1=nothing
                sum2=nothing
            end
        end

        ############## Diagonal terms update #############

        if i==10
            @bp
        end
        #Update GK(t+ϵ,t+ϵ) i.e GK(i+1,i+1) here  - needs Σₑᴿ on the i+1 block edges  i.e.
        for k=1:length(V_ph)
            Gᴷmatrix[k][i+1,i+1] = im*G₀ᴿ(k,i+1,i)*Gᴷmatrix[k][i,i+1]+ (h/2)*G₀ᴿ(k,i+1,i)* (RKₑ(k,i,i+1) + KAₑ(k,i,i+1))
        end
        #initialie Gr(i+1,i+1) = -im from before i.e. outside this loop

        #updating Σ𝑝ᴿ[i+1,i+1],Σ𝑝ᴷ[i+1,i+1]
        for k=1:length(V_ph)
            sum=0
            for p=1:length(V_ph)
                sum = sum -im*(λ𝑐^2/2)*( Gᴿmatrix[sumBZ1d(k,p)][i+1,i+1]*Gᴷmatrix[p][i+1,i+1] + Gᴷmatrix[sumBZ1d(k,p)][i+1,i+1]* conj(Gᴿmatrix[p][i+1,i+1]) )
            end
            Σ𝑝ᴿ[k][i+1,i+1]=sum
            sum=nothing
        end
        #- this will be used to update the DR[i+1,i+1]

        #updating Σ𝑝ᴷ(i+1,i+1)
        for k=1:length(V_ph)
            sum=0
            for p=1:length(V_ph)
                sum = sum - im*(λ𝑐^2/2)*( Gᴿmatrix[sumBZ1d(k,p)][i+1,i+1]*conj(Gᴿmatrix[p][i+1,i+1]) - Gᴷmatrix[sumBZ1d(k,p)][i+1,i+1]*conj(Gᴷmatrix[p][i+1,i+1])  )
            end
            Σ𝑝ᴷ[k][i+1,i+1] = sum
            sum=nothing
        end

        if i==10
            @bp
        end
        #Update DK(t+ϵ,t+ϵ) here, D̄(i,i) block is calculated already
        for k=1:length(V_ph)
            D̄ᴷmatrix[k][i+1,i] = ω𝑝(k)^2 * D₀ᴿ(k,i+1,i) * Dᴷmatrix[k][i,i] - D̄₀ᴿ(k,i+1,i) * D̄ᴷmatrix[k][i,i] + (h/2)*(  D̄₀ᴿ(k,i+1,i+1)* RK(k,i+1,i) + D̄₀ᴿ(k,i+1,i)* RK(k,i,i) + D̄₀ᴿ(k,i+1,i+1)* KA(k,i+1,i) + D̄₀ᴿ(k,i+1,i)* KA(k,i,i) )
            D̄ᴷmatrix[k][i,i+1] = conj( D̄ᴷmatrix[k][i+1,i] )
            Dᴷmatrix[k][i+1,i+1] = D̄₀ᴿ(k,i+1,i) * Dᴷmatrix[k][i,i+1] + D₀ᴿ(k,i+1,i) * D̄ᴷmatrix[k][i,i+1] + (h/2)*( D₀ᴿ(k,i+1,i)* RK(k,i,i+1) + D₀ᴿ(k,i+1,i)* KA(k,i,i+1) )
        end

        #updating sigmElR, sigmaElK here (given that now you've all Gfuncs at the i+1,i+1)
        for k=1:length(V_ph)
            sum=0
            for p=1:length(V_ph)
                sum = sum + im*( λ𝑐^2/2 )*( Gᴿmatrix[sumBZ1d(k,p)][i+1,i+1] * Dᴷmatrix[negative(p)][i+1,i+1] + Gᴷmatrix[sumBZ1d(k,p)][i+1,i+1] * Dᴿmatrix[negative(p)][i+1,i+1] )
            end
            Σₑᴿ[k][i+1,i+1]= sum
            sum=nothing
        end

        for k=1:length(V_ph)
            sum=0
            for p=1:length(V_ph)
                sum = sum+im*(λ𝑐^2/2)*( 0 + Gᴷmatrix[sumBZ1d(k,p)][i+1,i+1] * Dᴷmatrix[negative(p)][i+1,i+1]  )
            end
            Σₑᴷ[k][i+1,i+1]=sum
            sum=nothing
        end
        println(i)
    end
end


@run f4()


g(x) = x^3

f(x) = g(x)+x

f(2)
