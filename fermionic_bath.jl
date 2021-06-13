## Code to test for thermalisation using the fermionic sector of the Non-self_consistent code

## Try to use maximum parts from the non-self consistent code

######### Model/simulation parameters
using SpecialFunctions
using Plots
using QuadGK
using BenchmarkTools

tₑ=0.2 #electron bandwidth= = σ
tB =-1 #bath bandwidth
a1=1  #lattice constant
λ= 1 #system bath coupling
Tₑ=0# electron temperature
Temp_bath =0 #the bath temperature
μ = -1  # chemical potential of the electron
μbath = -1 #chem potential of the bath

#time-simulation parameters
h= 0.08 #the time spacing
Time_max = 50 #the net time
N𝑡= Int64(Time_max/h) #


#volume parameters
sitenum = 20 #gives the no. of sites in the lattice
a2=2*π*(1/(sitenum*a1)) #reciprocal space lattice constant
V_ph = collect(-0.5*sitenum*a2:a2:0.5*a2*(sitenum+1))




#%% Electron Definitions


function ϵₑ(k)
    return tₑ*(V_ph[k]^2/2)+0.2              #(1-cos(V_ph[k]*a1))
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


##% Convolution function definition

function Fₑ(k,t₁,t₂)
    if t₁>t₂
        return sum(t->Σₑᴿ[t₁,t]*Gᴿmatrix[k][t,t₂]*h, collect(t₂:t₁))
    elseif t₁==t₂
        return 0
    else
        return "You're convoluting in the opposite direction. Possible error at RR/electron conv"
    end
end

function RKₑ(k,t1,t2) #∫₀ᵗ Σₑᴿ⋅Dᴷ
    if t1>1
        return sum(t->Σₑᴿ[t1,t]*Gᴷmatrix[k][t,t2]*h, collect(1:t1))
    else
        return 0
    end
end

function KAₑ(k,t1,t2) #∫₀⋅Dᴿ
    if t2>1
        return sum( t->Σₑᴷ[t1,t]*conj(Gᴿmatrix[k][t2,t] )*h, collect(1:t2) )
    else
        return 0
    end
end


#%% Bath Functions

# what's my convention of fourier transform??????
#SETTLE THIS

Σᴿ(t1,t2) =-im*(λ^2)*2*π*(1/tB)*besselj1(2*tB*abs(t1-t2)*h)/(abs(t1-t2)*h)*(t1>t2)
J(ω) = (2/tB)*sqrt( 1- (ω/(2*tB))^2 )
#x=collect(-2*tB:0.01:2*tB)
#plot(x,J.(x).*tanh.(x.-μ))
#Σᴿ(10,1)

σᴷ(ω) = -im*(λ^2)*J(ω)*tanh((ω-μbath)/(2*Temp_bath))


function Σᴷ(t1,t2)
    dω=(1/500)*4*tB
    steps = collect(-2*tB:dω:2*tB)
    #filter!(e->e!=0,steps)
    result=0
    for ω in steps
        result = result + dω*σᴷ(ω)*exp(-im*ω*(t1-t2)*h)
    end
    return (result/(2*π))
end

#%% Matrix Initializations

Gᴿmatrix = Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

Gᴷmatrix = Array{Array{ComplexF64,2},1}(undef,length(V_ph)+2)

Σₑᴿ = Array{ComplexF64,2}(undef,N𝑡+5,N𝑡+5)

Σₑᴷ = Array{ComplexF64,2}(undef,N𝑡+5,N𝑡+5)

#%% Initialization


matinit = function ()
    for i=1:length(V_ph)+2
        Gᴿmatrix[i] = Array{ComplexF64,2}(undef,N𝑡+5,N𝑡+5)
        Gᴷmatrix[i] = Array{ComplexF64,2}(undef,N𝑡+5,N𝑡+5)
    end
end


matinit()


### Code to update the Sigma R, Sigma K matrix

for i=1:N𝑡
    Σₑᴿ[i,1] = Σᴿ(i,1)
end

for j=2:N𝑡
    for i=j:N𝑡
        Σₑᴿ[i,j] = Σₑᴿ[i-j+1,1]
    end
end


for i=1:N𝑡
    Σₑᴷ[i,1] = Σᴷ(i,1)
end

for i=1:N𝑡
    Σₑᴷ[1,i] = Σᴷ(1,i)
end

for j=1:N𝑡
    for i=j:N𝑡
        Σₑᴷ[i,j] = Σₑᴷ[i-j+1,1]
    end
end

for j=1:N𝑡
    for i=j:N𝑡
        Σₑᴷ[j,i] = Σₑᴷ[1,i-j+1]
    end
end




Tₑ
Temp_bath




######### G Initialization ##########
boxinitindex=1

boxinit=function()

    for k =1:length(V_ph)
        for i=1:N𝑡
            Gᴿmatrix[k][i,i] = -im #exactly true           ## Gr(t,t)≂̸0
        end
    end



    ######## Box Initialization ############

    #GF Initialization

    for k=1:length(V_ph)
        for i=1:boxinitindex
            for j=1:boxinitindex
                Gᴿmatrix[k][i,j] = G₀ᴿ(k,i,j)
                Gᴷmatrix[k][i,j] =  G₀ᴷ(k,i,j,Tₑ,μ)

            end
        end
    end

end

boxinit()



#%% Evolution equations

matinit()
boxinit()


N𝑡
testrange =100

####### Gr evolution #########
for i=boxinitindex:testrange     ### The diagonal value #should probably start from 2
    # Update GR, GK edges
    for k = 1 : length(V_ph)
        for j=1:i
            Gᴿmatrix[k][i+1,j] = im*G₀ᴿ(k,i+1,i)*Gᴿmatrix[k][i,j] + (h/2)* G₀ᴿ(k,i+1,i)*(Fₑ(k,i,j)) + (h/2)* G₀ᴿ(k,i+1,i+1)*(Fₑ(k,i+1,j))
        end
    end
end

######### GK evolution ##############

for i=boxinitindex:testrange     ### The diagonal value #should probably start from 2

    # Update GR, GK edges
    for k = 1 : length(V_ph)
        for j=1:i
            Gᴷmatrix[k][i+1,j] = im*G₀ᴿ(k,i+1,i)*Gᴷmatrix[k][i,j]+ (h/2)*G₀ᴿ(k,i+1,i)* (RKₑ(k,i,j) + KAₑ(k,i,j))
            Gᴷmatrix[k][j,i+1] = - conj(Gᴷmatrix[k][i+1,j]) # iGᴷ is hermitian  ⟹ iGᴷ(1,2) = conj((iGᴷ(2,1)) ⟹ Gᴷ(1,2) = - conj(Gᴷ(2,1))
        end
    end

    ############## Diagonal terms update #############
    #Update GK(t+ϵ,t+ϵ) i.e GK(i+1,i+1) here  - needs Σₑᴿ on the i+1 block edges  i.e.
    for k=1:length(V_ph)
        Gᴷmatrix[k][i+1,i+1] = im*G₀ᴿ(k,i+1,i)*Gᴷmatrix[k][i,i+1]+ (h/2)*G₀ᴿ(k,i+1,i)* (RKₑ(k,i,i+1) + KAₑ(k,i,i+1))
    end
end


#%% Testing

b2=[]
boxinitindex
testrange
for m=1:length(V_ph)
    b = Array{ComplexF64}(undef,testrange)
    for i=1:testrange
        #b[i] = Gᴿmatrix[m][i,1]
        b[i] = Gᴷmatrix[m][i,i]#(imag(Gᴷmatrix[m][i,i])+1)*0.5
        #b[i] = Gᴷmatrix[m][i,i]
        #println(testrange)
    end
    push!(b2,b)
    println(m)
end

b2
ser = collect(1:testrange)
plot(ser,real.(b2),legend=false,title = "Tₑ = $(Tₑ), Tbath = $(Temp_bath), μbath = $(μbath), μelectron = $(μ)",lw= 3)#lw= 3


disp_electron= []
for m=1:length(V_ph)
    push!(disp_electron,ϵₑ(m))
end
disp_electron

x=collect(-π/a1:a2:π/a1)
scatter(x,disp_electron,label="system")
scatter!(x, 2*-tB*cos.(x),title = "Dispersion of system & bath, Te = $(Tₑ), Tbath = $(Temp_bath)",label="Bath")
plot!(collect(-π/a1:1e-2:π/a1),μ.+collect(-π/a1:1e-2:π/a1).*0,label = "μ_electron = $(μ)")
plot!(collect(-π/a1:1e-2:π/a1),μbath.+collect(-π/a1:1e-2:π/a1).*0,label = "μ_bath = $(μbath)")


b3=nothing

b3=[]

for m=1:length(V_ph)
    push!(b3,G₀ᴷ(m,1,1,0,μ))
end

plot(imag.(b3))

using PyPlot
