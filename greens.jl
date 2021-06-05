
Phonon definitions

### Disperion relation
function ω𝑝(k)
    return t𝑝*abs(sin(k*a1))
end
ω𝑝(10)

### Definition of Bare D_0, Dzerobar, and D_zero_K
function D₀ᴿ(k,t1,t2)
    if t1>t2
        return (-(1)*sin(ω𝑝(k)*(t1-t2)*h))/(ω𝑝(k))     # the equal to case shall give 0
    else
        return 0
    end
end

D₀ᴿ(10,1,2)
-sin(ω𝑝(10)*h)/ω𝑝(10)


function D̄₀ᴿ(k,t1,t2)
    if t1>=t2                            ### What does D̄ do at equal times? produce 1? What if it rigorously doesn't hold?
        return (-1*cos(ω𝑝(k)*(t1-t2)*h))
    else
        return 0
    end                     #remember this is only true if t1>t2
end

D̄₀ᴿ(10,1,1)
cos(ω𝑝(10)*h)


function D₀ᴷ(k,t,t1,Tphonon)
    a= (-im)*(cos(ω𝑝(k)*(t-t1)*h) * coth(ω𝑝(k)*0.5/(Tphonon)) )* (1/ω𝑝(k))
    return a
end

D₀ᴷ(10,2,1,1)

function D̄₀ᴷ(k,t,t1,Tphonon)
    return im*sin(ω𝑝(k)*(t-t1)*h)*coth(ω𝑝(k)*0.5/(Tphonon))
end

D̄₀ᴷ(1,2,1,1)

#%%
Electron Definitions

function ϵₑ(k)
    return tₑ*(1-cos(k*a1))
end


function G₀ᴿ(k,t1,t2)
    if t1>=t2
        return -im*exp(-im*ϵₑ(k)*(t1-t2)*h)
    else
        return 0
    end
end

G₀ᴿ(10,1e4,1e4) #prints 0 for t1<t2

function G₀ᴷ(k,t1,t2,Telectron,μ)
    return -im*tanh((ϵₑ(k)-μ)/(2*Telectron))*exp(-im*ϵ(k)*(t1-t2)*h)
end
