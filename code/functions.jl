"""
TODO
"""
function timestep(N::Float64, P::Float64, p; parasite_dyn=specialist_dyn)
    Nt = host_dyn(N, P, p)
    Pt = parasite_dyn(N, P, p)
    return (Nt, Pt)
end


"""
***Generalist parasitoids population size(Initial population size of female parasitoids) ***

- `h`: Saturation number of parasitoids
- `b`: Rate of appoaching h
- `N`: Initial host population size 

Return : `P`: Generalist parasitoids population size
"""
function generalist_dyn(N::Float64, P::Float64, p)
    Pt = p.h*(1-exp(-N/p.B))
    return Pt
end

"""
***Specialist parasitoids population size(Population size of female parasitoids at next generation) ***

- `c`: Number of parasitoids emerging from each host parasitized
- `N`: Initial host population size 
- `P`: Initial parasitoid population size

Return : `Pt`: Specialist parasitoids population size
"""
function specialist_dyn(N::Float64, P::Float64, p)
    pescape = escape_probability(N,P,p)
    Pt = p.c*N*(1-pescape)
    return Pt
end

"""
***Insect population at next generation***

- `N`: Initial host population size 
- `P`: Initial parasitoid population size
- `F`: Finite rate of increase of the host population
- `D`: Density independent mortality (as a probability of survival)

Return : `Nt`: Population size at the next generation
"""
function host_dyn(N::Float64,P::Float64, p)
    pescape = escape_probability(N,P)
    Nt = p.F*N*pescape*p.D
    return Nt
end

"""
***Probability of escaping mortality from natural ennemies***

- `N`: Initial host population size 
- `P`: Initial parasitoid population size
- `a`: Searching efficiency (per capita)
- `m`: Extent of clumping of the parasitoid attacks
- `th`: Handling time (as a proportion of the total time)

Return : `pescape`: Probability of escaping mortality from natural ennemies
"""
function escape_probability(N::Float64,P::Float64, p)
    num = p.a*P
    den = p.m*(1+p.a*p.th*N)
    pescape = (1+num/den)^-p.m
    return pescape
end

"""
***Per capita searching efficiency***

- `N`: Host population size
- `P`: Parasitoid population size
- `S`: Survivors from parasitism

Return : `A`: Per capita searching efficiency at generation t
"""
function efficiency(N::Float64,P::Float64, S::Float64)
    prop = N/S
    A = 1/P * log(prop)
    return A
end

"""
***Mortality per generation***

- `N`: Host population size
- `S`: Survivors from parasitism

Return : `k`: mortality per generation
"""
function mortality(N::Float64,S::Float64)
    k = log10(N/S)
    return k
end
