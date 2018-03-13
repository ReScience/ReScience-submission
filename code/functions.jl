# list of functions

"""
***Generalist parasitoids population size(Initial population size of female parasitoids) ***

- `h`: Saturation number of parasitoids
- `b`: Rate of appoaching h
- `N`: Initial host population size 

Return : `P`: Generalist parasitoids population size
"""
function Generalistdynamic(h::Float64,b::Float64, N::Float64)
    P = h*(1-exp(-N/B))
    return P
end

"""
***Specialist parasitoids population size(Population size of female parasitoids at next generation) ***

- `c`: Number of parasitoids emerging from each host parasitized
- `N`: Initial host population size 
- `P`: Initial parasitoid population size

Return : `Pt`: Specialist parasitoids population size
"""
function Specialistdynamic(c::Float64,N::Float64, P::Float64)
    pescape = prob(N,P)
    Pt = c*N*(1-pescape)
    return Pt
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
function prob(N::Float64,P::Float64, a::Float64, m::Float64, th::Float64)
    num = a*P
    den = m*(1+a*th*N)
    pescape = (1+num/den)^-m
    return pescape
end

"""
***Insect population at next generation***

- `N`: Initial host population size 
- `P`: Initial parasitoid population size
- `F`: Finite rate of increase of the host population
- `D`: Density independent mortality (as a probability of survival)

Return : `Nt`: Population size at the next generation
"""
function hostdynamic(N::Float64,P::Float64, F::Float64, D::Float64)
    pescape = prob(N,P)
    Nt = F*N*pescape*D
    return Nt
end
