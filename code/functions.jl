# list of functions

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
