# list of functions
"""
***Insect population at next generation***

- `N`: Initial host population size 
- `P`: Initial parasitoid population size
- `F`: Finite rate of increase of the host population
- `D`: Density independent mortality (as a probability of survival)

Return : `Nt`: Population size at the next generation
"""

function hostdynamic(N::Float64,P::Float64, F::Float64, D::Float64)
    pescape = prop(N,P)
    Nt = F*N*pescape*D
    return Nt
end
