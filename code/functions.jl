# list of functions

"""
***Generalist parasitoids population size(Initial population size of female parasitoids) ***

- `h`: Saturation number of parasitoids
- `b`: Rate of appoaching h
- `N`: Initial host population size 

Return : `P`: Generalist parasitoids population size
"""
function Generalistdynamic(h::Float64,b::Float64, N::Float64)
    P = h*(1-e^(-N/B))
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
