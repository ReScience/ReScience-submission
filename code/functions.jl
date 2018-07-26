"""
This function (used internally by simulation) returns a random value of x, drawn
from a truncated normal (between 0 and 100).
"""
function generate(x, σ)
    @assert σ >= 0.0
    return σ > 0.0 ? rand(TruncatedNormal(x, σ, 0.0, 100.0)) : x
end

"""
***Simulations for the host population and specialist parasitoid population dynamics for the first 50 generations***

- `N`: Initial host population size
- `P`: Initial parasitoid population size
- `t`: Number of generations
- `f`: Type of parasite chosen to to the simulations (specialist or generalist)
- `F`: Finite rate of increase of the host population
- `D`: Density independent mortality (as a probability of survival)
- `c`: Number of parasitoids emerging from each host parasitized
- `a`: Searching efficiency (per capita)
- `h`: Saturation number of parasitoids
- `b`: Rate of appoaching h
- `th`: Handling time (as a proportion of the total time)
- `m`: Extent of clumping of the parasitoid attacks

Return : `dynamics`: Matrix containing the host and parasite density for the initial populations and for every generation simulated. Time is rows.
"""
function simulation(N::Float64, P::Float64; t::Int64=50, f=specialist_dyn, F=4.0, D=0.5, c=1.0, a=0.5, h=10.0, b=25.0, th=0.0, m=0.2, D_std=0.0, a_sd=0.0, h_sd=0.0, c_sd=0.0)
    # Matrix to store the output
    dynamics = zeros(Float64, (t+2,7))
    dynamics[1,2] = N
    dynamics[1,3] = P
    dynamics[1,4:7] = [D, a, h, c]
    # Iterations
    for current_time in 1:(t+1)
        current_D=generate(D, D_std)
        current_a=generate(a, a_sd)
        current_h=generate(h, h_sd)
        current_c=generate(c, c_sd)
        current_p = @NT(F=F, D=current_D, c=current_c, a=current_a, h=current_h, b=b, th=th, m=m)
        N_next, P_next = timestep(dynamics[current_time,2], dynamics[current_time,3], current_p; parasite_dyn=f)
        # Population sizes
        dynamics[current_time+1,1] = current_time
        dynamics[current_time+1,2] = N_next
        dynamics[current_time+1,3] = P_next
        # Real parameter values
        dynamics[current_time,4] = current_D
        dynamics[current_time,5] = current_a
        dynamics[current_time,6] = current_h
        dynamics[current_time,7] = current_c
    end
    p = @NT(F=F, D=D, c=c, a=a, h=h, b=b, th=th, m=m)
    # Return
    return dynamics[1:(end-1),:], p
end

"""
This function updates the population sizes internally during simulation.
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

Return : `Pt`: Generalist parasitoids population size
"""
function generalist_dyn(N::Float64, P::Float64, p)
    Pt = p.h*(1-exp(-N/p.b))
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
    pescape = escape_probability(N,P,p)
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
***Mortality per generation***

- `N`: Host population size
- `P`: Parasitoid population size
- `p`: parameters list

Return : `kvalue`: mortality per generation
"""
function mortality(N::Float64,P::Float64, p)
    pescape = escape_probability(N,P,p)
    S = pescape * N
    kvalue = log10(N/S)
    return kvalue
end

"""
***Per capita searching efficiency***

- `N`: Host population size
- `P`: Parasitoid population size
- `p`: parameters list

Return : `A`: Per capita searching efficiency at generation t
"""
function efficiency(N::Float64,P::Float64, p)
    pescape = escape_probability(N,P,p)
    S = pescape * N
    prop = N/S
    A = 1.0/P * log(prop)
    return A
end

function kvalue_by_generation(x, p)
   t, N, P, D, a, h, c = x
   current_p = @NT(F=p.F, D=D, c=c, a=a, h=h, b=p.b, th=p.th, m=p.m)
   return mortality(N, P, current_p)
end
