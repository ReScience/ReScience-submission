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

Return : `dynamics`: Matrix containing the host and parasite density for the initial populations and for every generation simulated
"""
function simulation(N::Float64, P::Float64; t::Int64=50, f=specialist_dyn, F=4.0, D=0.5, c=1.0, a=0.5, h=10.0, b=25.0, th=0.0, m=0.2, D_sd=0.0, a_sd=0.0, h_sd=0.0, c_sd=0.0)
    # Matrix to store the output
    dynamics = zeros(Float64, (t+1,3))
    dynamics[1,2] = N
    dynamics[1,3] = P
    # Iterations
    for current_time in 1:t
        current_D=randn()*D_sd+D
        current_a=randn()*a_sd+a
        current_h=randn()*h_sd+h
        current_c=randn()*c_sd+c
        current_p = @NT(F=F, D=current_D, c=current_c, a=current_a, h=current_h, b=b, th=th, m=m)
        N_next, P_next = timestep(dynamics[current_time,2], dynamics[current_time,3], current_p; parasite_dyn=f)
        dynamics[current_time+1,1] = current_time
        dynamics[current_time+1,2] = N_next
        dynamics[current_time+1,3] = P_next
        end
    p = @NT(F=F, D=D, c=c, a=a, h=h, b=b, th=th, m=m)
    # Return
    return dynamics, p
end

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
