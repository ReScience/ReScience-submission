using NamedTuples
using Plots

include("functions.jl")

# Simulation with specialist natural enemies with m=0.2
simulation1 = simulation(50.0, 25.0, m=0.2, F=4.0, D=0.5, c=1.0, a=0.5, th = 0.0)
# Fig 3a)
plot(simulation1[:,2:3], xlabel = "Generations",
                         labels = ["N", "P"],
                         ylims = (0, 50),
                         ylabel = "Populations")

# Simulation with specialist natural enemies with m=0.8
simulation2 = simulation(50.0, 25.0, m=0.8, F=4.0, D=0.5, c=1.0, a=0.5, th = 0.0)
# Fig 3b)
plot(simulation2[:,2:3], xlabel = "Generations",
                         labels = ["N", "P"],
                         ylims = (0, 50),
                         ylabel = "Populations")
