using NamedTuples
using Plots
include("functions.jl")
# Simulation with specialist natural enemies with stochasticity on D
sim1, params1 = simulation(50.0, 25.0, m=0.2, F=4.0, D=0.5, c=1.0, a=0.5, th= 0.0, m=0.5, f=specialist_dyn, D_sd=0.1)
# Fig 5a)
fig5a = plot(sim1[:,1],sim1[:,2], label="Hosts", frame=:origin, lw=3)
plot!(fig5a, sim1[:,1],sim1[:,3], label="Parasites", lw=3)
xlabel!(fig5a, "Generation")
ylabel!(fig5a, "Population size")

# Simulation with specialist natural enemies with stochasticity on c
sim2, params2 = simulation(50.0, 25.0, m=0.2, F=4.0, D=0.5, c=1.0, a=0.5, th= 0.0, m=0.5, f=specialist_dyn, c_sd=0.1)
# Fig 5b)
fig5b = plot(sim2[:,1],sim2[:,2], label="Hosts", frame=:origin, lw=3)
plot!(fig5b, sim2[:,1],sim2[:,3], label="Parasites", lw=3)
xlabel!(fig5b, "Generation")
ylabel!(fig5b, "Population size")

# Simulation with specialist natural enemies with stochasticity on a
sim3, params3 = simulation(50.0, 25.0, m=0.2, F=4.0, D=0.5, c=1.0, a=0.5, th= 0.0, m=0.5, f=specialist_dyn, a_sd=0.1)
# Fig 5c)
fig5c = plot(sim3[:,1],sim3[:,2], label="Hosts", frame=:origin, lw=3)
plot!(fig5c, sim3[:,1],sim3[:,3], label="Parasites", lw=3)
xlabel!(fig5c, "Generation")
ylabel!(fig5c, "Population size")

plot(fig5a, fig5b, fig5c, layout=(1,3), size=(900,900))
