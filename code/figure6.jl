using NamedTuples
using Plots
include("functions.jl")
# Simulation with specialist natural enemies with stochasticity on D
sim1, params1 = simulation(50.0, 25.0, m=0.5, F=4.0, D=0.5, h=10.0, b=25.0, a=0.5, th= 0.0, m=0.5, f=generalist_dyn, D_sd=0.1)
# Fig 5a)
fig6a = plot(sim1[:,1],sim1[:,2], label="Hosts", frame=:origin, lw=3)
plot!(fig6a, sim1[:,1],sim1[:,3], label="Parasites", lw=3)
xlabel!(fig6a, "Generation")
ylabel!(fig6a, "Population size")

# Simulation with specialist natural enemies with stochasticity on c
sim2, params2 = simulation(50.0, 25.0, m=0.5, F=4.0, D=0.5, h=10.0, b=25.0, a=0.5, th= 0.0, m=0.5, f=generalist_dyn, h_sd=5)
# Fig 5b)
fig6b = plot(sim2[:,1],sim2[:,2], label="Hosts", frame=:origin, lw=3)
plot!(fig6b, sim2[:,1],sim2[:,3], label="Parasites", lw=3)
xlabel!(fig6b, "Generation")
ylabel!(fig6b, "Population size")

# Simulation with specialist natural enemies with stochasticity on a
sim3, params3 = simulation(50.0, 25.0, m=0.5, F=4.0, D=0.5, h=10.0, b=25.0, a=0.5, th= 0.0, m=0.5, f=generalist_dyn, a_sd=0.1)
# Fig 5c)
fig6c = plot(sim3[:,1],sim3[:,2], label="Hosts", frame=:origin, lw=3)
plot!(fig6c, sim3[:,1],sim3[:,3], label="Parasites", lw=3)
xlabel!(fig6c, "Generation")
ylabel!(fig6c, "Population size")

plot(fig6a, fig6b, fig6c, layout=(1,3), size=(900,900))
