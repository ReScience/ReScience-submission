using NamedTuples
using Plots
include("functions.jl")
include("seed.jl")

# Simulation with specialist natural enemies with stochasticity on D
sim1, params1 = simulation(50.0, 25.0, m=0.5, F=4.0, D=0.5, h=10.0, b=25.0, a=0.5, th= 0.0, m=0.5, f=generalist_dyn, D_sd=0.1)
# Fig 5a)
fig6a = plot(sim1[:,1],sim1[:,2], label="Hosts", frame=:origin, lw=3)
plot!(fig6a, sim1[:,1],sim1[:,3], label="Parasites", lw=3)
xlabel!(fig6a, "Generation")
ylabel!(fig6a, "Population size")

kval1 = mapslices((r) -> kvalue_by_generation(r, params1), sim1, 2)
fig6d = plot(log10.(sim1[:,2]), kval1,
   m=:circle, msc=:white, msw=2, ms=4, mc=:grey,
   lc=:grey, lw=2,
   leg=false, frame=:origin,
   xlims=(-0.5, 2.0),
   ylims=(0, 1.5))
xlabel!(fig6d, "Host density (log 10)")
ylabel!(fig6d, "k-value")

# Simulation with specialist natural enemies with stochasticity on c
sim2, params2 = simulation(50.0, 25.0, m=0.5, F=4.0, D=0.5, h=10.0, b=25.0, a=0.5, th= 0.0, m=0.5, f=generalist_dyn, h_sd=5)
# Fig 5b)
fig6b = plot(sim2[:,1],sim2[:,2], label="Hosts", frame=:origin, lw=3)
plot!(fig6b, sim2[:,1],sim2[:,3], label="Parasites", lw=3)
xlabel!(fig6b, "Generation")
ylabel!(fig6b, "Population size")

kval2 = mapslices((r) -> kvalue_by_generation(r, params2), sim2, 2)
fig6e = plot(log10.(sim2[:,2]), kval2,
   m=:circle, msc=:white, msw=2, ms=4, mc=:grey,
   lc=:grey, lw=2,
   leg=false, frame=:origin,
   xlims=(-0.5, 2.0),
   ylims=(0, 1.5))
xlabel!(fig6e, "Host density (log 10)")
ylabel!(fig6e, "k-value")

# Simulation with specialist natural enemies with stochasticity on a
sim3, params3 = simulation(50.0, 25.0, m=0.5, F=4.0, D=0.5, h=10.0, b=25.0, a=0.5, th= 0.0, m=0.5, f=generalist_dyn, a_sd=0.1)
# Fig 5c)
fig6c = plot(sim3[:,1],sim3[:,2], label="Hosts", frame=:origin, lw=3)
plot!(fig6c, sim3[:,1],sim3[:,3], label="Parasites", lw=3)
xlabel!(fig6c, "Generation")
ylabel!(fig6c, "Population size")

kval3 = mapslices((r) -> kvalue_by_generation(r, params3), sim3, 2)
fig6f = plot(log10.(sim3[:,2]), kval3,
   m=:circle, msc=:white, msw=2, ms=4, mc=:grey,
   lc=:grey, lw=2,
   leg=false, frame=:origin,
   xlims=(-0.5, 2.0),
   ylims=(0, 1.5))
xlabel!(fig6f, "Host density (log 10)")
ylabel!(fig6f, "k-value")

plot(fig6a, fig6b, fig6c, fig6d, fig6e, fig6f, layout=(2,3), size=(1200,900))
savefig("figure6.pdf")
