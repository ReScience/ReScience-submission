using NamedTuples
using Plots
include("functions.jl")
include("seed.jl")

# Simulation with generalist natural enemies
sim, params = simulation(50.0, 25.0, m=0.5, F=4.0, D=0.5, h=10.0, b=25.0, a=0.5,
 th= 0.0, m=0.5, f=generalist_dyn)
# Fig 4a)
fig4a = plot(sim[:,1],sim[:,2], label="Hosts", frame=:origin, lw=3)
plot!(fig4a, sim[:,1],sim[:,3], label="Parasites", lw=3)
xlabel!(fig4a, "Generation")
ylabel!(fig4a, "Population size")

kval = mapslices((r) -> kvalue_by_generation(r, params), sim, 2)
fig4b = plot(log10.(sim[:,2]), kval,
   m=:circle, msc=:white, msw=2, ms=4, mc=:grey,
   lc=:grey, lw=2,
   leg=false, frame=:origin,
   xlims=(-0.5, 2.0),
   ylims=(0, 1.5))
xlabel!(fig4b, "Host density (log 10)")
ylabel!(fig4b, "k-value")

plot(fig4a, fig4b, layout=(2,1), size=(500,900))
savefig("article/figures/figure4.pdf")
