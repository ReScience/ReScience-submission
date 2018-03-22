using NamedTuples
using Plots

include("functions.jl")

function kvalue_by_generation(x, p)
   t, N, P = x
   return mortality(N, P, p)
end

# Simulation with specialist natural enemies with m=0.2
sim1, params1 = simulation(50.0, 25.0, m=0.2, F=4.0, D=0.5, c=1.0, a=0.5, th = 0.0)
# Fig 3a)
fig3a = plot(sim1[:,1],sim1[:,2], label="Hosts", frame=:origin, lw=3)
plot!(fig3a, sim1[:,1],sim1[:,3], label="Parasites", lw=3)
xlabel!(fig3a, "Generation")
ylabel!(fig3a, "Population size")

kval1 = mapslices((r) -> kvalue_by_generation(r, params1), sim1, 2)
fig3c = plot(log10.(sim1[:,2]), kval1,
   m=:circle, msc=:white, msw=2, ms=4, mc=:grey,
   lc=:grey, lw=2,
   leg=false, frame=:origin,
   xlims=(-0.5, 2.0),
   ylims=(0, 1.5))
xlabel!(fig3c, "Host density (log 10)")
ylabel!(fig3c, "k-value")


# Simulation with specialist natural enemies with m=0.8
sim2, params2 = simulation(50.0, 25.0, m=0.8, F=4.0, D=0.5, c=1.0, a=0.5, th = 0.0)
# Fig 3b)
fig3b = plot(sim2[:,1],sim2[:,2], label="Hosts", frame=:origin, lw=3, leg=false)
plot!(fig3b, sim2[:,1],sim2[:,3], label="Parasites", lw=3)
xlabel!(fig3b, " ")
ylabel!(fig3b, " ")

kval2 = mapslices((r) -> kvalue_by_generation(r, params2), sim2, 2)
fig3d = plot(log10.(sim2[:,2]), kval2,
   m=:circle, msc=:white, msw=2, ms=4, mc=:grey,
   lc=:grey, lw=2,
   leg=false, frame=:origin,
   xlims=(-0.5, 2.0),
   ylims=(0, 1.5))
xlabel!(fig3d, " ")
ylabel!(fig3d, " ")

plot(fig3a, fig3b, fig3c, fig3d, layout=(2,2), size=(900,900))
savefig("figure3.pdf")
