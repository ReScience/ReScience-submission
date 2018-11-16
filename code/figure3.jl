include("seed.jl")

sim1, params1 = simulation(50.0, 25.0; m=0.2)
# Fig 3a)
fig3a = plot(sim1[:,1], sim1[:,2], label="Hosts", frame=:origin, lw=3)
plot!(fig3a, sim1[:,1],sim1[:,3], label="Parasites", lw=3)
xlabel!(fig3a, "Generation")
ylabel!(fig3a, "Population size")
#annotate!(fig3a, 5, 45, text("(a)"))

kval1 = mapslices((r) -> kvalue_by_generation(r, params1), sim1; dims=2)
fig3c = plot(log10.(sim1[:,2]), kval1,
   m=:circle, msc=:white, msw=2, ms=4, mc=:grey,
   lc=:grey, lw=2,
   leg=false, frame=:box,
   xlims=(1.2, 1.8),
   ylims=(0.2, 0.4))
xlabel!(fig3c, "Host density (log 10)")
ylabel!(fig3c, "k-value")
#annotate!(fig3c, -0.25, 1.35, text("(c)"))

sim2, params2 = simulation(50.0, 25.0; m=0.8)
# Fig 3b)
fig3b = plot(sim2[:,1], sim2[:,2], label="Hosts", frame=:origin, lw=3, leg=false)
plot!(fig3b, sim2[:,1],sim2[:,3], label="Parasites", lw=3)
xlabel!(fig3b, " ")
ylabel!(fig3b, " ")
#annotate!(fig3b, 5, 45, text("(b)"))

kval2 = mapslices((r) -> kvalue_by_generation(r, params2), sim2; dims= 2)
fig3d = plot(log10.(sim2[:,2]), kval2,
   m=:circle, msc=:white, msw=2, ms=4, mc=:grey,
   lc=:grey, lw=2,
   leg=false, frame=:box,
   xlims=(-0.5, 2.0),
   ylims=(0.0, 1.2))
xlabel!(fig3d, " ")
ylabel!(fig3d, " ")
#annotate!(fig3d, -0.25, 1.35, text("(d)"))

plot(fig3a, fig3b, fig3c, fig3d, ann = (:top_left, :auto), layout=(2,2), size=(900,900), margin=5mm)
savefig("article/figures/figure3.pdf")
