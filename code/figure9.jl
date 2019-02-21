include("seed.jl")

# Simulation with specialist natural enemies with stochasticity on h
sim, params = simulation(25.0, 8.0; m=0.5, h=5.0, f=generalist_dyn, h_sd=5)
# Fig 5b)
fig9a = plot(sim[:,1], sim[:,2], label="Hosts", frame=:origin, lw=2, leg=false, ylims=(0, 90), xlims=(0,50))
plot!(fig9a, sim[:,1],sim[:,3], label="Parasites", lw=2)
xlabel!(fig9a, "Generation")
annotate!(fig9a, generate_legend_position((0,50),(0,90); shift=0.4)..., text("h = 5 ± 5"))
annotate!(fig9a, generate_legend_position((0,50),(0,90))..., text("(a)"))

kval = mapslices((r) -> kvalue_by_generation(r, params), sim; dims=2)
fig9b = plot(log10.(sim[:,2]), kval,
   m=:circle, msc=:white, msw=2, ms=4, mc=:grey,
   lc=:grey, lw=2,
   leg=false, frame=:box,
   xlims=(1.0, 2.0),
   ylims=(0, 0.75))
xlabel!(fig9b, "Host density (log 10)")
annotate!(fig9b, generate_legend_position((1,2),(0,0.75))..., text("(b)"))

sim_R2 = cor(vec(log10.(sim[:,2])), vec(kval))^2
sim_a, sim_b = linreg(vec(log10.(sim[:,2])), vec(kval))

plot(fig9a, fig9b, layout=(2,1), size=(500,900), margin=5mm)
savefig("../article/figures/figure9.pdf")
