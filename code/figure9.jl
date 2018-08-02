include("seed.jl")

# Simulation with specialist natural enemies with stochasticity on h
sim, params = simulation(25.0, 8.0, F=4.0, D=0.5, h=5.0, b=25.0, a=0.5, th= 0.0, m=0.5, f=generalist_dyn, h_sd=5)
# Fig 5b)
fig9a = plot(sim[:,1], sim[:,2], label="Hosts", frame=:origin, lw=3, leg=false, ylims=(0, 90))
plot!(fig9a, sim[:,1],sim[:,3], label="Parasites", lw=3)
xlabel!(fig9a, "Generation")
annotate!(fig9a, 5, 85, text("(a)"))
annotate!(fig9a, 25, 80, text("h = 5 ± 5"))

kval = mapslices((r) -> kvalue_by_generation(r, params), sim, 2)
fig9b = plot(log10.(sim[:,2]), kval,
   m=:circle, msc=:white, msw=2, ms=4, mc=:grey,
   lc=:grey, lw=2,
   leg=false, frame=:origin,
   xlims=(-0.5, 2.0),
   ylims=(0, 0.6))
xlabel!(fig9b, "Host density (log 10)")
annotate!(fig9b, -0.25, 0.55, text("(b)"))

sim_R2 = cor(vec(log10.(sim2[:,2])), vec(kval2))^2
sim_a, sim_b = linreg(vec(log10.(sim2[:,2])), vec(kval2))

plot(fig3a, fig9b, layout=(2,1), size=(500,900), margin=5mm)
savefig("article/figures/figure9.pdf")
