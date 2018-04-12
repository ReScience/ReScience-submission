include("seed.jl")

N = 5000

cor_a = zeros(Float64, N)
cor_D = zeros(Float64, N)
cor_c = zeros(Float64, N)

for i in eachindex(cor_a)
    sd = i == 1 ? 0.0 : 0.5
    sim, params = simulation(50.0, 25.0, m=0.2, F=4.0, D=0.5, c=1.0, a=0.5, th= 0.0, m=0.5, f=specialist_dyn, a_sd=sd)
    kval = mapslices((r) -> kvalue_by_generation(r, params), sim, 2)
    correlation = cor(vec(log10.(sim[:,2])), vec(kval))
    cor_a[i] = correlation
end

for i in eachindex(cor_D)
    sd = i == 1 ? 0.0 : 0.5
    sim, params = simulation(50.0, 25.0, m=0.2, F=4.0, D=0.5, c=1.0, a=0.5, th= 0.0, m=0.5, f=specialist_dyn, D_sd=sd)
    kval = mapslices((r) -> kvalue_by_generation(r, params), sim, 2)
    correlation = cor(vec(log10.(sim[:,2])), vec(kval))
    cor_D[i] = correlation
end

for i in eachindex(cor_c)
    sd = i == 1 ? 0.0 : 0.5
    sim, params = simulation(50.0, 25.0, m=0.2, F=4.0, D=0.5, c=1.0, a=0.5, th= 0.0, m=0.5, f=specialist_dyn, c_sd=sd)
    kval = mapslices((r) -> kvalue_by_generation(r, params), sim, 2)
    correlation = cor(vec(log10.(sim[:,2])), vec(kval))
    cor_c[i] = correlation
end

pl_a = density(cor_a[2:end], xlim=(0,1.0), ylim=(0,10), frame=:origin, c=:black, fill=(0, :grey, 0.2), leg=false)
vline!(pl_a, [first(cor_a)], lw=2, ls=:dot, c=:grey)

pl_D = density(cor_D[2:end], xlim=(0,1.0), ylim=(0,10), frame=:origin, c=:black, fill=(0, :grey, 0.2), leg=false)
vline!(pl_D, [first(cor_D)], lw=2, ls=:dot, c=:grey)
ylabel!(pl_D, "Density")

pl_c = density(cor_c[2:end], xlim=(0,1.0), ylim=(0,10), frame=:origin, c=:black, fill=(0, :grey, 0.2), leg=false)
vline!(pl_c, [first(cor_c)], lw=2, ls=:dot, c=:grey)
xlabel!(pl_c, "Correlation")

plot(pl_D, pl_c, pl_a, layout=(1,3), size=(1200,400), margin=5mm)
savefig("article/figures/figure7.pdf")
