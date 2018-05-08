include("WY94.jl")

k=[50,150]
G = linspace(0.0, 0.3, 12)
out = zeros((length(G), 5))
for i in eachindex(G)
   g = G[i]
   out[i,1] = g
   avg_spe = 0.0
   avg_gen = 0.0
   nrep = 10
   for rep in 1:nrep
      o = WY94(50, 50, 50, g=g, T=100, K1=k, K2=k)
      m = mean(o[end-100:end,:],1)
      avg_spe += (m[2]+m[3])/2
      avg_gen += m[4]
   end
   out[i,2] = avg_spe/nrep
   out[i,3] = avg_gen/nrep
end

plot(out[:,1], out[:,2], m=:diamond, c=:black, lab="Specialists", leg=false)
plot!(out[:,1], out[:,3], m=:circle, c=:black, lab="Generalist", leg=false)
xaxis!("Proportion of non-adaptive choice")
yaxis!("Average density")

savefig("figure06.pdf")
