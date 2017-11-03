include("WY94.jl")

k2 = 100
K2=[0.5,1.5].*k2
K = logspace(2,3,8)
out = zeros((length(K), 4))
for i in eachindex(K)
   K1 = [0.5,1.5].*K[i]
   out[i,1] = K[i]/k2
   avg_spe_1, avg_spe_2, avg_gen = 0.0, 0.0, 0.0
   nrep = 10
   for rep in 1:nrep
      o = WY94(50, 50, 50, T=100, K1=K1, K2=K2, linked=false)
      m = mean(o[end-100:end,:],1)
      avg_spe_1 += m[2]
      avg_spe_2 += m[3]
      avg_gen += m[4]
   end
   out[i,2] = avg_spe_1/nrep
   out[i,3] = avg_spe_2/nrep
   out[i,4] = avg_gen/nrep
end

plot(out[:,1], out[:,2], m=:utriangle, c=:black, lab="Specialist 1")
plot!(out[:,1], out[:,3], m=:dtriangle, c=:black, lab="Specialist 2")
plot!(out[:,1], out[:,4], m=:circle, c=:black, lab="Generalist", leg=false)
xaxis!("Ratio of habitat quality", :log10)
yaxis!("Average density", :log10, (10, 1000))

savefig("figure04.pdf")
