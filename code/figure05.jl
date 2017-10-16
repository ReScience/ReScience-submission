include("WY94.jl")

k=[10,190]
ks=[50,150]
B = linspace(0.3, 1.0, 12)
out = zeros((length(B), 5))
for i in eachindex(B)
   b = B[i]
   out[i,1] = b
   avg_spe, avg_spe_s = 0.0, 0.0
   avg_gen, avg_gen_s = 0.0, 0.0
   nrep = 10
   for rep in 1:nrep
      o = WY94(50, 50, 50, b=b, T=100, K1=k, K2=k, linked=false)
      os = WY94(50, 50, 50, b=b, T=100, K1=ks, K2=ks, linked=false)
      m = mean(o[end-100:end,:],1)
      ms = mean(os[end-100:end,:],1)
      avg_spe += (m[2]+m[3])/2
      avg_gen += m[4]
      avg_spe_s += (ms[2]+ms[3])/2
      avg_gen_s += ms[4]
   end
   out[i,2] = avg_spe/nrep
   out[i,3] = avg_gen/nrep
   out[i,4] = avg_spe_s/nrep
   out[i,5] = avg_gen_s/nrep
end

plot(out[:,1], out[:,2], m=:diamond, c=:black, lab="Specialists")
plot!(out[:,1], out[:,3], m=:circle, c=:black, lab="Generalist")
plot!(out[:,1], out[:,4], m=:diamond, c=:black, ls=:dash, lab="Specialists", leg=false)
plot!(out[:,1], out[:,5], m=:circle, c=:black, ls=:dash, lab="Generalist", leg=false)
xaxis!("Ability of generalist")
yaxis!("Average density")

savefig("figure05.pdf")
