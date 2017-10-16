include("WY94.jl")

R = linspace(0.0, 100.0, 15)
out = zeros((length(R), 5))
for i in eachindex(R)
   K1 = [-R[i], R[i]].+100
   K2 = [-R[i], R[i]].+100
   out[i,1] = R[i]*2
   avg_spe, avg_spe_l = 0.0, 0.0
   avg_gen, avg_gen_l = 0.0, 0.0
   nrep = 10
   for rep in 1:nrep
      o = WY94(50, 50, 50, T=100, K1=K1, K2=K2, linked=false)
      ol = WY94(50, 50, 50, T=100, K1=K1, K2=K2, linked=true)
      m = mean(o[end-100:end,:],1)
      ml = mean(ol[end-100:end,:],1)
      avg_spe += (m[2]+m[3])/2
      avg_gen += m[4]
      avg_spe_l += (ml[2]+ml[3])/2
      avg_gen_l += ml[4]
   end
   out[i,2] = avg_spe/nrep
   out[i,3] = avg_gen/nrep
   out[i,4] = avg_spe_l/nrep
   out[i,5] = avg_gen_l/nrep
end

plot(out[:,1], out[:,2], m=:diamond, c=:black, ls=:dash, lab="Specialists")
plot!(out[:,1], out[:,3], m=:circle, c=:black, ls=:dash, lab="Generalist")
plot!(out[:,1], out[:,4], m=:diamond, c=:black, lab="Specialists", leg=false)
plot!(out[:,1], out[:,5], m=:circle, c=:black, lab="Generalist", leg=false)
xaxis!("Range of variation")
yaxis!("Average density")

savefig("figure03.pdf")
