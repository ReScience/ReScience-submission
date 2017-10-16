include("WY94.jl")

o1 = WY94(0, 0, 3, T=15)
o2 = WY94(3, 3, 270, T=76)[1:5:71,:]

p1d = plot(o1[:,1], o1[:,4], m=:circle, c=:black, leg=false)
p1p = plot(o1[:,1], o1[:,end], m=:circle, c=:black, leg=false)
yaxis!(p1p, [0, 1], "Habitat preference")
yaxis!(p1d, [0, 270], "Density")
xaxis!(p1p, "Time")

p2d = plot(o2[:,1], o2[:,4], m=:circle, c=:black, leg=false)
plot!(p2d, o2[:,1], o2[:,2], m=:utriangle, c=:black, leg=false)
plot!(p2d, o2[:,1], o2[:,3], m=:dtriangle, c=:black, leg=false)

p2p = plot(o2[:,1], o2[:,end], m=:circle, c=:black, leg=false)
plot!(p2p, o2[:,1], o2[:,5], m=:utriangle, c=:black, leg=false)
plot!(p2p, o2[:,1], o2[:,6], m=:dtriangle, c=:black, leg=false)

yaxis!(p2p, [0, 1])
yaxis!(p2d, [0, 270])
xaxis!(p2p, "Time")

plot(p1d, p2d, p1p, p2p)

savefig("figure02.pdf")
