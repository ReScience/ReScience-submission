include("WY94.jl")

gs = 25

b = linspace(0.3, 0.9, gs)
kvar = linspace(0.0, 100.0, gs)

o = zeros(Float64, (gs, gs))

for i in eachindex(b)
    for j in eachindex(kvar)
        hstore = zeros(10)
        k = [100.0, 100.0].+[-kvar[j], kvar[j]]
        for re in 1:length(hstore)
            s = WY94(3,3,3; K1=k, K2=k, b=b[i])
            hstore[re] = mean(pielou(s[:,2],s[:,3],s[:,4]))
        end
        o[i,j] = mean(hstore)
    end
end

plot(2.0.*kvar, b, o, lt=:surface, c=:viridis, leg=false, zlim=[0,1], zlab="Pielou's evenness")
xaxis!("Range of variation")
yaxis!("Ability of generalist")

savefig("figureSurf.pdf")
