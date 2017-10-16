using Plots
using StatsBase

pyplot()

function habitat_preference(N,K;trials=50)
   p = [1.0, 0.0, rand()] # Because why not?
   TN = sum(N)
   for iteration in 1:trials
      l, m, n = sample(1:3, 3, replace=false)
      K1, K2 = K[l,1], K[l,2]
      NUM = (K1+K2)*N[n]*p[n] + (K1+K2)*N[m]*p[m] - K1*N[l] - K1*N[m] - K1*N[n]
      DEN = (K1+K2)*N[l]
      RATIO = minimum([maximum([-(NUM/DEN), 0.0]), 1.0])
      p[l] = DEN > 0.0 ? RATIO : 0.0
   end
   return p
end

function generate_K(k1, k2, a, b, sync)
   rnd1 = rand()
   rnd2 = sync ? 1 - rnd1 : rand()
   K1 = k1[1] + rnd1*(k1[2]-k1[1])
   K2 = k2[1] + rnd2*(k2[2]-k2[1])
   return [K1 a*K2; a*K1 K2; b*K1 b*K2]
end

function distribute_across(N1, N2, N3, g, p)
   n = [N1, N2, N3]
   choice = hcat(vec(g/2+(1-g).*p), vec(g/2+(1-g).*(1.-p)))
   return choice .* n
end

function WY94(N1, N2, N3; a=0.1, b=0.9, K1=[200,200], K2=[100,100], T=100, g=0.0, r=1.3, linked=true)
   output = zeros((T+1,7))
   K = generate_K(K1, K2, a, b, linked)
   p = habitat_preference([N1, N2, N3], K)
   output[1,2:4] = [N1, N2, N3]
   output[1,5:end] = p
   for i in 1:T
      K = generate_K(K1, K2, a, b, linked)
      p = habitat_preference([N1, N2, N3], K)
      Nij = distribute_across(N1, N2, N3, g, p)
      Wij = exp.(r.*(1.-sum(Nij, 1)./K))
      Nijt = Nij.*Wij
      N1, N2, N3 = vec(sum(Nijt, 2))
      output[(i+1),1] = i
      output[(i+1),2:4] = [N1, N2, N3]
      output[(i+1),5:end] = p
   end
   return output
end

function shannon(n1, n2, n3)
   t = n1 .+ n2 .+ n3
   p1, p2, p3 = n1./t, n2./t, n3./t
   return -(p1 .* log.(p1) .+ p2 .* log.(p2) .+ p3 .* log.(p3))./log(3)
end
