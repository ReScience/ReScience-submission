srand(42)

using Plots
using StatsBase

pyplot()

"""
**Habitat preference**

- `N` is an array with ``N_1, N_2, N_3``
- `K` is a matrix with habitat carrying capacities
- `trials=50` is an optional keyword argument to determine how many iterations should be done

**Returns** an array with ``p_1, p_2, p_3``.
"""
function habitat_preference(N,K;trials=50)
   p = [1.0, 0.0, rand()]
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

"""
**Generate the K matrix**

- `k1` is an array with ``K_{1,min}, K_{1,max}``
- `k2` is an array with ``K_{2,min}, K_{2,max}``
- `a` and `b` are the performances of the specialists and generalists
- `sync` is a Boolean to determine whether the environments are linked or independant

**Returns** a 3 x 2 matrix (species x habitat) of carrying capacity
"""
function generate_K(k1, k2, a, b, sync)
   @assert k1[2] >= k1[1]
   @assert k2[2] >= k2[1]
   @assert a <= 1.0
   @assert b <= 1.0
   @assert a >= 0.0
   @assert b >= a
   rnd1 = rand()
   rnd2 = sync ? 1 - rnd1 : rand()
   K1 = k1[1] + rnd1*(k1[2]-k1[1])
   K2 = k2[1] + rnd2*(k2[2]-k2[1])
   return [K1 a*K2; a*K1 K2; b*K1 b*K2]
end

"""
**Distribute species across habitats**

- `N1`, `N2`, `N3` are the population densities
- `g` is the proportion chosing randomly
- `p` is an array with ``p_1, p_2, p_3``

**Returns** a 3 x 2 matrix (species x habitat) of population densities after habitat selection
"""
function distribute_across(N1, N2, N3, g, p)
   @assert 0.0 <= g
   @assert g <= 1.0
   for i in eachindex(p)
      @assert 0.0 <= p[i] <= 1.0
   end
   n = [N1, N2, N3]
   choice = hcat(vec(g/2+(1-g).*p), vec(g/2+(1-g).*(1.-p)))
   return choice .* n
end


"""
**Performs one run of the model**

- `N1`, `N2`, `N3` are the population densities
- `a=0.1` and `b=0.9` are the performances of the specialists and generalists
- `K1=[200,200]` is an array with ``K_{1,min}, K_{1,max}``
- `K2=[100,100]` is an array with ``K_{2,min}, K_{2,max}``
- `g=0.0` is the proportion chosing randomly
- `r=1.3` is the growth rate
- `T=100` is the number of generations
- `linked=true` is a Boolean to determine whether the environments are linked or independant

**Returns** a matrix with the time, the three population densities, and the three habitat choice values, one generation per line.
"""
function WY94(N1, N2, N3; a=0.1, b=0.9, K1=[200,200], K2=[100,100], T=100, g=0.0, r=1.3, linked=true)
   @assert 0.0 <= g
   @assert g <= 1.0
   @assert a <= 1.0
   @assert b <= 1.0
   @assert a >= 0.0
   @assert b >= a
   @assert 0.0 < r
   @assert K1[2] >= K1[1]
   @assert K2[2] >= K2[1]
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

"""
**Calculates Pielou's evenness for three populations**
"""
function pielou(n1, n2, n3)
   t = n1 .+ n2 .+ n3
   p1, p2, p3 = n1./t, n2./t, n3./t
   return -(p1 .* log.(p1) .+ p2 .* log.(p2) .+ p3 .* log.(p3))./log(3)
end
