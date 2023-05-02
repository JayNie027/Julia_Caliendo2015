
using Printf
using LinearAlgebra
# This file is replicate the multi country and multi sector model from Caliendo 'welfare analysis for NAFTA'
# We assume two countries and two sectors with sectoral linkages
# to simplify, assume trade is balance in each sector. Since usually we assume exogenous trade deficit, one can add that in trade equlibrium condition

include("function.jl")
include("solve_eq.jl")

N = 2;
J = 2;

# assume productivity 
A = ones(N,J);
# each row is country, and each column is sector 

# labor 
L = 10 .* ones(N)
# L1 = sum(L, dims=2)
# L[end,end] = 2 # country 2 has 2 labor units at sector 2

# sectoral linkage: assume identical for countries 
gamma = 0.5 .* ones(J,J)     # each row is supply sector, and each column is demand sector, so column sum is 1
# each row is k, each col is j, gamma 12 means sector 2 use intermediates from sector 1
# this construction is consistent with input ouput table， kj (supply, demand)
gamma = repeat(gamma, 1, 1, N)

# consumption share nj
alpha = 0.5 .* ones(N,J) # consumption share on each sector, each row is n, each column is sector j

# trade cost 
tau = ones(N, N, J)
for n in 1:N
    for i in 1:N
        for j in 1:J
            if n != i
                tau[n, i, j] = 1.0
            end
        end
    end
end    

# value added share, nj 
beta = 0.5 .* ones(N,J)    

Bnj = ones(N, J)
for n in 1:N
    for k in 1:J
        for j in 1:J
            Bnj[n, j] *= (beta[n, j]) .^ (-1 .* beta[n, j]) .* ((1 .- beta[n, j]) .* gamma[k, j, n])^(-1 .* (1 .- beta[n, j]) .* gamma[k, j, n])
        end
    end
end
Bnj # constant in marginal cost function, not that important

# define trade elasticity, assume both sector equals to 4 
theta = 4.0 .* ones(N)


# defined all exogenous variables
#------------------------------------------------------------------------------------------------------------------ 


# w0 = ones(N) 
w0 = [1.0, 7.0] # initial guess for wage 
# free mobie labor so equal wage across sectors
P0 = 0.5 .* ones(N, J)      # initial guess for price
w, P, bts, X, EX, IM, = compute_wage(w0::Array{Float64,1},
                                            P0::Array{Float64,2},
                                            A::Array{Float64,2},
                                            tau::Array{Float64,3},
                                            theta::Array{Float64,1},
                                            gamma::Array{Float64,3},
                                            N::Int64,
                                            J::Int64,
                                            beta::Array{Float64,2},
                                            alpha::Array{Float64,2})

w
PC = ones(N)
for n in 1:N
    for j in 1:J 
        PC[n] *=(P[n,j]./ alpha[n,j])^(alpha[n,j])
    end
end

w./PC
