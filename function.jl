function compute_price_tradeshares(w::Array{Float64,1},
                                   P::Array{Float64,2},
                                   A::Array{Float64,2},
                                    tau::Array{Float64,3},
                                    theta::Array{Float64,1},
                                    gamma::Array{Float64,3},
                                    N::Int64,
                                    J::Int64,
                                    beta::Array{Float64,2},
                                    alpha::Array{Float64,2})


    u = Array{Float64,2}(undef,N,J)
    bts = Array{Float64, 3}(undef, N,N,J)
    P = copy(P0)
    Tp = copy(P0)

    tol = 1.0e-4 # Tolerance for convergence
    dif = 10.0 # Initialize difference from prices
    iter = 0
    while (dif>tol) & (iter<10000)
        
    iter += 1

    P = copy(Tp)

        for n in 1:N
            for j in 1:J
                comp = ones(N,J)
                    for k in 1:J 

                    comp[n,j] *= P[n,k] .^ gamma[k,j,n]

                    end

                u[n,j] = Bnj[n,j] .* ((w[n]).^(beta[n,j]) .* (comp[n,j]).^(1.0 .-beta[n,j]))

            end
        end
# print(u)

        for n in 1:N
            for j in 1:J 
                denom = zeros(N, J)

                for i in 1:N

                    denom[n,j] += A[i,j] .* ( u[i,j] .* tau[n,i,j]).^(-theta[j])

                end

                Tp[n,j] = denom[n,j] .^(-1.0./theta[j])

                for i in 1:N 

                    bts[n,i,j] = tau[n, i, j].^(-theta[j]) .* A[i,j] .* u[i,j].^(-theta[j])./ denom[n,j]
                end

            end
        end

    pmax = Tp .- P
    dif = maximum(abs.(pmax)) 

    # @printf " %d iterations. price diff: %f \n" iter dif

    end 

        return P, bts, u

end

#---------------------------------------------------------------------------------------------------