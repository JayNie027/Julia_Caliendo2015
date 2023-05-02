
function compute_wage(w::Array{Float64,1},
    P::Array{Float64,2},
    A::Array{Float64,2},
    tau::Array{Float64,3},
    theta::Array{Float64,1},
    gamma::Array{Float64,3},
    N::Int64,
    J::Int64,
    beta::Array{Float64,2},
    alpha::Array{Float64,2})

    P = Array{Float64,2}(undef, N, J)
    bts = Array{Float64,3}(undef, N, N, J)
    u = Array{Float64,2}(undef, N, J)
    EX = Array{Float64,1}(undef, N)
    IM = Array{Float64,1}(undef, N)
    w = copy(w0)
    Tw = copy(w0)

    X = Array{Float64,2}(undef, N, J)
    I_temp = Array{Float64,2}(undef, N, J)


    tol = 1.0e-4 # Tolerance for convergence
    dif = 10.0 # Initialize difference from trade balance
    iter = 0
    while (dif > tol) & (iter < 10000)

        iter += 1

        w = copy(Tw)


        P, bts, u = compute_price_tradeshares(w::Array{Float64,1},
                                              P::Array{Float64,2},
                                              A::Array{Float64,2},
                                              tau::Array{Float64,3},
                                              theta::Array{Float64,1},
                                              gamma::Array{Float64,3},
                                              N::Int64,
                                              J::Int64,
                                                beta::Array{Float64,2},
                                                alpha::Array{Float64,2})

        bts_t = permutedims(bts, (2, 1, 3))
        bts_temp = repeat(bts_t, inner=(J, 1, 1))
        bts_temp1 = permutedims(bts_temp, (1, 3, 2))
        bts_m = reshape(bts_temp1, N * J, N * J)

        # make 3D 2*4*2 matrix, than stack third dimension
        beta_m = repeat(beta, 1,1, J)
        beta_m1 = permutedims(beta_m, (2,3,1))
        temp = repeat(gamma .* beta_m1 , 1, 2, 1)
        temp1 = permutedims(temp, (1, 3, 2))
        gamma_m = reshape(temp1, N * J, N * J)

        # I 
        I_temp = bts_m .* gamma_m
        temp_I = Matrix(I, N * J, N * J) .- I_temp
        temp_I1 = inv(temp_I)
        # for n in 1:N*J
        #     for i in 1:N*J
        #         if temp_I1[n, i] < 0.0001
        #             temp_I1[n, i] = 0.0001
        #         end
        #     end
        # end


        VA = repeat(w, inner=J) .* repeat(L, inner=J)
        X = temp_I1 * (reshape(transpose(alpha), N * J, 1) .* VA)
        X = Matrix(transpose(reshape(X, J, N)))

        EX = Array{Float64,1}(undef, N)
        IM = Array{Float64,1}(undef, N)
        for n in 1:N
            EX[n] = 0.0
            IM[n] = 0.0
            for i in 1:N
                for j in 1:J
                    if n != i
                        EX[n] += X[i, j] .* bts[i, n, j] # country n's export
                        IM[n] += X[n, j] .* bts[n, i, j]
                    end
                end
            end
        end

        # define excess demand (net export to GDP)
        Zw = (EX .- IM) ./ (w .* L)
        dif = maximum(abs.(Zw))

        # update the wages and guess
        scl = 0.1
        Tw = w .* (1.0 .+ scl .* Zw)
        #   if mod(iter,100)==0
        @printf " %d iterations, Excess demand: %f\n" iter dif

        #   end

    end

    return w, P, bts, X, EX, IM

end