using LinearAlgebra
using SparseArrays
using ExponentialUtilities
using Plots

function hamiltonian(L::Int, J::Float64, J_imp::Float64)
    center = L ÷ 2
    H = zeros(Complex{Float64}, L, L)
    for i in 1:L-1
        if i == center - 1 || i == center
            bond = J_imp
        else
            bond = J
        end
        H[i, i+1] = bond/2
        H[i+1, i] = bond/2
    end
    return H
end

function time_evolution_op(H, t::Float64)
    evals = eigen(H).values
    evecs = eigen(H).vectors
    U_t = evecs * Diagonal(exp.(-im * evals * t)) * evecs'
    return U_t, evals, evecs
end

function initial_correlation_matrix(H, μ_L, μ_R, L::Int)
    beta=1e6
    _, evals, evecs = time_evolution_op(H, 0.0)
    fL = 1 ./ (1 .+ exp.(beta .* (evals .- μ_L)))
    fR = 1 ./ (1 .+ exp.(beta .* (evals .- μ_R)))
    CL = evecs * Diagonal(fL) * evecs'
    CR = evecs * Diagonal(fR) * evecs'
    C0 = zeros(Complex{Float64}, L, L)
    j0 = L ÷ 2 - 1
    C0[1:j0+1, 1:j0+1] = CL[1:j0+1, 1:j0+1]
    C0[j0+1:end, j0+1:end] = CR[j0+1:end, j0+1:end]
    return C0
end

### Dynamics of fermions

"""
times = 0:0.1:T_steady
onsite_occupations = []
l = 50 ### site of interest
evecs = eigen(H).vectors
evals = eigen(H).values
for t in times
    U_t = evecs * Diagonal(exp.(-im * evals * t)) * evecs'
    C_t = U_t * C0 * U_t'
    push!(onsite_occupations, real(C_t[l, l]))
end
"""

### Dynamics of current
"""
times = 0:0.1:T_steady
current = []
l = 50 ### site of interest
evecs = eigen(H).vectors
evals = eigen(H).values
for t in times
    U_t = evecs * Diagonal(exp.(-im * evals * t)) * evecs'
    C_t = U_t * C0 * U_t'
    push!(current, 2 * J_imp * imag(C_t[l, l+1]))
end
"""

