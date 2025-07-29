using LinearAlgebra
using SparseArrays
using ExponentialUtilities
using Plots
using LaTeXStrings

# Create a free-fermion Hamiltonian with a weak-link at sites L/2-1, L/2 and L/2+1

function Hamiltonian(L::Int, Jimp::Float64)
    H = zeros(Complex{Float64}, L, L)
    for i in 1:(L ÷ 2 - 1)
        H[i, i + 1] = 0.5
        H[i + 1, i] = 0.5
    end
    for i in (L ÷ 2 + 1):L
        H[i, i - 1] = 0.5
        H[i - 1, i] = 0.5
    end
    H[L ÷ 2 - 1, L ÷ 2] = Jimp ÷ 2
    H[L ÷ 2, L ÷ 2 - 1] = Jimp ÷ 2
    H[L ÷ 2, L ÷ 2 + 1] = Jimp ÷ 2
    H[L ÷ 2 + 1, L ÷ 2] = Jimp ÷ 2
    return H
end

## Time evolution operator
"""
function time_evolution_op(H, t::Float64)
    U = exp(-im * t * H)
    return U
end
"""

function time_evolution_op(H, t::Float64)
    eig = eigen(H)
    V = eig.vectors
    D = Diagonal(exp.(-im * t * eig.values))
    U = V * D * V'
end

function initial_state(L::Int, defect_loc::Int)
    initial_state = zeros(Complex{Float64}, L)
    for i in 1:(defect_loc - 1)
        initial_state[i] = 1.0
    end
    initial_state = initial_state / norm(initial_state)
    return initial_state
end

function evolve_state(L::Int, Jimp::Float64, defect_loc::Int, T::Float64)

    # Build Hamiltonian and initial state
    H = Hamiltonian(L, Jimp)
    state = initial_state(L, defect_loc)

    # Time evolution setup
    psi = [];
    push!(psi, state)

    # Time evolution operator  

    for t in 0.1:0.1:T
        U = time_evolution_op(H, t)
        state = U * state
        state /= norm(state)
        push!(psi, state)
    end

    return psi
end

function defect_bond_op(L::Int, J_imp::Float64)
    op = zeros(Complex{Float64}, L, L)

    for i in Int(L/2-1):Int(L/2+1)
        op[i, i+1] = 1.0 * im * J_imp
        op[i+1, i] = -1.0 * im * J_imp
    end
    return op
end