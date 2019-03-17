# imports
using DelimitedFiles
using LinearAlgebra
using Printf

# Code
function idx2(i,j)
    if i >= j
        div((i*(i+1)),2) + j
    else
        div((j*(j+1)),2) + i
    end
end

function idx4(i, j, k, l)
    a = i-1
    b = j-1
    c = k-1
    d = l-1
    idx2(idx2(a, b), idx2(c, d)) + 1
end

function main()
    # load data
    convergence_DM = readdlm("../../data/convergence_DM.txt")[1]
    convergence_E = readdlm("../../data/convergence_E.txt")[1]
    S = readdlm("../../data/S.txt")
    T = readdlm("../../data/T.txt")
    V = readdlm("../../data/V.txt")
    eri = readdlm("../../data/eri.txt")
    E_nuc = readdlm("../../data/E_nuc.txt")[1]
    iteration_max = readdlm("../../data/iteration_max.txt", Int)[1]
    num_ao = readdlm("../../data/num_ao.txt", Int)[1]
    num_elec_alpha = readdlm("../../data/num_elec_alpha.txt", Int)[1]
    num_elec_beta = readdlm("../../data/num_elec_beta.txt", Int)[1]
    iteration_max = readdlm("../../data/iteration_max.txt", Int)[1]

    # loop variables
    iteration_num = 0
    E_total = 0
    E_elec = 0.0
    iteration_E_diff = 0.0
    iteration_rmsc_dm = 0.0
    converged = false
    exceeded_iterations = false

    D = zeros(num_ao, num_ao)
    
    s, L = eigen(Symmetric(S))
    X = zeros(num_ao, num_ao)
    s = [1.0/sqrt(i) for i in s]
    X = Diagonal(s)
    X = L * X * transpose(L)

    H = T + V
    while (!converged & !exceeded_iterations)
        # store last iteration and increment counters
        iteration_num += 1
        E_elec_last = E_elec
        D_last = deepcopy(D)
        # form G matrix
        G = zeros(num_ao, num_ao)
        for i in 1:num_ao
            for j in 1:num_ao
                for k in 1:num_ao
                    for l in 1:num_ao
                        G[i, j] += D[k, l] * ((2.0*(eri[idx4(i, j, k, l)])) - (eri[idx4(i, k, j, l)]))
                    end
                end
            end
        end
        # build fock matrix
        F = H + G
        F_prime = X * F * X
        # solve the eigenvalue problem
        E_orbitals, C_prime = eigen(Symmetric(F_prime))
        C = X * C_prime
        # compute new density matrix
        D = zeros(num_ao, num_ao)
        for i in 1:num_ao
            for j in 1:num_ao
                for k in 1:num_elec_alpha
                    D[i, j] += C[i, k] * C[j, k]
                end
            end
        end
        # calculate electronic energy
        E_elec = sum(D .* (H + F))
        # calculate energy change of iteration
        iteration_E_diff = abs(E_elec - E_elec_last)
        # rms change of density matrix
        iteration_rmsc_dm = sqrt(sum((D - D_last)^2))
        if((abs(iteration_E_diff) < convergence_E) & (iteration_rmsc_dm < convergence_DM))
            converged = true
        end
        if(iteration_num == iteration_max)
            exceeded_iterations = true
        end
    end
    E_total = E_elec + E_nuc
end

E_total = main()
@printf("%20.15f\n",E_total)
