using LinearAlgebra, Plots

# Define the transition matrix P
P = [0.5 0.3 0.2;
     0.2 0.7 0.1;
     0.3 0.3 0.4]

# Define state spaces for Zt and Xt
Z_states = ["z1", "z2", "z3"]  
X_states = 0:5                 

# policy function 
function policy(Xt, Zt)
    if Zt == "z1"
        return 0
    elseif Zt == "z2"
        return Xt
    elseif Zt == "z3" && Xt <= 4
        return Xt + 1
    elseif Zt == "z3" && Xt == 5
        return 3
    end
end

# joint transition matrix for (Xt, Zt)
function joint_transition_matrix(P, X_states, Z_states, policy)
    num_X = length(X_states)
    num_Z = length(Z_states)
    joint_P = zeros(num_X * num_Z, num_X * num_Z)

    for i in 1:num_X, j in 1:num_Z
        current_state_index = (j - 1) * num_X + i
        for k in 1:num_Z
            new_X = policy(X_states[i], Z_states[j])
            new_X_index = findall(x -> x == new_X, X_states)[1]
            next_state_index = (k - 1) * num_X + new_X_index
            joint_P[current_state_index, next_state_index] = P[j, k]
        end
    end

    return joint_P
end

joint_P = joint_transition_matrix(P, X_states, Z_states, policy)
println("Joint Transition Matrix:")
println(joint_P)

# stationary distribution
function stationary_distribution(joint_P)
    vals, vecs = eigen(joint_P')
    stationary_vec = vecs[:, findall(x -> isapprox(x, 1.0, atol=1e-6), vals)[1]]
    return normalize(stationary_vec, 1)
end

stationary_dist = stationary_distribution(joint_P)
println("Stationary Distribution:")
println(stationary_dist)

# Marginal distribution for Xt
function marginal_X_distribution(stationary_dist, X_states, Z_states)
    num_X = length(X_states)
    num_Z = length(Z_states)
    marginal_X = zeros(num_X)

    for i in 1:num_X
        for j in 1:num_Z
            marginal_X[i] += stationary_dist[(j - 1) * num_X + i]
        end
    end

    return marginal_X
end

marginal_X = marginal_X_distribution(stationary_dist, X_states, Z_states)
println("Marginal Distribution of Xt:")
println(marginal_X)

# Expected value of Xt
function expected_value_X(marginal_X, X_states)
    return sum(marginal_X .* X_states)
end

expected_X = expected_value_X(marginal_X, X_states)
println("Expected Value of Xt:", expected_X)
