# Problem 3: Convergence in the Neoclassical Growth Model

using Plots, Parameters

@with_kw struct NGMProblem
    
    β = 0.95 
    α = 0.3 
    δ = 0.05 
    γ = 2.0 

    f = k -> k^α 
    u = γ == 1 ? c -> log(c) : c -> c ^ (1-γ) / (1-γ)   
    k_star = ((β^(-1) - 1 + δ) / α) ^(1/(α-1)) 

    k_min = 0.75 * k_star 
    k_max = 1.25 * k_star 
    
    n = 100 
    k_grid = range(k_min,stop=k_max,length=n) 
    
end

function T(v, model) 
    @unpack n, k_grid, β, α, δ, f, u = model

    v_new = zeros(n)
    reward = zeros(n, n)
    σ = zeros(n)

    for (k_index, k) in enumerate(k_grid) 
        for (k_next_index, k_next) in enumerate(k_grid) 

            c = k^α - k_next + (1-δ)*k 
            if c > 0
                reward[k_index, k_next_index] = u(c) + β * v[k_next_index]
            else
                reward[k_index, k_next_index] = -Inf
            end
        end 

        v_new[k_index], k_next_index_opt = findmax(reward[k_index, :])
        σ[k_index] = k_grid[k_next_index_opt]
    end
    
    return v_new, σ
end

function vfi(model; maxiter=1000, tol=1e-8) # value function iteration
    @unpack n, k_grid, β, α, δ, f, u = model
    v_init = zeros(n)
    err = tol + 1.0
    iter = 1
    v = v_init
    v_history = [v_init]
    σ = zeros(n)

    while err > tol && iter < maxiter
        v_new, σ = T(v, model)
        err = maximum(abs.(v_new - v)) 
        push!(v_history, v_new)
        v = v_new
        iter += 1
    end

    return v, σ, iter, err, v_history
end

# Task 1: Table showing periods to close half the gap for different γ values
function convergence_table(γ_values, β=0.95, α=0.3, δ=0.05, k0_ratio=0.5)
    results = []
    for γ in γ_values
        model = NGMProblem(β=β, α=α, δ=δ, γ=γ)
        v, σ, _, _, _ = vfi(model)

        k_path = [model.k_grid[1]]
        for t in 2:1000
            k_path = push!(k_path, σ[findfirst(x -> x == k_path[end], model.k_grid)])
            if model.k_star - k_path[end] < 0.5 * (model.k_star - model.k_grid[1])
                push!(results, (γ, t))
                break
            end
        end
    end

    return results
end

# Task 2: Four-panel plot for different γ values
function plot_convergence(γ_values, β=0.95, α=0.3, δ=0.05, k0_ratio=0.5)
    k_star = ((β^(-1) - 1 + δ) / α)^(1 / (α - 1))
    k0 = k0_ratio * k_star
    times = 50

    plt1 = plot(title="Capital Over Time", xlabel="Time (Years)", ylabel="Capital", titlefontsize=10, guidefontsize=8, legend=:bottomright)
    plt2 = plot(title="Output Over Time", xlabel="Time (Years)", ylabel="Output", titlefontsize=10, guidefontsize=8, legend=false)
    plt3 = plot(title="Investment to Output Ratio", xlabel="Time (Years)", ylabel="Investment/Output", titlefontsize=10, guidefontsize=8, legend=false)
    plt4 = plot(title="Consumption to Output Ratio", xlabel="Time (Years)", ylabel="Consumption/Output", titlefontsize=10, guidefontsize=8, legend=false)

    for γ in γ_values
        model = NGMProblem(β=β, α=α, δ=δ, γ=γ)
        v, σ, _, _, _ = vfi(model)

        # Initialize paths
        k_index = argmin(abs.(model.k_grid .- k0))
        k_path = [model.k_grid[k_index]]
        output_path, invest_path, cons_path = [], [], []

        for t in 1:times
            # Update paths based on policy function
            output = model.k_grid[k_index]^α
            investment = σ[k_index] - (1 - δ) * model.k_grid[k_index]
            consumption = output - investment

            push!(output_path, output)
            push!(invest_path, investment / output)
            push!(cons_path, consumption / output)

            # Update capital index for the next period
            k_index = argmin(abs.(model.k_grid .- σ[k_index]))
            push!(k_path, model.k_grid[k_index])
        end

        # Debug investment values
        println("Investment path for γ = $γ: ", invest_path)

        # Add line for this γ to each plot
        plot!(plt1, 1:times, k_path[1:times], label="γ = $γ")
        plot!(plt2, 1:times, output_path)
        plot!(plt3, 1:times, invest_path)
        plot!(plt4, 1:times, cons_path)
    end

    # Combine plots
    plot(plt1, plt2, plt3, plt4, layout=(2, 2))
end



# Example Usage
γ_values = [0.5, 1, 2]
results = convergence_table(γ_values)
println("Convergence Table:", results)
plot_convergence(γ_values)
savefig("convergence_plots.png")

#= 
## Results
Convergence Table: Any[(0.5, 5), (1.0, 8), (2.0, 11)]

## Comments
Considering a low-income country starting with low capital, the time to reach the steady state 
and the behavior during the transition depend on risk aversion (γ) and model parameters:

- The speed of convergence is primarily influenced by risk aversion, as supported by the results. 
  Lower γ prioritizes investment, leading to faster capital accumulation, while higher γ slows this process 
  by favoring consumption.

- The plots align with the expected behavior: during the transition, the economy exhibits high investment 
  and low consumption in the early stages, driving rapid growth in capital and output. 
  As the economy approaches the steady state, investment stabilizes to offset depreciation, 
  and consumption increases.
=#