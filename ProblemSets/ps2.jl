# Problem 1: iterative solver for nonlinear equations

using Plots, LinearAlgebra , Roots

function iterative_solver(f, x0, α; ε=1e-6, maxiter=1000)
    # store iterations
    x_history = [x0] # points
    residuals = Float64[]
    
    x_current = x0
    
    for i in 1:maxiter
        # Calculate g(x) = f(x) + x
        g_x = f(x_current) + x_current
        
        # Calculate next iteration with dampening
        x_next = (1 - α) * g_x + α * x_current
        
        # Calculate residual
        residual = abs(x_next - x_current)
        push!(residuals, residual)
        
        push!(x_history, x_next)
        
        # Check convergence condition
        if residual < ε *(1+abs(x_current)) #
            return (
                flag = 0,  
                solution = x_next,
                value = f(x_next),
                error = abs(x_next - g_x),
                history = x_history,
                residuals = residuals
            )
        end
        
        x_current = x_next
    end
    
    # Reached maxiter we return NaN
    return (
        flag = 1,
        solution = NaN,
        value = NaN,
        error = NaN,
        history = x_history,
        residuals = residuals
    )
end

# Test functions
f1(x) = (x + 1)^(1/3) - x
h1(x) =  x^3 -x -1    


# Test cases

println("Testing with f(x) = (x + 1)^(1/3) - x and x0 = 1.0")
for α in range(0.0, stop=0.9, step=0.1)
    println("\nα = $α:")
    result = iterative_solver(f1, 1.0, α)
    println("Solution: ")
    if result.flag == 0
        println("Converged in $(length(result.history)) iterations")
        println("Solution: $(result.solution)")
        println("Error: $(result.error)")
    else
        println("Failed to converge")
    end
end

println("Testing with h(x) = x^3 -x -1 and x0 = 1.0")
for α in range(0.0, stop=0.9, step=0.1)
    println("\nα = $α:")
    result = iterative_solver(h1, 100.0, α)
    println("Solution: ")
    if result.flag == 0
        println("Converged in $(length(result.history)) iterations")
        println("Solution: $(result.solution)")
        println("Error: $(result.error)")
    else
        println("Failed to converge")
    end
end


result = iterative_solver(f1, 1.0, 0)

result = iterative_solver(h1, 100.5, 0.5)




# Problem 2: Linear algebra
# Ax = b

using LinearAlgebra, Printf, PrettyTables


function exact_solution(α, β)
    # Get exact solution Ax=b by solving system backwards

    # from the equations:
    x5 = 1.0
    x4 = x5
    x3 = x4
    x2 = x3
    # From first equation: x1 + (-1)x2 + (α-β)x4 + βx5 = α
    x1 = α + x2 - (α-β)*x4 - β*x5
    
    return [x1, x2, x3, x4, x5]
end




function solve_system(α, β)
    # Returns exact solution, numerical solution, condition number, and relative residual.
    
    A = Float64[
        1 -1  0  α-β  β;
        0  1 -1  0    0;
        0  0  1 -1    0;
        0  0  0  1   -1;
        0  0  0  0    1
    ]
    
    b = Float64[α, 0, 0, 0, 1]
    
    # Get exact and numerical solution
    x_exact = exact_solution(α, β)
    x_numerical = A \ b
    
    # Calculate relative residuals
    rel_residual = norm(A*x_numerical - b) / norm(b)
    
    cond_num = cond(A)
    
    return (
        exact = x_exact,
        numerical = x_numerical,
        residual = rel_residual,
        cond = cond_num
    )
end



# Create analysis table

b  = 0.1
β_values = Float64[10.0^i for i in 0:12]

x1_exact = Float64[]
x1_numerical = Float64[]
cond_numbers = Float64[]
residuals = Float64[]

for β in β_values
    result = analyze_system(α, β)
    push!(x1_exact, result.exact[1])
    push!(x1_numerical, result.numerical[1])
    push!(cond_numbers, result.cond)
    push!(residuals, result.residual)
end


header = ["β", "x1 (exact)", "x1 (numerical)", "cond(A)", "relative residual"]
data = hcat(β_values, x1_exact, x1_numerical, cond_numbers, residuals)

pretty_table(data;
    header = header,
    formatters = (ft_printf("%.0e", [1]), 
                    ft_printf("%.6f", [2,3]),
                    ft_printf("%.2e", [4,5])))

                    
# Problem 3 : Internal Rate of Return  (IRR)
using Roots
# NPV function
function NPV(r,c)
    return sum(c[t]/ (1+r)^t for t in 1:length(c))
end

# IRR solver
function internal_rate(c)
    f(r) = NPV(r,c)
    root = find_zero(f,(0.0, 1.0))
    return root 
end

# Example usage of the functions created 
c = [-5.0, 0.0, 2.5, 5.0]
irr = internal_rate(c)
println("The internal rate of return (irr) is : ",irr)

# Problem 4: CES Production Function

using Plots
using Optim

# CES production function
function production_function(x1, x2, α, σ)
    return (α * x1^(σ - 1) / σ + (1 - α) * x2^(σ - 1) / σ)^(σ / (σ - 1))
end

# Contour plot for the production function
function plot_production_function(α, σ)
    x1_vals = LinRange(0.1, 10, 100)
    x2_vals = LinRange(0.1, 10, 100)
    Z = [production_function(x1, x2, α, σ) for x1 in x1_vals, x2 in x2_vals]
    contour(x1_vals, x2_vals, Z, title="CES Production Function", xlabel="x1", ylabel="x2")
    
    # Save the plot for production function
    savefig("CES_production_function.png")  
end
# Example usage for production function plot
plot_production_function(0.5, 0.25)

# Cost function and optimization
function cost_function(w1, w2, y, α, σ)
    # Define the objective function (cost function)
    objective(x) = w1 * x[1] + w2 * x[2]  
    
    # Optimization using the COBYLA method (no need for gradient)
    result = optimize(objective, [1.0, 1.0], LowerBound(0.1), UpperBound(10.0), method = :cobyla)
    return result
end

# Plot demand and cost for different values of σ
function plot_demand_and_cost(w1, w2, y, α, σ_vals)
    plot_data = []
    for σ in σ_vals
        result = cost_function(w1, w2, y, α, σ)
        println("Result for σ = ", σ, ": ", result)
        push!(plot_data, result)
    end
    
    # Plotting cost and demand for different σ values
    σ_vals = [0.25, 1.0, 4.0]
    demands = [result.minimizer[1] for result in plot_data]  
    costs = [result.minimum for result in plot_data]  
    plot(σ_vals, demands, label="Demand for Input 1", xlabel="σ", ylabel="Demand (x1)", title="Cost and Input Demand vs σ", linewidth=2)
    plot!(σ_vals, costs, label="Cost", ylabel="Cost", linewidth=2)
    
    # Save the plot for demand and cost
    savefig("cost_and_demand.png")  
end

# Example: plot cost and input demand for different values of σ
plot_demand_and_cost(1.0, 1.0, 1.0, 0.5, [0.25, 1, 4])