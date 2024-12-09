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


using Plots, NLopt

# CES production function
function production_function(x1, x2, α, σ)
    if σ == 1
        return x1^α * x2^(1-α)  # Cobb-Douglas case
    else
        return (α * x1^((σ-1)/σ) + (1-α) * x2^((σ-1)/σ))^(σ/(σ-1))
    end
end

# Part 1: Function to create contour plot
function plot_production_function(α, σ, x1_range=0.1:0.1:10, x2_range=0.1:0.1:10)
    Z = [production_function(x1, x2, α, σ) for x1 in x1_range, x2 in x2_range]
    contour(x1_range, x2_range, Z, 
           title="σ = $σ",
           xlabel="x1",
           ylabel="x2",
           fill=true)
end

# Part 2: Create plot using the previous function
α = 0.5
x1_range = 0.1:0.1:100
x2_range = 0.1:0.1:100

p1 = plot_production_function(α, 0.25, x1_range, x2_range);
p2 = plot_production_function(α, 1.0, x1_range, x2_range);
p3 = plot_production_function(α, 4.0, x1_range, x2_range);

plot(p1, p2, p3, layout=(1,3), size=(1000, 350), 
margin=5Plots.mm,
plot_title="CES Production Functions",
plot_titlevspan=0.1)


# Part 3: Cost minimization function
function minimize_cost(α, σ, w1, w2, y)
    function objective(x, gradx1)
        if length(grad) > 0
            grad[1] = w1
            grad[2] = w2
        end
        return w1*x[1] + w2*x[2]
    end
    
    function constraint(xx1, gradx1)
        if length(grad) > 0
            if σ == 1
                grad[1] = α * x[1]^(α-1) * x[2]^(1-α)
                grad[2] = (1-α) * x[1]^α * x[2]^(-α)
            else
                denom = x[1]^((σ-1)/σ) + (1-α) * x[2]^((σ-1)/σ)
                grad[1] = α * ((σ-1)/σ) * x[1]^((σ-1)/σ-1)
                grad[2] = (1-α) * ((σ-1)/σ) * x[2]^((σ-1)/σ-1)
            end
        end
        return production_function(x[1], x[2], α, σ) - y
    end

    opt = Opt(:LD_SLSQP, 2)
    opt.lower_bounds = [1e-6, 1e-6]
    opt.upper_bounds = [100.0, 100.0]
    opt.ftol_rel = 1e-6
    opt.maxeval = 1000
    
    opt.min_objective = objective
    equality_constraint!(opt, constraint)
    
    (minf, minx, ret) = optimize(opt, [1.0, 1.0])
    return minf, minx
end

# Part 4: Plot cost function and input demands
α = 0.5
w2 = 1.0
y = 1.0
w1_range = 0.1:0.1:5

σ_values = [0.25, 1.0, 4.0]
colors = [:blue, :red, :green]

costs = zeros(length(w1_range), 3)
x1_demands = zeros(length(w1_range), 3)
x2_demands = zeros(length(w1_range), 3)


for (j, σ) in enumerate(σ_values)
    for (i, w1) in enumerate(w1_range)
        try
            cost, x = minimize_cost(α, σ, w1, w2, y)
            costs[i,j] = cost
            x1_demands[i,j] = x[1]
            x2_demands[i,j] = x[2]
        catch e
            costs[i,j] = NaN
            x1_demands[i,j] = NaN
            x2_demands[i,j] = NaN
        end
    end
end

# Create  plot
p1 = plot(title="Cost Function");
p2 = plot(title="x1 Demand");
p3 = plot(title="x2 Demand");

for (j, σ) in enumerate(σ_values)
    plot!(p1, w1_range, costs[:,j], label="σ=$σ", color=colors[j])
    plot!(p2, w1_range, x1_demands[:,j], label="σ=$σ", color=colors[j])
    plot!(p3, w1_range, x2_demands[:,j], label="σ=$σ", color=colors[j])
end

for p in [p1, p2, p3]
    xlabel!(p, "w1")
end

plot(p1, p2, p3, layout=(1,3), size=(1200,400))



