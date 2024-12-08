# Problem 1: iterative solver for nonlinear equations

using Plots, LinearAlgebra

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
