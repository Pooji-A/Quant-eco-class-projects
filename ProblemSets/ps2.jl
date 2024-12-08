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




