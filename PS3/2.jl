using Distributions, Plots

function create_job_search_model(;
    n=100, # wage grid size
    w_min=10.0, # lowest wage
    w_max=60.0, # highest wage
    a=200, # wage distribution parameter
    b=100, # wage distribution parameter
    β=0.96, # discount factor
    c=10.0, # unemployment compensation
    p=0.0  # probability of losing job
    )
    w_vals = collect(LinRange(w_min, w_max, n))
    ϕ = pdf(BetaBinomial(n-1, a, b))
    return (; n, w_vals, ϕ, β, c, p)
end

function solve_VE(w, β, p, expected_VU)
    return (w + β*p*expected_VU)/(1 - β*(1-p))
end

# Bellman operator 
function T_operator!(VU_new, VU, model)
    (; n, w_vals, ϕ, β, c, p) = model
    expected_VU = sum(VU .* ϕ)
    
    for i in 1:n
        VE_w = solve_VE(w_vals[i], β, p, expected_VU)
        VU_new[i] = max(VE_w, c + β*expected_VU)
    end
end

function solve_separation_model(model; tol=1e-8, maxiter=1000)
    (; n, w_vals, ϕ, β, c, p) = model
    
    VU = zeros(n)
    VU_new = similar(VU)
    
    error = tol + 1.0
    iter = 0
    
    while error > tol && iter < maxiter
        copyto!(VU_new, VU)
        T_operator!(VU_new, VU, model)
        error = maximum(abs.(VU_new - VU))
        copyto!(VU, VU_new)
        iter += 1
    end
    
    expected_VU = sum(VU .* ϕ)
    VE = [solve_VE(w, β, p, expected_VU) for w in w_vals]
    reservation_wage = w_vals[findfirst(VE .≥ VU)]
    accept_prob = sum(ϕ[findall(w_vals .≥ reservation_wage)])
    exp_duration = 1/accept_prob
    
    return reservation_wage, accept_prob, exp_duration, VE, VU
end

# Test the special case p = 0
model_class = create_job_search_model(p=0.0)
w_star_class, q_class, dur_class, VE_class, VU_class = solve_separation_model(model_class)

w_test = model_class.w_vals[50]
VE_analytical = w_test/(1-model_class.β)
VE_computed = VE_class[50]
println("Test p=0 case:")
println("VE analytical: ", VE_analytical)
println("VE computed: ", VE_computed)

p_grid = 0.0:0.02:0.5
results = []

for p in p_grid
    model = create_job_search_model(p=p)
    res = solve_separation_model(model)
    push!(results, res)
end

w_stars = [r[1] for r in results]
accept_probs = [r[2] for r in results]
durations = [r[3] for r in results]

p1 = plot(p_grid, w_stars,
    label="Reservation Wage",
    xlabel="Probability of Losing Job (p)",
    ylabel="w*",
    title="Reservation Wage vs Probability of Losing Job")

p2 = plot(p_grid, accept_probs,
    label="Acceptance Probability",
    xlabel="Probability of Losing Job (p)",
    ylabel="q",
    title="Job Acceptance Probability vs Probability of Losing Job")

p3 = plot(p_grid, durations,
    label="Expected Duration",
    xlabel="Probability of Losing Job (p)",
    ylabel="Periods",
    title="Expected Unemployment Duration vs Probability of Losing Job")

plot(p1, p2, p3, layout=(3,1), size=(800,900))