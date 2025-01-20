using Statistics

struct OrchidParams
    X::Float64        # Maximum willingness to pay
    c::Float64        # Base cost
    f::Float64        # Mental cost per vendor visit 
    q::Float64        # Probability of finding orchid
    pmin::Float64     # Minimum price
    pmax::Float64     # Maximum price
    N::Int64         # Number of vendors
end

function solve_bellman_equation(params::OrchidParams)
    
    # prices with 0.1 increments
    prices = params.pmin:0.1:params.pmax
    price_prob = ones(length(prices)) / length(prices)
    
    # Value functions
    v = zeros(params.N + 1)      # Main value function
    vT = zeros(params.N + 1)     # Value of terminating
    vB = zeros(params.N + 1, length(prices))  # Value of buying
    vA = zeros(params.N + 1)     # Value of approaching
    
    # Policy functions (b)
    σ_approach = zeros(Bool, params.N + 1)
    σ_buy = zeros(Bool, params.N + 1, length(prices))
    
    # Backward induction
    for n in params.N:-1:0
        
        # Terminal values
        vT[n+1] = -params.c * n - params.f * n
        
        # Value of buying at each price
        for (p_idx, p) in enumerate(prices)
            vB[n+1, p_idx] = params.X - p - params.c * n - params.f * n
        end
        
        # Value of approaching another vendor
        if n < params.N
            EV_found = sum(max.(vB[n+2,:], vT[n+2]) .* price_prob)
            v_not_found = v[n+2]
            vA[n+1] = -params.f + params.q * EV_found + (1-params.q) * v_not_found
        end
        
        # Optimal value and policies
        v[n+1] = max(vA[n+1], vT[n+1]) # a) value function
        σ_approach[n+1] = vA[n+1] > vT[n+1]
        
        # Buy policy
        for (p_idx, p) in enumerate(prices)
            continue_value = max(vA[n+1], vT[n+1])
            buy_value = params.X - p - params.c * n - params.f * n
            σ_buy[n+1, p_idx] = buy_value >= continue_value
        end
    end
    
    # c) Function to calculate buy probability and expected price
    function calc_probability_and_price(n::Int64)
        if n > params.N
            return 0.0, 0.0
        end
        
        # Probability of reaching this state
        p_reach = (1-params.q)^(n-1) * params.q
        
        # Probability of buying at each price
        p_buy_price = p_reach * price_prob .* σ_buy[n,:]
        
        # Total probability of buying at this n
        total_prob = sum(p_buy_price)
        
        # Expected price conditional on buying
        exp_price = total_prob > 0 ? sum(p_buy_price .* prices) / total_prob : 0.0
        
        return total_prob, exp_price
    end
    
    # Analyze price acceptance behavior over n
    price_acceptance_thresholds = fill(Inf, params.N + 1)
    for n in 1:params.N+1
        # For each n, find the maximum acceptable price
        max_acceptable = -Inf
        for (p_idx, p) in enumerate(prices)
            if σ_buy[n, p_idx]
                max_acceptable = max(max_acceptable, p)
            end
        end
        if max_acceptable > -Inf
            price_acceptance_thresholds[n] = max_acceptable
        end
    end
    
    return v, σ_approach, σ_buy, calc_probability_and_price, price_acceptance_thresholds
end


# Example usage
params = OrchidParams(
    50.0,   # Maximum willingness to pay
    0.5,    # base cost
    0.5,   # f is undefined in the problem, so I set it to 0.5
    0.15,   # q
    10.0,   # pmin
    100.0,  # pmax
    50      # N
)

# Get results
v, σ_approach, σ_buy, prob_calc, price_thresholds = solve_bellman_equation(params)

# Calculate total probability of buying and expected price
total_prob = 0.0
weighted_price_sum = 0.0
weighted_n_sum = 0.0

for n in 1:params.N
    prob, exp_price = prob_calc(n)
    total_prob += prob
    weighted_price_sum += prob * exp_price
    weighted_n_sum += n * prob
end

# Calculate final metrics
prob_buy = total_prob
exp_price = weighted_price_sum / total_prob
exp_vendors = weighted_n_sum / total_prob

# Analyze results
println("a) Probability of buying: ", round(prob_buy, digits=3))
println("b) Expected price if bought: ", round(exp_price, digits=2))
println("c) Expected number of vendors approached: ", round(exp_vendors, digits=2))
println("d) As Basil approaches more vendors, his costs increase and he would be less willing to agree for a higher price. This is because he has to pay the mental cost of approaching more vendors, and the willingness to pay is defined by the maximum willingness to pay minus the costs and the price of the orquid. This can be seen in the price acceptance thresholds over time.")

println("\nPrice acceptance thresholds over time:")
println("Vendor | Max Acceptable Price")
println("------------------------")
for n in 1:params.N+1
    if price_thresholds[n] < Inf
        println("$(n-1)\t| $(round(price_thresholds[n], digits=1))")
    end
end

println("\nDetailed probabilities and prices by vendor:")
for n in 1:params.N
    prob, exp_price = prob_calc(n)
    if prob > 0
        println("n = $n: Probability = $(round(prob, digits=3)), Expected price = $(round(exp_price, digits=2))")
    end
end