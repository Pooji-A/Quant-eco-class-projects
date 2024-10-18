using Plots
using Statistics
using DelimitedFiles


function odd_or_even(n)
    if iseven(n)
        println("Even")
    else
        println("Odd")
    end
end


odd_or_even(6)


function compare_three(a, b, c)
    ispositive = x -> x > 0
    iszero = x -> x == 0
    if all(map(ispositive, [a, b, c]))
        println("All numbers are positive")
    elseif all(map(iszero, [a, b, c]))
        println("All numbers are zero")
    else
        println("At least one number is not positive")
    end

end


compare_three(1,2,3)
compare_three(0,0,0)
compare_three(0,0,-4)


function my_factorial(n)
    result = 1
    for i in 1:n
        result *= i
    end
    return result
end

my_factorial(5)
my_factorial(7)


function count_positives(arr)
    counter = 0
    for num in arr
        if num > 0
            counter += 1
        end
    end
    return counter
end

println(count_positives([1, -3, 4, 7, -2, 0]))
println(count_positives([-5, -10, 0, 6]))



function plot_powers(n)
    plot(legend=true)
    for i in 1:n
        pow(x) = x^i
        plot!(pow,-10:0.1:10)
    end
    p=plot!(legend=true)
    savefig(p,"plot_powers_output.png")
end    

plot_powers(2)


function standard_deviation(x)
    mu_hat = sum(x)/length(x)
    d = x.- mu_hat
    squared_d = d.^ 2
    variance = sum(squared_d)/(length(x) - 1)
    return sqrt(variance)
end

println(standard_deviation([1,2,3,4,5]))
println(standard_deviation([5,10,15]))
println(standard_deviation(collect(2:7)))


data = readdlm("dataset.txt", '\t', Float64)

#println(data)

earnings = data[:,1]
education = data[:,2]
hours_worked = data[:,3]

p1=scatter(education, earnings, xlabel="Education", ylabel="Earnings", title="Relationship between education and earnings",color = "green")
p2=scatter(hours_worked, earnings, xlabel="Hours Worked", ylabel="Earnings", title="Relationship between hours worked and earnings",color = "red")

plot(p1)
savefig(p1,"education_vs_earnings.png")

plot(p2)
savefig(p2,"hoursworked_vs_earnings.png")

corr_education_earnings = cor(education, earnings)
corr_hoursworked_earnings = cor(hours_worked, earnings)

println("Correlation between earnings and education:",corr_education_earnings)
println("Correlation between hours worked and earnings:",corr_hoursworked_earnings)

