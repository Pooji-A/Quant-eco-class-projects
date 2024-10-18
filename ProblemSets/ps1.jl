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

count_positives([1, -3, 4, 7, -2, 0])
count_positives([-5, -10, 0, 6]) 

