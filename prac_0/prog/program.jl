using Base
using Graphs

function approx_ln(a, deg::Int)
    ret = zero(a)

    point = a - one(a)
    n = one(a)

    for i in 1:deg
        val = point^i
        val /= n
        val *= (-1)^(i+1)
        ret += val
    end

    return ret
end

test = approx_ln(3, 4)

println(test)
println(log(3))