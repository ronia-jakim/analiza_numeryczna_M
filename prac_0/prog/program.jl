using Base
#using Graphs

function approx_ln(x::Float64, deg::Int)
    point = Float64(2.72)
    # x - one(x)
    fac = 1.0
    k = zero(x)

    ret = zero(x)

    for i in 0:deg
        val = zero(x)
        if k != zero(x)
            val = (x-point)^k
            val /= k * (point^k)
        else
            val = 1
            val /= point^k
        end
        k += one(x)
        val *= (-1)^k
        ret += val
    end

    return ret
end

t = approx_ln(Float64(3), 4)
println(t, " ", log(3))

t = approx_ln(Float64(3), 8)
println(t, " ", log(3))

t = approx_ln(Float64(3), 16)
println(t, " ", log(3))