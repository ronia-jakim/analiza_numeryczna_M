using Base
setprecision(BigFloat, 1024)
function fact(x::Int)
    res = 1.0
    for i in 1::x
        res *= i
    end
end

module chudnowsky
# https://en.wikipedia.org/wiki/Chudnovsky_algorithm



    function calc(iterations::Int)
        res = BigFloat(0)
        for q in 0:iterations
            sign = 1 
            if q % 2 
                sign = -1
            end
            a = fact(6 * q) / (fact(3 * q) * fact(q))
            b = (545140134 * q + 13591409) / ((640320)^(3 * q + 3/2))
            res += sign * a * b
        end

        return 1.0/(res * 12.0)
    end

end

module ramanujan
# https://en.wikipedia.org/wiki/Srinivasa_Ramanujan#Mathematical_achievements

    function calc(iterations::Int)
        res = BigFloat(0)
        for q in 0:iterations
            a = fact(4 * q) * (1103 + 26390 * k)
            b = fact(q) ^ 4 * 396 ^ (4 * k)
            rest += a / b
        end
        return 1.0 / ((2 * sqrt(2) * res) / 9801)
    end

end

module montecarlo
    function experiment()
        x = rand()
        y = rand()
        d = x * x + y * y
        if d <= 1
            return true
        else 
            return false
        end
    end

    function calc(iterations::Int)
        inside = BigFloat(0)
        one = BigFloat(1)
        for i in 1::iterations
            if experiment()
                inside += one
            end
        end
        return BigFloat(4) * inside / BigFloat(iterations)
    end

end 

module taylor 



    function remainder(iterations::Int)
        if (n + 1) % 2 == 0
            BigFloat(1) / (BigFloat(2) * BigFloat(iterations) + BigFloat(3))
        else 
            BigFloat(-1) / (BigFloat(2) * BigFloat(iterations) + BigFloat(3))
        end
    end

    function calc(iterations::Int)
        res = BigFloat(0) 
        for i in 0::iterations
            el = BigFloat(1) / (BigFloat(2 * i) + BigFloat(1))
            if i % 2 == 0
                res += el
            else 
                res -= el
            end
        end
        return res
    end

end

module Viete


    function calc(iterations::Int)
        ak = BigFloat(0)
        pot = BigFloat(2)
        for i in 1::(iterations-1)
            ak = sqrt(BigFloat(2) + ak)
            pot = pot * BigFloat(2)
        end
        return pot * sqrt(BigFloat(2) - ak)
    end
    
end