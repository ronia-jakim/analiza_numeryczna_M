using Base
setprecision(BigFloat, 1000000)
using Plots


module chudnowsky
# https://en.wikipedia.org/wiki/Chudnovsky_algorithm
    function fact(x::Int)
        res = BigFloat(1.0)
        for i in 1:x
            res *= BigFloat(i)
        end
        return res
    end


    function calc(iterations::Int)
        c = BigFloat(426880) * sqrt(BigFloat(10005))
        l = BigFloat(13591409)
        x = BigFloat(1)
        m = BigFloat(1)
        k = BigFloat(-6)
        res = l * m / x
        for q in 0:iterations
            l = l + BigFloat(545140134)
            x = x * BigFloat(-262537412640768000)
            k = k + BigFloat(12)
            m = m * ((k)^(BigFloat(3))  - BigFloat(16) * k) / (BigFloat(q+1) ^ (BigFloat(3)))
            res += l * m / x
        end
        return c / res
    end

    function calc_steps(iterations::Int)
        
        c = BigFloat(426880) * sqrt(BigFloat(10005))
        l = BigFloat(13591409)
        x = BigFloat(1)
        m = BigFloat(1)
        k = BigFloat(-6)
        res = l * m / x
        results = Vector{BigFloat}([res])
        for q in 0:iterations
            l = l + BigFloat(545140134)
            x = x * BigFloat(-262537412640768000)
            k = k + BigFloat(12)
            m = m * ((k)^(BigFloat(3))  - BigFloat(16) * k) / (BigFloat(q+1) ^ (BigFloat(3)))
            res += l * m / x
            append!(results, c/res)
        end
        return results
    end

end

module ramanujan
# https://en.wikipedia.org/wiki/Srinivasa_Ramanujan#Mathematical_achievements
    function fact(x::Int)
        res = BigFloat(1.0)
        for i in 1:x
            res *= BigFloat(i)
        end
        return res
    end

    function calc(iterations::Int)
        res = BigFloat(0)
        for q in 0:iterations
            a = fact(4 * q) * (BigFloat(1103) + BigFloat(26390) * q)
            b = fact(q) ^ BigFloat(4) * BigFloat(396) ^ (BigFloat(4) * q)
            res += a / b
        end
        return BigFloat(1.0) / ((BigFloat(2) * sqrt(BigFloat(2)) * res) / BigFloat(9801))
    end

    function calc_steps(iterations::Int)
        res = BigFloat(0)
        results = Vector{BigFloat}([])
        cons = BigFloat(2) * sqrt(BigFloat(2)) / BigFloat(9801)
        for q in 0:iterations
            a = fact(4 * q) * (BigFloat(1103) + BigFloat(26390) * q)
            b = fact(q) ^ BigFloat(4) * BigFloat(396) ^ (BigFloat(4) * q)
            res += a / b
            append!(results, BigFloat(1.0) / (res * cons))
        end
        return results
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
        for i in 1:iterations
            if experiment()
                inside += one
            end
        end
        return BigFloat(4) * inside / BigFloat(iterations)
    end

    function calc_steps(iterations::Int)
        inside = BigFloat(0)
        one = BigFloat(1)
        results = Vector{BigFloat}([])
        for i in 1:iterations
            if experiment()
                inside += one
            end
            append!(results, BigFloat(4) * inside / BigFloat(iterations))
        end
        return results
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

    function calc_steps(iterations::Int)
        results = Vector{BigFloat}([])
        res = BigFloat(0) 
        for i in 0:iterations
            el = BigFloat(1) / (BigFloat(2 * i) + BigFloat(1))
            if i % 2 == 0
                res += el
            else 
                res -= el
            append!(results, res * BigFloat(4))
            end
        end
        return results
    end

    function calc(iterations::Int)
        res = BigFloat(0) 
        for i in 0:iterations
            el = BigFloat(1) / (BigFloat(2 * i) + BigFloat(1))
            if i % 2 == 0
                res += el
            else 
                res -= el
            end
        end
        return res * BigFloat(4)
    end

end

module viete


    function calc(iterations::Int)
        ak = BigFloat(0)
        pot = BigFloat(2)
        for i in 1:iterations
            ak = sqrt(BigFloat(2) + ak)
            pot = pot * BigFloat(2)
        end
        return pot * sqrt(BigFloat(2) - ak)
    end

    function calc_steps(iterations::Int)
        ak = BigFloat(0)
        pot = BigFloat(2)
        results = Vector{BigFloat}([])
        for i in 1:iterations
            ak = sqrt(BigFloat(2) + ak)
            pot = pot * BigFloat(2)
            append!(results, pot * sqrt(BigFloat(2) - ak))
        end
        return results
    end

end


module gauss_legendre

    function  calc(iterations::Int)
        a = BigFloat(1)
        b = BigFloat(1) / sqrt(BigFloat(2))
        t = BigFloat(1) / BigFloat(4)
        p = BigFloat(1)
        for i in 1:iterations
            an = (a + b) / BigFloat(2)
            b = sqrt(a * b)
            t = t - p * (a - an) * (a - an)
            p = BigFloat(2) * p
            a = an
            
        end
        return (a + b) * (a + b) / (BigFloat(4) * t)
    end

    function  calc_steps(iterations::Int)
        a = BigFloat(1)
        b = BigFloat(1) / sqrt(BigFloat(2))
        t = BigFloat(1) / BigFloat(4)
        p = BigFloat(1)
        results = Vector{BigFloat}([])
        for i in 1:iterations
            an = (a + b) / BigFloat(2)
            b = sqrt(a * b)
            t = t - p * (a - an) * (a - an)
            p = BigFloat(2) * p
            a = an
            append!(results, (a + b) * (a + b) / (BigFloat(4) * t))
            
        end
        return results
    end

end

function abs_error(res::BigFloat, correct::BigFloat)
    return abs(res - correct)
end

function rel_error(abse::BigFloat, correct::BigFloat)
    return abs(abse / correct)
end


function calc_abs(lst::Vector{BigFloat}, bib::BigFloat)
    res = Vector{BigFloat}([])
    for x in lst
        append!(res, [abs_error(x, bib)])
    end
    return res
end

function calc_rel(lst::Vector{BigFloat}, bib::BigFloat)
    res = Vector{BigFloat}([])
    for x in lst
        append!(res, [rel_error(x, bib)])
    end
    return res
end
function print_err_graph(ab::Vector{BigFloat}, re::Vector{BigFloat})
    x = 1:length(ab)
    p = plot(x, ab, title = "Błąd przybliżenia w zależności\n od iteracji", xlabel = "liczba iteracji", label = "błąd bezwzgledny")
    plot!(p, x, re, label="błąd względny")
    # plot!(p, x, re, label="błąd względny")
    savefig(p, "err_plot")
    display(p)
end


function print_rel_err_log_graph(rel::Vector{BigFloat})
    relog = map(x -> log(abs(x)), rel)
    x = 1:length(relog)
    p = plot(x, relog, title = "Logarytm z wartości bezwględnej: błędu względnego \n w zależności od liczby iteracji", xlabel = "liczba iteracji", label = "log |błąd wzgledny|")
    display(p)
    savefig(p, "rel_log_plot")
end

# fajne pi https://julialang.org/blog/2017/03/piday/ 
function main()
    results = gauss_legendre.calc_steps(50)
    ab = calc_abs(results, BigFloat(pi))
    rel = calc_rel(ab, BigFloat(pi))

    #println(results)
    print_rel_err_log_graph(rel)
    print_err_graph(ab, rel)
end

main()