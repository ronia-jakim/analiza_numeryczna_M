using Base
setprecision(BigFloat, 17000)
using Plots


module geometry3

    function calc(iterations::Int)
        L = BigFloat(6) * sqrt(BigFloat(3))
        l = BigFloat(3) * sqrt(BigFloat(3))
        for i in 1:iterations
            L = BigFloat(2) * (L * l) / (L + l)
            l = sqrt(L * l)
        end
        return (L + l) / BigFloat(4)
    end

    function calc_steps(iterations::Int)
        L = BigFloat(6) * sqrt(BigFloat(3))
        l = BigFloat(3) * sqrt(BigFloat(3))
        results = Vector{BigFloat}([])
        for i in 1:iterations
            L = BigFloat(2) * (L * l) / (L + l)
            l = sqrt(L * l)
            append!(results, (L + l) / BigFloat(4))
        end
        return results
    end
end

module geometry4

    function calc(iterations::Int)
        L = BigFloat(8) 
        l = BigFloat(4) * sqrt(BigFloat(2))
        for i in 1:iterations
            L = BigFloat(2) * (L * l) / (L + l)
            l = sqrt(L * l)
        end
        return (L + l) / BigFloat(4)
    end

    function calc_steps(iterations::Int)
        L = BigFloat(8) 
        l = BigFloat(4) * sqrt(BigFloat(2))
        results = Vector{BigFloat}([])
        for i in 1:iterations
            L = BigFloat(2) * (L * l) / (L + l)
            l = sqrt(L * l)
            append!(results, (L + l) / BigFloat(4))
        end
        return results
    end
end

module chudnowsky

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
        return d <= 1
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
            end
            append!(results, res * BigFloat(4))
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
        a = sqrt(BigFloat(2))
        res = BigFloat(1)
        for i in 1:iterations
            res = res * a / BigFloat(2)
            a = sqrt(BigFloat(2) + a)
        end
        return BigFloat(2) / res
    end

    function calc_steps(iterations::Int)
        a = sqrt(BigFloat(2))
        res = BigFloat(1)
        results = Vector{BigFloat}([])
        for i in 1:iterations
            res = res * a / BigFloat(2)
            a = sqrt(BigFloat(2) + a)
            append!(results, BigFloat(2) / res)
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


function log_error_graph_gen(iterations::Int, method, name::String, file_name::String, prec::BigFloat)
    results = method(iterations)
    ab = calc_abs(results, BigFloat(pi))
    rel = calc_rel(ab, BigFloat(pi))
    relog = map(x -> -1 * log(abs(x)), rel)
    x = 1:length(relog)
    prec_lim = Vector{BigFloat}([])
    for x in 1:length(relog)
        append!(prec_lim, [prec])
    end
    p = plot(x, [relog prec_lim], title = "Wykres zbieżności metody:\n" * name, xlabel = "liczba iteracji", label = ["ilość cyfr znaczących" "granica precyzji"], legend=:outerbottom)
    savefig(p, file_name  * "_log_error.png")
    display(p)
end

function estimate_precision(iterations::Int, method)
    results = method(iterations)
    ab = calc_abs(results, BigFloat(pi))
    rel = calc_rel(ab, BigFloat(pi))
    relog = map(x -> -1 * log(abs(x)), rel)
    maks = BigFloat(0)
    for i in 1:length(relog)
        maks = max(maks, relog[i])
    end
    return maks
end

function convergence_experiment(func, iterations::Int, p::Int, name::String, file_name::String)
    results = func(iterations)
    ab = calc_abs(results, BigFloat(pi))
    ratios = Vector{BigFloat}([])
    for i in 1:(iterations-1)
        ratio = ab[i+1] / (ab[i] ^ (BigFloat(p)))
        append!(ratios, ratio)
    end
    x = 1:length(ratios)
    p = plot(x, ratios, title="Iloraz kolejnych błędów bezwzglednych metodą:\n" * name, label="iloraz", xlabel = "iteracja")
    savefig(p, file_name  * "_error_ratio.png")
    display(p)
end



function main()
    z = estimate_precision(10000, viete.calc_steps)
    log_error_graph_gen(10000, geometry3.calc_steps, "Przybliżania wielokątami", "geo3", z)
    log_error_graph_gen(10000, montecarlo.calc_steps, "Monte carlo", "monte_carlo", z)
    log_error_graph_gen(10000, taylor.calc_steps, "Wielomian Taylora", "taylor", z)
    log_error_graph_gen(21, gauss_legendre.calc_steps, "Gauss-Legendre'a", "gauss_legendre", z)
    log_error_graph_gen(450, chudnowsky.calc_steps, "Algorytm Chudnovskych", "chudnowsky", z)
    log_error_graph_gen(10000, viete.calc_steps, "Algorytm Viete'a", "viete", z)
    log_error_graph_gen(740, ramanujan.calc_steps, "Wzór Srinivasa Ramanujana", "ramanujan", z)

    convergence_experiment(taylor.calc_steps,10000, 1, "Taylor'a", "taylor")
    convergence_experiment(chudnowsky.calc_steps,450, 1, "Chudnovskych", "chudnowsky")
    convergence_experiment(viete.calc_steps,8000, 1, "Viete'a", "viete")
    convergence_experiment(geometry3.calc_steps,1000, 1, "Przybliżania wielokątami", "geo3")
    convergence_experiment(ramanujan.calc_steps, 450, 1, "Srinivasa Ramanujana", "ramanujan")

    convergence_experiment(gauss_legendre.calc_steps, 21, 2, "Gauss-Legendre'a", "gauss_legendre")
    convergence_experiment(montecarlo.calc_steps, 1000, 1, "Monte carlo", "monte_carlo")

end

main()