using Base

module program
    using Plots

    function approx_ln(x::Float64, deg::Int)
        ret = zero(x)
        for i in 1:deg
            temp = x^i
            temp /= i
            if i % 2 == 0
                temp *= -1
            end
            ret += temp
        end

        return ret
    end

    function abs_error(res::Float64, ln::Float64)
        return abs(res - ln)
    end

    function rel_error(err::Float64, ln::Float64)
        return err / abs(ln)
    end

    function calc_val(x::Float64, max_deg::Int)
        res = Vector{Float64}([])

        for deg in 1:max_deg
            append!(res, [approx_ln(x-1, deg)])
        end

        return res
    end

    function calc_abs(lst::Vector{Float64}, bib::Float64)
        res = Vector{Float64}([])
        for x in lst
            append!(res, [abs_error(x, bib)])
        end
        return res
    end

    function calc_rel(lst::Vector{Float64}, bib::Float64)
        res = Vector{Float64}([])
        for x in lst
            append!(res, [rel_error(x, bib)])
        end
        return res
    end

    function main()
        b = log(0.5)

        println("BIBLIOTECZNIE: ", b)

        vv = program.calc_val(0.5, 8)

        ab = program.calc_abs(vv, b)

        re = program.calc_rel(ab, b)

        println("obliczono: ", vv)
        println("bezwzglednie: ", ab)
        println("wzglednie: ", re)

        print_err_graph(ab, re)
        print_val_graph(vv, b)
    end

    function print_err_graph(ab::Vector{Float64}, re::Vector{Float64})
        x = 1:length(ab)
        p = plot(x, ab, title = "Błąd przybliżenia w zależności od stopnia szeregu Taylora", xlabel = "stopień wielomianu", label = "błąd bezwzgledny")
        plot!(p, x, re, label="błąd względny")
        # plot!(p, x, re, label="błąd względny")
        savefig(p, "err_plot")
        display(p)
    end

    function print_val_graph(vv::Vector{Float64}, ln::Float64)
        x = 1:length(vv)
        p = plot(x, vv, title = "Wynik przybliżenia w zależności od stopnia szeregu Taylora", xlabel = "stopień wielomianu", label = "wartość przybliżona")

        maruda = Vector{Float64}([])
        for i in x
            append!(maruda, [ln])
        end

        plot!(p, x, maruda, label="wartość biblioteczna")
        # plot!(p, x, re, label="błąd względny")
        savefig(p, "val_plot")
        display(p)
    end

    function calc_err_deg(deg::Int, b::Float64, x::Float64)
        vv = program.calc_val(x, deg)

        ab = program.calc_abs(vv, b)
        re = program.calc_rel(ab, b)

        print_err_graph(ab, re)
    end

    function calc_val_deg(deg::Int, b::Float64, x::Float64)
        vv = program.calc_val(x, deg)

        print_val_graph(vv, b)
    end

end

#program.main()