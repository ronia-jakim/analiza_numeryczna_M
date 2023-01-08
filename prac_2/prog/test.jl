function householder!(x)
    x[1] = x[1] + sign(x[1]) .* norm(x)
    x ./= norm(x);
end

function tridiag_qr(T)
    Q = eye(size(T)...)
    R = copy(T)

    for i in 1:(size(R, 1) - 1)
        u = householder!(R[i:i+1, i])
        M = u * u'

        for j in 1:size(R, 2)
            # 2 optimized matrix multiplications
            # equivalent to R[i:i+1, j] -= 2 .* M * R[i:i+1, j] 
            tmp = 2 .* (M[1, 1] .* R[i, j] + M[1, 2] .* R[i+1, j])
            R[i+1, j] -= 2 .* (M[2, 1] .* R[i, j] + M[2, 2] .* R[i+1, j])
            R[i,j] -= tmp

            # similar to Q[i:i+1, j] -= 2 .* M * R[i:i+1, j], except all transposed
            tmp = 2 .* (M[1, 1] .* Q[j, i] + M[1, 2] .* Q[j, i+1])
            Q[j, i+1] -= 2 .* (M[2, 1] .* Q[j, i] + M[2, 2] .* Q[j, i+1])
            Q[j, i] -= tmp
        end
    end
    Q, R
end

function rand_tridiag(size)
    full(Tridiagonal(rand(size-1), rand(size), rand(size-1)))
end


function main()
    SIZE = 1500
    TRIALS = 20

    subsup = rand(SIZE - 1)
    diagonal = rand(SIZE)

    tridiag = Tridiagonal(subsup, diagonal, subsup)
    T = full(tridiag)

    for i=1:TRIALS
        println("$(@elapsed tridiag_qr(T))")
        #println("$(@elapsed qr(T))")
    end
end

function test()
    T = rand_tridiag(5)
    show(qr(T))
    println()
    show(tridiag_qr(T))
end

main()