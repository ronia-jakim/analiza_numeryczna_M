using LinearAlgebra


module gauss
  function gauss_elimination(matrix, n)
    for i in 1:(n-1)
      pivot = matrix[i, i]
      for j in (i+1):n
        base = matrix[j, i] / pivot
        matrix[j,:] = matrix[j,:] - (base .* matrix[i, :])
      end
    end
    return matrix
  end

  function back_substitution(A, n)
    b = A[:, end]
    b[end] /= A[end, end-1]
    for i in n-1:-1:1
      pivot = A[i,i]
      b[i] -= sum(A[i,i+1:end-1] .* b[i+1:end])
      b[i] /= pivot
    end
    return b
  end

  function calc_gauss(A, b) 
    n = size(A, 1)
    matrix = zeros(Float64, n, n+1)
    for i in 1:n
      for j in 1:n
        matrix[i, j] = A[i, j]
      end
    end
    for i in 1:n
      matrix[i, n+1] = b[i]
    end
    nA = gauss_elimination(matrix, n)
    return back_substitution(nA, n)
  end

end

module householder
  using LinearAlgebra
  function qrfact(A)
    m,n = size(A)
    Qt = diagm(ones(m))
    R = float(copy(A))
    for k in 1:n
        z = R[k:m,k]
        w = [ -sign(z[1])*norm(z) - z[1]; -z[2:end] ]
        nrmw = norm(w)
        if nrmw < eps() continue; end    # nie chcemi dzielić przez 0, ale czemu?
        v = w / nrmw;
        
        # aplikowanie refleksi do wierszy
        for j in 1:n
            R[k:m,j] -= v*( 2*(v'*R[k:m,j]) )
        end
        # aplikowanie refleksi do kolumn
        for j in 1:m
            Qt[k:m,j] -= v*( 2*(v'*Qt[k:m,j]) )
        end
    end
    return Qt',triu(R)
  end

  function backsub(U,b)

    n = size(U,1)
    x = zeros(n)
    x[n] = b[n]/U[n,n]
    for i = n-1:-1:1
        s = sum( U[i,j]*x[j] for j=i+1:n )
        x[i] = ( b[i] - s ) / U[i,i]
    end
    return x
  end

  function calc_householder(A, b)
    Q,R = qrfact(A)
    c = Q' .* b
    x = backsub(R, c)
    return x

  end
end


A = [1.0 2.0 3.0; 2.0 7.0 5.0; 1.0 4.0 9.0]
b = [1.0 6.0 -3.0]
#print(householder.calc_householder(A, b))
#print(gauss.calc_gauss(A, b))

Q, R = householder.qrfact(A)

println(Q)
println(R)

println("====")
# QRx = b

