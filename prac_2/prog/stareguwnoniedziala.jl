module house2
  using LinearAlgebra
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

  function calc_householder(A, b)
    n,m = size(A)
    Q,R = qr(A)
    Q = transpose(Q)
    matrix = zeros(n,n+1)
    for i in 1:n
      for j in 1:n
        matrix[i,j] = Q[i,j]
      end
    end

    for i in 1:n
      matrix[i,n+1] = b[i]
    end
    return back_substitution(matrix,n)
  end
end