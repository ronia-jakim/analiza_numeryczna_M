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
    matrix = zeros(BigFloat, n, n+1)
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
    Qt = diagm(ones(BigFloat, m))
    R = float(copy(A))
    for k in 1:n
        z = R[k:m,k]
        w = [ -sign(z[1])*norm(z) - z[1]; -z[2:end] ]
        nrmw = norm(w)
        if nrmw < eps() continue; end    # nie chcemi dzielić przez 0, ale czemu?
        v = w / nrmw;
        
        # aplikowanie refleksi do wierszy
        for j in 1:n
            R[k:m,j] -= v*( BigFloat(2)*(v'*R[k:m,j]) )
        end
        # aplikowanie refleksi do kolumn
        for j in 1:m
            Qt[k:m,j] -= v*( BigFloat(2)*(v'*Qt[k:m,j]) )
        end
    end
    return Qt',triu(R)
  end

  function backsub(U,b)

    n = size(U,1)
    x = zeros(BigFloat, n)
    x[n] = b[n]/U[n,n]
    for i = n-1:-1:1
        s::BigFloat = sum( U[i,j]*x[j] for j=i+1:n )
        x[i] = ( b[i] - s ) / U[i,i]
    end
    return x
  end

  function calc_householder(A::Array{BigFloat}, b::Vector{BigFloat})
    n,m = size(A)
    Q,R = qrfact(A)
    Q = transpose(Q)
    c = zeros(BigFloat, n)
    for i in 1:n
      for j in 1:n
        c[i] += Q[i,j] * b[j]
      end
    end
    x = backsub(R, c)
    return x

  end
end
#histrioa zbugowana


function gen_test(n, x)
  A = rand(BigFloat, n, n)
  b = zeros(BigFloat,n)
  for i in 1:n
    for j in 1:n
      b[i] += A[i, j] * x[j]
    end
  end
  return A, b, x
end

function gen_b(A, x)
  n = size(A, 1)
  b = zeros(BigFloat, n)
  for i in 1:n
    for j in 1:n
      b[i] += A[i, j] * x[j]
    end
  end
  return b
end

function calc_error(A, b, x)
  nb = gen_b(A, x)
  return -log(norm(nb .- b))
end

function matrix_norm(A)
  #n,m = size(A)
  #x = ones(BigFloat, n)
  #x = gen_b(A, x)
  #return -log(norm(x))
  return -log(norm(A))
end

function calc_house_err1(A)
  Q,R = householder.qrfact(A)
  T = A - (Q * R)
  return matrix_norm(T)
end

function calc_house_err2(A)
  Q,R = householder.qrfact(A)
  T = Q' * Q - I
  return matrix_norm(T)
end

function calc_house_err3(A)
  Q,R = householder.qrfact(A)
  T = Q' * A - R
  return matrix_norm(T)
end

function calc_house_err4(A)
  Q,R = householder.qrfact(A)
  T = A * inv(R) - Q
  return matrix_norm(T)
end

setprecision(64)

