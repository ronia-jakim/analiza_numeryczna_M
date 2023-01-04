
module gauss

  function gauss(matrix, n)
    for i in 1:(n-1)
      for j in (i+1):n
        ratio = matrix[j, i] / matrix[i, i]
        for k in 1:(n+1)
          matrix[j, k] = matrix[j, k] - ratio * matrix[i, k]
        end
      end
    end
    return matrix
  end

  function back_substitution(matrix, n)
    solutions = fill(0.0, n)
    solutions[n] = matrix[n, n+1] / matrix[n, n]
    for i in (n-1):1
      solutions[i] = matrix[i, n+1]
      for j in (i+1):n
        solutions[i] = solutions[i] - matrix[i, j] * solutions[j]
      end
      solutions[i] = solutions[i] / matrix[i, i]
    end
    return solutions
  end

end

