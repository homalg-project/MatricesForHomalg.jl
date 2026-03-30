using Singular
using MatricesForHomalg
using Test

@testset "SyzygiesOfRows for Singular matrices" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    # A simple 3x2 matrix
    A = Singular.zero_matrix(R, 3, 2)
    A[1, 1] = x
    A[1, 2] = y
    A[2, 1] = x^2
    A[2, 2] = y^2
    A[3, 1] = x * y
    A[3, 2] = x * y

    S = SyzygiesOfRows(A)

    # S * A must be zero
    @test iszero(S * A)

    # The syzygy matrix should have 3 columns (matching rows of A)
    @test Singular.ncols(S) == 3
end

@testset "SyzygiesOfColumns for Singular matrices" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    # A 2x3 matrix
    A = Singular.zero_matrix(R, 2, 3)
    A[1, 1] = x
    A[1, 2] = x^2
    A[1, 3] = x * y
    A[2, 1] = y
    A[2, 2] = y^2
    A[2, 3] = x * y

    S = SyzygiesOfColumns(A)

    # A * S must be zero
    @test iszero(A * S)

    # The syzygy matrix should have 3 rows (matching columns of A)
    @test Singular.nrows(S) == 3
end

@testset "SyzygiesOfRows - identity matrix has no syzygies" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    I = Singular.identity_matrix(R, 3)
    S = SyzygiesOfRows(I)

    @test iszero(S)
end

@testset "SyzygiesOfColumns - identity matrix has no syzygies" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    I = Singular.identity_matrix(R, 3)
    S = SyzygiesOfColumns(I)

    @test iszero(S)
end

@testset "Relative SyzygiesOfRows for Singular matrices" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    A = Singular.zero_matrix(R, 2, 2)
    A[1, 1] = x
    A[1, 2] = y
    A[2, 1] = x^2
    A[2, 2] = y^2

    N = Singular.zero_matrix(R, 1, 2)
    N[1, 1] = x
    N[1, 2] = y

    K = SyzygiesOfRows(A, N)

    # K * A should be in the row span of N,
    # i.e., K * A + L * N = 0 for some L,
    # equivalently: each row of K*A is a multiple of the row of N.
    @test Singular.ncols(K) == 2  # same as nrows(A)
end

@testset "Relative SyzygiesOfColumns for Singular matrices" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    A = Singular.zero_matrix(R, 2, 2)
    A[1, 1] = x
    A[1, 2] = x^2
    A[2, 1] = y
    A[2, 2] = y^2

    N = Singular.zero_matrix(R, 2, 1)
    N[1, 1] = x
    N[2, 1] = y

    K = SyzygiesOfColumns(A, N)

    # A * K should be in the column span of N
    @test Singular.nrows(K) == 2  # same as ncols(A)
end

@testset "SyzygiesOfRows via ideal syzygies" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    # Row matrix (1x3) - equivalent to ideal syzygy
    A = Singular.zero_matrix(R, 1, 3)
    A[1, 1] = x^2 * y
    A[1, 2] = x * y^2
    A[1, 3] = x^2 * y^2

    S = SyzygiesOfRows(A)

    # the result should be trivial: a 0x1 matrix since 1 row has no left kernel
    # (a single row vector's left kernel is trivially 0 or identity)
    @test iszero(S * A)
end

@testset "Consistency: SyzygiesOfRows = transpose of SyzygiesOfColumns of transpose" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    A = Singular.zero_matrix(R, 3, 2)
    A[1, 1] = x
    A[1, 2] = y
    A[2, 1] = x^2
    A[2, 2] = y^2
    A[3, 1] = x * y
    A[3, 2] = x * y

    S_rows = SyzygiesOfRows(A)
    S_cols = SyzygiesOfColumns(Singular.transpose(A))

    # S_rows * A = 0 and transpose(A) * S_cols = 0
    @test iszero(S_rows * A)
    @test iszero(Singular.transpose(A) * S_cols)
end

@testset "BasisOfColumns for Singular matrices" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    A = Singular.zero_matrix(R, 2, 3)
    A[1, 1] = x;  A[1, 2] = x^2; A[1, 3] = x * y
    A[2, 1] = y;  A[2, 2] = y^2; A[2, 3] = x * y
    B = BasisOfColumns(A)

    # B should have same number of rows as A
    @test Singular.nrows(B) == 2
    # Column rank should be at most ncols(A)
    @test Singular.ncols(B) <= Singular.ncols(A)
    # Each column of A should be reducible to zero modulo B
    @test iszero(DecideZeroColumns(A, B))
end

@testset "BasisOfRows for Singular matrices" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    A = Singular.zero_matrix(R, 3, 2)
    A[1, 1] = x;  A[2, 1] = x^2; A[3, 1] = x * y
    A[1, 2] = y;  A[2, 2] = y^2; A[3, 2] = x * y
    B = BasisOfRows(A)

    @test Singular.ncols(B) == 2
    @test Singular.nrows(B) <= Singular.nrows(A)
    @test iszero(DecideZeroRows(A, B))
end

@testset "DecideZeroColumns for Singular matrices" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    A = Singular.zero_matrix(R, 2, 2)
    A[1, 1] = x; A[1, 2] = y
    A[2, 1] = y; A[2, 2] = x

    # B is a column in the span of A
    B = Singular.zero_matrix(R, 2, 1)
    B[1, 1] = x + y; B[2, 1] = x + y  # = col1 + col2
    @test iszero(DecideZeroColumns(B, A))

    # Identity should reduce to zero modulo itself
    I = Singular.identity_matrix(R, 2)
    @test iszero(DecideZeroColumns(I, I))
end

@testset "DecideZeroRows for Singular matrices" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    A = Singular.zero_matrix(R, 2, 2)
    A[1, 1] = x; A[1, 2] = y
    A[2, 1] = y; A[2, 2] = x

    @test iszero(DecideZeroRows(A, A))
end

@testset "LeftDivide (2-arg) for Singular matrices" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    A = Singular.zero_matrix(R, 2, 2)
    A[1, 1] = x; A[1, 2] = y
    A[2, 1] = y; A[2, 2] = x

    B = Singular.zero_matrix(R, 2, 1)
    B[1, 1] = x^2 + y^2; B[2, 1] = R(2) * x * y

    X = LeftDivide(A, B)
    @test X != "fail"
    @test A * X == B
end

@testset "LeftDivide (2-arg) fail case" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    A = Singular.zero_matrix(R, 2, 2)
    A[1, 1] = x; A[1, 2] = y
    A[2, 1] = y; A[2, 2] = x

    B = Singular.zero_matrix(R, 2, 1)
    B[1, 1] = x^3 + R(1); B[2, 1] = y^3 + R(1)

    @test LeftDivide(A, B) == "fail"
end

@testset "RightDivide (2-arg) for Singular matrices" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    A = Singular.zero_matrix(R, 2, 2)
    A[1, 1] = x; A[1, 2] = y
    A[2, 1] = y; A[2, 2] = x

    B = Singular.zero_matrix(R, 1, 2)
    B[1, 1] = x^2 + y^2; B[1, 2] = R(2) * x * y

    X = RightDivide(B, A)
    @test X != "fail"
    @test X * A == B
end

@testset "LeftDivide (3-arg) for Singular matrices" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    # A*X + L*Y = B with known solution
    A = Singular.zero_matrix(R, 2, 2)
    A[1, 1] = x; A[1, 2] = y
    A[2, 1] = y; A[2, 2] = x

    L = Singular.zero_matrix(R, 2, 1)
    L[1, 1] = R(1); L[2, 1] = R(0)

    # B = A*[x;y] + L*[1] = [x^2+y^2+1; 2xy]
    B = Singular.zero_matrix(R, 2, 1)
    B[1, 1] = x^2 + y^2 + R(1); B[2, 1] = R(2) * x * y

    X = SafeLeftDivide(A, B, L)
    @test Singular.nrows(X) == 2  # ncols(A) rows
    @test Singular.ncols(X) == 1  # ncols(B) cols
    # B - A*X should be in column span of L
    @test iszero(DecideZeroColumns(B - A * X, L))
end

@testset "RightDivide (3-arg) for Singular matrices" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])

    A = Singular.zero_matrix(R, 2, 2)
    A[1, 1] = x; A[1, 2] = y
    A[2, 1] = y; A[2, 2] = x

    L = Singular.zero_matrix(R, 1, 2)
    L[1, 1] = R(1); L[1, 2] = R(0)

    # B = [x;y]^T * A + [1] * L = [x^2+y^2+1, 2xy]
    B = Singular.zero_matrix(R, 1, 2)
    B[1, 1] = x^2 + y^2 + R(1); B[1, 2] = R(2) * x * y

    X = SafeRightDivide(B, A, L)
    @test Singular.nrows(X) == 1
    @test Singular.ncols(X) == 2
    # B - X*A should be in row span of L
    @test iszero(DecideZeroRows(B - X * A, L))
end
