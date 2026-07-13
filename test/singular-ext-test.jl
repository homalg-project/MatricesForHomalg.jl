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

@testset "Ring constructors for Singular" begin
    Z = HomalgRingOfIntegersInSingular()
    @test Z isa Singular.Integers

    Q = HomalgFieldOfRationalsInSingular()
    @test Q isa Singular.Rationals

    Z3 = HomalgRingOfIntegersInSingular(3)
    @test Z3 isa Singular.N_ZnRing
end

@testset "getindex syntax for Singular rings" begin
    R1 = Singular.ZZ["x"]
    @test R1 isa Singular.PolyRing
    @test length(Singular.gens(R1)) == 1

    R2 = Singular.QQ["x", "y"]
    @test R2 isa Singular.PolyRing
    @test length(Singular.gens(R2)) == 2
end

@testset "RingName for Singular rings" begin
    @test RingName(Singular.QQ) == "Q"
    @test RingName(Singular.ZZ) == "Z"

    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    @test RingName(R) == "Q[x,y]"

    S = R / [x^2 + y^2 - 1]
    @test RingName(S) == "Q[x,y] / ( x^2 + y^2 - 1 )"
end

@testset "HasHasInvariantBasisProperty and HasInvariantBasisProperty for PolyRing" begin
    R, _ = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    @test HasHasInvariantBasisProperty(R) == true
    @test HasInvariantBasisProperty(R) == true
end

@testset "Quotient ring via Base.:/" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    S = R / [x^2 + y^2 - 1]
    @test S isa Singular.PolyRing
    @test Singular.is_quotient_ring(S)
end

@testset "Indeterminates for Singular rings" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    gens = Indeterminates(R)
    @test length(gens) == 2
    @test gens[1] == x
    @test gens[2] == y
end

@testset "String parsing (R::PolyRing)(s::String)" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    @test R("x^2 + 2*y") == x^2 + 2*y
    @test R("x*y") == x * y
end

@testset "HomalgMatrix for Singular rings" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    A = HomalgMatrix([[x, y], [x^2, y^2]], 2, 2, R)
    @test NumberRows(A) == 2
    @test NumberColumns(A) == 2
    @test A[1, 1] == x
    @test A[2, 2] == y^2
end

@testset "HomalgZeroMatrix and HomalgIdentityMatrix for Singular rings" begin
    R, _ = Singular.polynomial_ring(Singular.QQ, ["x"])
    Z = HomalgZeroMatrix(2, 3, R)
    @test IsZero(Z)
    @test NumberRows(Z) == 2
    @test NumberColumns(Z) == 3

    I = HomalgIdentityMatrix(3, R)
    @test IsOne(I)
    @test NumberRows(I) == 3
    @test NumberColumns(I) == 3
end

@testset "RandomMatrix for Singular rings" begin
    R, _ = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    m = RandomMatrix(2, 3, R)
    @test NumberRows(m) == 2
    @test NumberColumns(m) == 3
    @test HomalgRing(m) == R
end

@testset "Negation of Singular matrix" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    A = HomalgMatrix([[x, y], [x^2, y^2]], 2, 2, R)
    B = -A
    @test B[1, 1] == -x
    @test B[2, 2] == -y^2
    @test IsZero(A + B)
end

@testset "NumberRows, NumberColumns, HomalgRing, TransposedMatrix for smatrix" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    A = HomalgMatrix([[x, y], [x^2, y^2]], 2, 2, R)
    @test NumberRows(A) == 2
    @test NumberColumns(A) == 2
    @test HomalgRing(A) == R
    At = TransposedMatrix(A)
    @test NumberRows(At) == 2
    @test NumberColumns(At) == 2
    @test At[1, 2] == x^2
    @test At[2, 1] == y
end

@testset "IsZero, IsOne, IsEmptyMatrix, IsSymmetricMatrix for smatrix" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    @test IsZero(HomalgZeroMatrix(2, 3, R))
    @test !IsZero(HomalgIdentityMatrix(2, R))
    @test IsOne(HomalgIdentityMatrix(3, R))
    @test !IsOne(HomalgZeroMatrix(2, 2, R))
    @test IsEmptyMatrix(HomalgZeroMatrix(0, 3, R))
    @test IsEmptyMatrix(HomalgZeroMatrix(2, 0, R))
    @test !IsEmptyMatrix(HomalgIdentityMatrix(2, R))
    sym = HomalgMatrix([[x, y], [y, x]], 2, 2, R)
    @test IsSymmetricMatrix(sym)
    nonsym = HomalgMatrix([[x, y], [x^2, y^2]], 2, 2, R)
    @test !IsSymmetricMatrix(nonsym)
end

@testset "CertainRows and CertainColumns for smatrix" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    A = HomalgMatrix([[x, y], [x^2, y^2]], 2, 2, R)
    R1 = CertainRows(A, [1])
    @test NumberRows(R1) == 1
    @test R1[1, 1] == x

    C2 = CertainColumns(A, [2])
    @test NumberColumns(C2) == 1
    @test C2[2, 1] == y^2

    # empty selection
    @test IsEmptyMatrix(CertainRows(A, []))
    @test IsEmptyMatrix(CertainColumns(A, []))
end

@testset "UnionOfRows and UnionOfColumns for smatrix" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    A = HomalgMatrix([[x, y], [x^2, y^2]], 2, 2, R)
    r1 = CertainRows(A, [1])
    r2 = CertainRows(A, [2])
    U = UnionOfRows(R, 2, [r1, r2])
    @test NumberRows(U) == 2
    @test U[1, 1] == x
    @test U[2, 2] == y^2

    c1 = CertainColumns(A, [1])
    c2 = CertainColumns(A, [2])
    V = UnionOfColumns(R, 2, [c1, c2])
    @test NumberColumns(V) == 2
    @test V[1, 2] == y
    @test V[2, 1] == x^2

    # empty list
    @test IsEmptyMatrix(UnionOfRows(R, 2, []))
    @test IsEmptyMatrix(UnionOfColumns(R, 2, []))
end

@testset "KroneckerMat for smatrix" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    M1 = HomalgMatrix([[x, R(0)], [R(0), y]], 2, 2, R)
    M2 = HomalgIdentityMatrix(2, R)
    K = KroneckerMat(M1, M2)
    @test NumberRows(K) == 4
    @test NumberColumns(K) == 4
    @test K[1, 1] == x
    @test K[3, 3] == y
    @test IsZero(K[1, 3])
end

@testset "HomalgDiagonalMatrix for Singular rings" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    D = HomalgDiagonalMatrix([x, y, x*y], R)
    @test NumberRows(D) == 3
    @test NumberColumns(D) == 3
    @test D[1, 1] == x
    @test D[2, 2] == y
    @test D[3, 3] == x * y
    @test IsZero(D[1, 2])
end

@testset "ReducedSyzygiesOfRows and ReducedSyzygiesOfColumns for smatrix" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    A = HomalgMatrix([[x, y], [x^2, y^2], [x*y, x*y]], 3, 2, R)
    S = ReducedSyzygiesOfRows(A)
    @test iszero(S * A)

    B = HomalgMatrix([[x, x^2, x*y], [y, y^2, x*y]], 2, 3, R)
    T = ReducedSyzygiesOfColumns(B)
    @test iszero(B * T)
end

@testset "Relative ReducedSyzygiesOfRows and ReducedSyzygiesOfColumns" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    A = HomalgMatrix([[x, y], [x^2, y^2]], 2, 2, R)
    N = HomalgMatrix([[x, y]], 1, 2, R)
    K = ReducedSyzygiesOfRows(A, N)
    @test iszero(DecideZeroRows(K * A, N))

    Ac = TransposedMatrix(A)
    Nc = TransposedMatrix(N)
    Kc = ReducedSyzygiesOfColumns(Ac, Nc)
    @test iszero(DecideZeroColumns(Ac * Kc, Nc))
end

@testset "SafeLeftDivide and SafeRightDivide (2-arg) for smatrix" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    A = HomalgMatrix([[x, R(0)], [R(0), y]], 2, 2, R)
    B = HomalgMatrix([[x^2, R(0)], [R(0), y^2]], 2, 2, R)
    X = SafeLeftDivide(A, B)
    @test A * X == B

    Y = SafeRightDivide(B, A)
    @test Y * A == B
end

@testset "LeftDivide and RightDivide (3-arg) for smatrix" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    A = HomalgMatrix([[x, R(0)], [R(0), y]], 2, 2, R)
    L = HomalgMatrix([[R(1), R(0)], [R(0), R(1)]], 2, 2, R)
    B = HomalgMatrix([[x^2, R(0)], [R(0), y^2]], 2, 2, R)
    X = LeftDivide(A, B, L)
    @test X != "fail"
    @test iszero(DecideZeroColumns(B - A * X, L))

    Y = RightDivide(B, A, L)
    @test Y != "fail"
    @test iszero(DecideZeroRows(B - Y * A, L))
end

@testset "RightDivide (2-arg) fail case" begin
    R, (x, y) = Singular.polynomial_ring(Singular.QQ, ["x", "y"])
    @test RightDivide(HomalgMatrix([[R(1)]], 1, 1, R), HomalgMatrix([[x], [x]], 2, 1, R)) == "fail"
end

@testset "AssignGeneratingVariables for Singular rings" begin
    Rt, _ = Singular.polynomial_ring(Singular.QQ, ["t_singular_test"])
    AssignGeneratingVariables(Rt)
    @test isdefined(Main, :t_singular_test)
end
