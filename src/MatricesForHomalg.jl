@doc """
Matrices for the homalg project
"""
module MatricesForHomalg

import Base: getindex

import Nemo
import Nemo: ZZ, QQ, matrix

export ZZ, QQ

## Constructors of homalg rings

TypeOfRingForHomalg = Nemo.NCRing

TypeOfFieldForHomalg = Nemo.Field

TypeOfRingElementForHomalg = Nemo.NCRingElement

include("Declarations.jl")

"""
    HomalgRingOfIntegers()

Returns the ring of (rational) integers.

```jldoctest
julia> ring = HomalgRingOfIntegers()
Integer ring
```
"""
function HomalgRingOfIntegers()::Nemo.ZZRing
    return Nemo.ZZ
end

"""
    HomalgFieldOfRationals()

Returns the field of rationals.

```jldoctest
julia> field = HomalgFieldOfRationals()
Rational field
```
"""
function HomalgFieldOfRationals()::Nemo.QQField
    return Nemo.QQ
end


"""
    RingName(ring)

Returns the name of the ring as a string.

```jldoctest
julia> ring = HomalgRingOfIntegers()
Integer ring

julia> RingName(ring)
"Z"

julia> field = HomalgFieldOfRationals()
Rational field

julia> RingName(field)
"Q"
```
"""
function RingName(field_of_rationals::Nemo.QQField)::String
    return "Q"
end

function RingName(ring_of_integers::Nemo.ZZRing)::String
    return "Z"
end

export HomalgRingOfIntegers, HomalgFieldOfRationals, RingName

## Constructors of homalg matrices

TypeOfMatrixForHomalg::Type = Union{Nemo.MatrixElem, Nemo.SetElem}

"""
    HomalgMatrix(L, r, c, R)

Construct a (r x c)-matrix over the ring R with L as list of entries

```jldoctest
julia> mat1 = HomalgMatrix([1,2,3,4,5,6], 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> mat2 = HomalgMatrix([[1,2,3],[4,5,6]], 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> mat1 == mat2
true
```
"""
function HomalgMatrix(entries, r, c, R)::TypeOfMatrixForHomalg
    return matrix(R, r, c, vcat(entries...))
end

"""
    HomalgIdentityMatrix(r, R)

Construct a (r x r)-identity matrix over the ring R

```jldoctest
julia> mat = HomalgIdentityMatrix(3, ZZ)
[1   0   0]
[0   1   0]
[0   0   1]
```
"""
function HomalgIdentityMatrix(r, R)::TypeOfMatrixForHomalg
    return Nemo.identity_matrix(R, Int(r))
end

"""
    HomalgZeroMatrix(r, c, R)

Construct a (r x c)-zero matrix over the ring R

```jldoctest
julia> mat = HomalgZeroMatrix(3, 2, ZZ)
[0   0]
[0   0]
[0   0]
```
"""
function HomalgZeroMatrix(r, c, R)::TypeOfMatrixForHomalg
    return Nemo.zero_matrix(R, r, c)
end

"""
    HomalgRowVector(L, c, R)

Construct a (1 x c)-matrix over the ring R with entries in the list L

```jldoctest
julia> mat = HomalgRowVector(1:5, 5, ZZ)
[1   2   3   4   5]
```
"""
function HomalgRowVector(entries, c, R)::TypeOfMatrixForHomalg
    return HomalgMatrix(entries, 1, c, R)
end

"""
    HomalgRowVector(L, R)

Construct a (1 x c)-matrix over the ring R with entries in the list L, where c = length(L)

```jldoctest
julia> mat = HomalgRowVector(1:5, ZZ)
[1   2   3   4   5]
```
"""
function HomalgRowVector(entries, R)::TypeOfMatrixForHomalg
    return HomalgRowVector(entries, length(entries), R)
end

"""
    HomalgColumnVector(L, r, R)

Construct a (r x 1)-matrix over the ring R with entries in the list L

```jldoctest
julia> mat = HomalgColumnVector(1:5, 5, ZZ)
[1]
[2]
[3]
[4]
[5]
```
"""
function HomalgColumnVector(entries, r, R)::TypeOfMatrixForHomalg
    return HomalgMatrix(entries, r, 1, R)
end

"""
    HomalgColumnVector(L, r, R)

Construct a (r x 1)-matrix over the ring R with entries in the list L, where r = length(L)

```jldoctest
julia> mat = HomalgColumnVector(1:5, ZZ)
[1]
[2]
[3]
[4]
[5]
```
"""
function HomalgColumnVector(entries, R)::TypeOfMatrixForHomalg
    return HomalgColumnVector(entries, length(entries), R)
end

"""
    HomalgDiagonalMatrix(L, R)

Construct a diagonal matrix over the ring R using the list L of diagonal entries

```jldoctest
julia> mat = HomalgDiagonalMatrix(1:5, ZZ)
[1   0   0   0   0]
[0   2   0   0   0]
[0   0   3   0   0]
[0   0   0   4   0]
[0   0   0   0   5]
```
"""
function HomalgDiagonalMatrix(diagonal_entries, R)::TypeOfMatrixForHomalg
    return Nemo.block_diagonal_matrix(map(a->HomalgMatrix([a],1,1,R), diagonal_entries))
end

"""
    R * mat

Rewrite the matrix mat over the ring R (if possible)

```jldoctest
julia> mat = HomalgMatrix([1,2,3,4,5,6], 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> QQ * mat
[1   2   3]
[4   5   6]

julia> qmat = QQ * mat
[1   2   3]
[4   5   6]

julia> ZZ * qmat == mat
true
```
"""
Base.:*(R, mat) = Nemo.change_base_ring(R, mat)

export HomalgMatrix, HomalgIdentityMatrix, HomalgZeroMatrix, HomalgRowVector, HomalgColumnVector, HomalgDiagonalMatrix

## Properties of homalg matrices

"""
    IsOne(mat)

Return true if mat is an identity matrix, otherwise false

```jldoctest
julia> mat = HomalgIdentityMatrix(3, ZZ)
[1   0   0]
[0   1   0]
[0   0   1]

julia> IsOne(mat)
true
```
"""
function IsOne(mat)::Bool
    return NumberRows(mat) == NumberColumns(mat) && Nemo.Nemo.isone(mat) ## the precondition fixes a bug in Nemo 0.48.0
end

"""
    IsZero(mat)

Return true if mat is a zero matrix, otherwise false

```jldoctest
julia> mat = HomalgZeroMatrix(3, 2, ZZ)
[0   0]
[0   0]
[0   0]

julia> IsZero(mat)
true
```
"""
function IsZero(mat)::Bool
    return Nemo.iszero(mat)
end

"""
    IsEmptyMatrix(mat)

Return true if mat does not contain any entry, otherwise false

```jldoctest
julia> mat = HomalgMatrix([], 0, 0, ZZ)
0 by 0 empty matrix

julia> IsEmptyMatrix(mat)
true
```
"""
function IsEmptyMatrix(mat)::Bool
    return isempty(mat)
end

"""
    IsSymmetricMatrix(mat)

Return true if the matrix mat is symmetric with respect to its main diagonal, otherwise false

```jldoctest
julia> mat = HomalgMatrix([1,2,3,2,4,5,3,5,6], 3, 3, ZZ)
[1   2   3]
[2   4   5]
[3   5   6]

julia> IsSymmetricMatrix(mat)
true
```
"""
function IsSymmetricMatrix(mat)::Bool
    return Nemo.is_symmetric(mat)
end

export IsOne, IsZero, IsEmptyMatrix, IsSymmetricMatrix

## Attributes of homalg matrices

"""
    HomalgRing(mat)

Return the ring underlying the matrix mat.

```jldoctest
julia> mat = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> HomalgRing(mat)
Integer ring

julia> mat = HomalgMatrix(1:6, 2, 3, QQ)
[1   2   3]
[4   5   6]

julia> HomalgRing(mat)
Rational field
```
"""
function HomalgRing(mat)
    return Nemo.base_ring(mat)
end

"""
    NumberRows(mat)

The number of rows of the matrix mat

```jldoctest
julia> mat = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> NumberRows(mat)
2
```
"""
function NumberRows(mat)::Int64
    return Nemo.nrows(mat)
end

"""
    NumberColumns(mat)

The number of columns of the matrix mat

```jldoctest
julia> mat = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> NumberColumns(mat)
3
```
"""
function NumberColumns(mat)::Int64
    return Nemo.ncols(mat)
end

"""
    TransposedMatrix(mat)

Return the transposed matrix of mat

```jldoctest
julia> mat = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> TransposedMatrix(mat)
[1   4]
[2   5]
[3   6]
```
"""
function TransposedMatrix(mat)::TypeOfMatrixForHomalg
    return Base.transpose(mat)
end

"""
    ConvertMatrixToRow(mat)

Unfold the matrix M row-wise into a row.

```jldoctest
julia> mat = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> ConvertMatrixToRow(mat)
[2   3   4   5   6   7]

```
"""
function ConvertMatrixToRow(mat)::TypeOfMatrixForHomalg
    return UnionOfColumns(HomalgRing(mat), 1, map(i -> CertainRows(mat, i:i), 1:NumberRows(mat)))
end

"""
    ConvertMatrixToColumn(mat)

Unfold the matrix M column-wise into a column.

```jldoctest
julia> mat = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> ConvertMatrixToColumn(mat)
[2]
[4]
[6]
[3]
[5]
[7]

```
"""
function ConvertMatrixToColumn(mat)::TypeOfMatrixForHomalg
    return UnionOfRows(HomalgRing(mat), 1, map(j -> CertainColumns(mat, [j]), 1:NumberColumns(mat)))
end

"""
    RowReducedEchelonForm(mat)

Return the reduced row-echelon form and the rank of the matrix mat.

```jldoctest
julia> mat = HomalgMatrix(reverse(1:9), 3, 3, ZZ)
[9   8   7]
[6   5   4]
[3   2   1]

julia> RowReducedEchelonForm(mat)
([3 0 -3; 0 1 2; 0 0 0], 2)

julia> mat = HomalgMatrix(reverse(1:9), 3, 3, QQ)
[9   8   7]
[6   5   4]
[3   2   1]

julia> RowReducedEchelonForm(mat)
([1 0 -1; 0 1 2; 0 0 0], 2)
"""
function RowReducedEchelonForm(mat::Nemo.ZZMatrix)::Tuple{TypeOfMatrixForHomalg, Int64}
    hnf = Nemo.hnf(mat)
    rank = Nemo.rank(hnf)
    return hnf, rank
end

function RowReducedEchelonForm(mat::Nemo.QQMatrix)::Tuple{TypeOfMatrixForHomalg, Int64}
    rank, rref = Nemo.rref(mat)
    return rref, rank
end

"""
    BasisOfRows(mat)

Return the triangular form of mat

```jldoctest
julia> mat = HomalgMatrix(1:9, 3, 3, ZZ)
[1   2   3]
[4   5   6]
[7   8   9]

julia> BasisOfRows(mat)
[1   2   3]
[0   3   6]

julia> mat = HomalgMatrix(1:9, 3, 3, QQ)
[1   2   3]
[4   5   6]
[7   8   9]

julia> BasisOfRows(mat)
[1   0   -1]
[0   1    2]
```
"""
function BasisOfRows(mat)::TypeOfMatrixForHomalg
    rref, rank = RowReducedEchelonForm(mat)
    return rref[1:rank, :]
end

"""
    BasisOfColumns(mat)

Return the triangular form of mat

```jldoctest
julia> mat = HomalgMatrix(1:9, 3, 3, ZZ)
[1   2   3]
[4   5   6]
[7   8   9]

julia> BasisOfColumns(mat)
[1   0]
[1   3]
[1   6]

julia> mat = HomalgMatrix(1:9, 3, 3, QQ)
[1   2   3]
[4   5   6]
[7   8   9]

julia> BasisOfColumns(mat)
[ 1   0]
[ 0   1]
[-1   2]
```
"""
function BasisOfColumns(mat)::TypeOfMatrixForHomalg
    return TransposedMatrix(BasisOfRows(TransposedMatrix(mat)))
end

"""
    ZeroRows(mat)

Return a (possibly empty) list of positive integers. The list of zero rows of the matrix A.

```jldoctest
julia> mat = HomalgMatrix(4:9, 3, 2, ZZ)
[4   5]
[6   7]
[8   9]

julia> ZeroRows(mat)
Int64[]

julia> mat = HomalgMatrix([0, 2, 6, 0, 0, 0], 3, 2, ZZ)
[0   2]
[6   0]
[0   0]

julia> ZeroRows(mat)
1-element Vector{Int64}:
 3

julia> mat = HomalgZeroMatrix(3,3,ZZ)
[0   0   0]
[0   0   0]
[0   0   0]

julia> ZeroRows(mat)
3-element Vector{Int64}:
 1
 2
 3
```
"""
function ZeroRows(mat)::Vector{Int64}
    return filter(i -> IsZero(CertainRows(mat, [i])), 1:NumberRows(mat))
end

"""
    ZeroColumns(mat)

Return a (possibly empty) list of positive integers. The list of zero columns of the matrix A.

```jldoctest
julia> mat = HomalgMatrix(4:9, 2, 3, ZZ)
[4   5   6]
[7   8   9]

julia> ZeroColumns(mat)
Int64[]

julia> mat = HomalgMatrix([0, 2, 6, 0, 0, 0], 2, 3, ZZ)
[0   2   6]
[0   0   0]

julia> ZeroColumns(mat)
1-element Vector{Int64}:
 1

julia> mat = HomalgZeroMatrix(3,3,ZZ)
[0   0   0]
[0   0   0]
[0   0   0]

julia> ZeroColumns(mat)
3-element Vector{Int64}:
 1
 2
 3
```
"""
function ZeroColumns(mat)::Vector{Int64}
    return filter(i -> IsZero(CertainColumns(mat, [i])), 1:NumberColumns(mat))
end

"""
    FirstZeroRow(mat)

Return a positive integer of the first zero row.

```jldoctest
julia> mat = HomalgMatrix(4:9, 3, 2, ZZ)
[4   5]
[6   7]
[8   9]

julia> FirstZeroRow(mat)
4

julia> mat = HomalgMatrix([0, 2, 6, 0, 0, 0], 3, 2, ZZ)
[0   2]
[6   0]
[0   0]

julia> FirstZeroRow(mat)
3

julia> mat = HomalgMatrix([0, 0, 6, 0, 0, 0], 3, 2, ZZ)
[0   0]
[6   0]
[0   0]

julia> FirstZeroRow(mat)
1
```
"""
function FirstZeroRow(mat)::Int64
    first_zero_row = findfirst(i -> IsZero(CertainRows(mat, [i])), 1:NumberRows(mat))
    if first_zero_row === nothing
        return 1+NumberRows(mat)
    end
    return first_zero_row
end

"""
    FirstZeroColumn(mat)

Return a positive integer of the first zero column.

```jldoctest
julia> mat = HomalgMatrix(4:9, 2, 3, ZZ)
[4   5   6]
[7   8   9]

julia> FirstZeroColumn(mat)
4

julia> mat = HomalgMatrix([0, 2, 0, 0, 0, 0], 2, 3, ZZ)
[0   2   0]
[0   0   0]

julia> FirstZeroColumn(mat)
1

julia> mat = HomalgZeroMatrix(3,3,ZZ)
[0   0   0]
[0   0   0]
[0   0   0]

julia> FirstZeroColumn(mat)
1
```
"""
function FirstZeroColumn(mat)::Int64
    first_zero_column = findfirst(i -> IsZero(CertainColumns(mat, [i])), 1:NumberColumns(mat))
    if first_zero_column === nothing
        return 1+NumberColumns(mat)
    end
    return first_zero_column
end

"""
    SyzygiesOfRows(A)

Return a homalg matrix.
Let R be the ring over which A is defined (R:= HomalgRing(A)). The matrix of row syzygies SyzygiesGeneratorsOfRows(A) is a matrix whose rows span the left kernel of A,
i.e. the R-submodule of the free left module R(1xNrRows(A)) consisting of all rows X satisfying XA=0

```jldoctest
julia> A = HomalgMatrix(4:9, 3, 2, ZZ)
[4   5]
[6   7]
[8   9]

julia> S = SyzygiesOfRows(A)
[1   -2   1]

julia> S*A
[0   0]
```
"""
function SyzygiesOfRows(A)::TypeOfMatrixForHomalg
    ring = HomalgRing(A)
    nr_rows = NumberRows(A)
    nr_cols = NumberColumns(A)

    ident_mat_a = HomalgIdentityMatrix(nr_rows, ring)

    temp_mat = UnionOfColumns(ring, nr_rows, [A, ident_mat_a])

    resulting_mat = BasisOfRows(temp_mat)

    BA = resulting_mat[:, 1 : nr_cols]

    s = FirstZeroRow(BA)

    return resulting_mat[s:nr_rows, nr_cols + 1 : nr_cols + nr_rows]
end

"""
    SyzygiesOfColumns(A)

Return a homalg matrix.
Let R be the ring over which A is defined (R:=HomalgRing(A)). The matrix of column syzygies SyzygiesGeneratorsOfColumns(A) is a matrix whose columns span the right kernel of A,
i.e. the R-submodule of the free right module R(NrColumns(A)x1) consisting of all columns X satisfying AX=0

```jldoctest
julia> A = TransposedMatrix(HomalgMatrix(4:9, 3, 2, ZZ))
[4   6   8]
[5   7   9]

julia> S = SyzygiesOfColumns(A)
[ 1]
[-2]
[ 1]

julia> A*S
[0]
[0]
```
"""
function SyzygiesOfColumns(A)::TypeOfMatrixForHomalg
    return TransposedMatrix(SyzygiesOfRows(TransposedMatrix(A)))
end

"""
    RowRankOfMatrix(mat)

Returns the row rank of the matrix mat as a nonnegative integer.

```jldoctest
julia> mat1 = HomalgMatrix(4:9, 3, 2, ZZ)
[4   5]
[6   7]
[8   9]

julia> RowRankOfMatrix(mat1)
2

julia> mat2 = HomalgIdentityMatrix(3, ZZ)
[1   0   0]
[0   1   0]
[0   0   1]

julia> RowRankOfMatrix(mat2)
3
```
"""
function RowRankOfMatrix(mat)::Int64
    RowReducedEchelonForm(mat)[2]
end

"""
    ColumnRankOfMatrix(mat)

Returns the column rank of the matrix mat as a nonnegative integer.

```jldoctest
julia> mat1 = TransposedMatrix(HomalgMatrix(4:9, 3, 2, ZZ))
[4   6   8]
[5   7   9]

julia> ColumnRankOfMatrix(mat1)
2

julia> mat2 = HomalgIdentityMatrix(3, ZZ)
[1   0   0]
[0   1   0]
[0   0   1]

julia> ColumnRankOfMatrix(mat2)
3
```
"""
function ColumnRankOfMatrix(mat)::Int64
    RowRankOfMatrix(mat)
end

"""
    StringDisplay(mat)

Returns the matrix mat as a string.

```jldoctest
julia> mat = HomalgMatrix(4:9, 3, 2, ZZ)
[4   5]
[6   7]
[8   9]

julia> StringDisplay(mat)
"[4 5; 6 7; 8 9]"
```
"""
function StringDisplay(mat)::String
    string(mat)
end

export HomalgRing, NumberRows, NumberColumns, TransposedMatrix, ConvertMatrixToRow, ConvertMatrixToColumn,
    RowReducedEchelonForm, BasisOfRows, BasisOfColumns, ZeroRows, ZeroColumns, FirstZeroRow, FirstZeroColumn,
    SyzygiesOfRows, SyzygiesOfColumns, RowRankOfMatrix, ColumnRankOfMatrix, StringDisplay,
    ReducedSyzygiesOfRows, ReducedSyzygiesOfColumns

## Operations of homalg matrices

"""
    UnionOfRows(R, nr_cols, list)

Return the matrices in list stacked, where all of them have same number of columns nr_cols.

```jldoctest
julia> UnionOfRows(ZZ, 3, [])
0 by 3 empty matrix

julia> mat = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> UnionOfRows(ZZ, 3, [mat, mat])
[1   2   3]
[4   5   6]
[1   2   3]
[4   5   6]
```
"""
function UnionOfRows(R, nr_cols, list)::TypeOfMatrixForHomalg
    if length(list) == 0
        return HomalgZeroMatrix(0, nr_cols, R)
    end

    return vcat(list...)
end

"""
    UnionOfColumns(R, nr_rows, list)

Return the matrices in list augmented, where all of them have same number of rows nr_rows.
.

```jldoctest
julia> UnionOfColumns(ZZ, 2, [])
2 by 0 empty matrix

julia> mat = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> UnionOfColumns(ZZ, 2, [mat, mat])
[1   2   3   1   2   3]
[4   5   6   4   5   6]
```
"""
function UnionOfColumns(R, nr_rows, list)::TypeOfMatrixForHomalg
    if length(list) == 0
        return HomalgZeroMatrix(nr_rows, 0, R)
    end

    return hcat(list...)
end

"""
    CertainColumns(mat, list)

Return the matrix of which the i-th column is the k-th column of the homalg matrix M, where k=list[i].

```jldoctest
julia> mat = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> CertainColumns(mat, [2, 2, 1])
[2   2   1]
[5   5   4]

julia> CertainColumns(mat, [])
2 by 0 empty matrix

julia> CertainColumns(mat, 4:3)
2 by 0 empty matrix
```
"""
function CertainColumns(mat, list)::TypeOfMatrixForHomalg
    if length(list) == 0
        return HomalgZeroMatrix(NumberRows(mat), 0, HomalgRing(mat))
    end
    return mat[:, map(x->Int(x), list)]
end

"""
    CertainRows(mat, list)

Return the matrix of which the i-th row is the k-th row of the homalg matrix M, where k=list[i].

```jldoctest
julia> mat = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> CertainRows(mat, [2, 2, 1])
[4   5]
[4   5]
[2   3]

julia> CertainRows(mat, [])
0 by 2 empty matrix

julia> CertainRows(mat, 4:3)
0 by 2 empty matrix
```
"""
function CertainRows(mat, list)::TypeOfMatrixForHomalg
    if length(list) == 0
        return HomalgZeroMatrix(0, NumberColumns(mat), HomalgRing(mat))
    end
    return mat[map(x->Int(x), list), :]
end

"""
    KroneckerMat(mat1, mat2)

Return the Kronecker (or tensor) product of the two homalg matrices mat1 and mat2.

```jldoctest
julia> mat1 = HomalgMatrix(1:6, 2, 3, ZZ)
[1   2   3]
[4   5   6]

julia> mat2 = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> KroneckerMat(mat1, mat2)
[ 2    3    4    6    6    9]
[ 4    5    8   10   12   15]
[ 6    7   12   14   18   21]
[ 8   12   10   15   12   18]
[16   20   20   25   24   30]
[24   28   30   35   36   42]

julia> KroneckerMat(mat2, mat1)
[ 2    4    6    3    6    9]
[ 8   10   12   12   15   18]
[ 4    8   12    5   10   15]
[16   20   24   20   25   30]
[ 6   12   18    7   14   21]
[24   30   36   28   35   42]

```
"""
function KroneckerMat(mat1, mat2)::TypeOfMatrixForHomalg
    return Nemo.kronecker_product(mat1, mat2)
end

"""
    SafeRightDivide(B, A, L)

Returns: a homalg matrix

Same as RightDivide, but asserts that the result is not fail.

```jldoctest
julia> A = HomalgMatrix(1:9, 3, 3, ZZ)
[1   2   3]
[4   5   6]
[7   8   9]

julia> B = HomalgMatrix([3, 5, 7, 13, 16, 19, 29, 33, 37], 3, 3, ZZ)
[ 3    5    7]
[13   16   19]
[29   33   37]

julia> L = HomalgMatrix(2:10, 3, 3, ZZ)
[2   3    4]
[5   6    7]
[8   9   10]

julia> X = SafeRightDivide(B, A, L)
[0   0   -2]
[0   0   -1]
[0   0    0]

julia> Y = HomalgMatrix([1, 3, 0, 0, 4, 0, -3, 7, 0], 3, 3, ZZ)
[ 1   3   0]
[ 0   4   0]
[-3   7   0]

julia> X*A+Y*L
[ 3    5    7]
[13   16   19]
[29   33   37]

julia> B = HomalgMatrix([3, 5, 7, 0, 16, 19, 0, 33, 37], 3, 3, ZZ)
[3    5    7]
[0   16   19]
[0   33   37]

julia> SafeRightDivide(B, A, L)
ERROR: Unable to solve linear system
```
"""
function SafeRightDivide(B, A, L)::TypeOfMatrixForHomalg
    ring = HomalgRing(A)
    nr_cols = NumberColumns(A)
    nr_rows_a = NumberRows(A)
    nr_rows_b = NumberRows(B)
    nr_rows_l = NumberRows(L)

    ident_mat_b = HomalgIdentityMatrix(nr_rows_b, ring)
    zero_mat_a = HomalgZeroMatrix(nr_rows_a, nr_rows_b, ring)
    zero_mat_l = HomalgZeroMatrix(nr_rows_l, nr_rows_b, ring)

    zero_mat_b = HomalgZeroMatrix(nr_rows_b, nr_rows_a, ring)
    ident_mat_a = HomalgIdentityMatrix(nr_rows_a, ring)
    zero_mat_l_2 = HomalgZeroMatrix(nr_rows_l, nr_rows_a, ring)

    union_rows_ident_zero = UnionOfRows(ring, nr_rows_b, [ident_mat_b, zero_mat_a, zero_mat_l])
    union_rows_zero_ident = UnionOfRows(ring, nr_rows_a, [zero_mat_b, ident_mat_a, zero_mat_l_2])
    union_rows_a_b_l = UnionOfRows(ring, nr_cols, [B, A, L])

    union_mat = UnionOfColumns(ring, nr_rows_a + nr_rows_b + nr_rows_l, [union_rows_ident_zero, union_rows_a_b_l, union_rows_zero_ident])
    
    temp_mat = RowReducedEchelonForm(union_mat)[1]
    b′ = temp_mat[1:nr_rows_b, (nr_rows_b + 1) : (nr_rows_b + nr_cols)]

    if b′ != 0
        error("Unable to solve linear system")
    end
    return -temp_mat[1:nr_rows_b, (nr_rows_b + nr_cols +1) : NumberColumns(temp_mat)]
end

"""
    RightDivide(B, A, L)

Returns: a homalg matrix or fail

Let B, A and L be matrices having the same number of columns and defined over the same ring. The matrix RightDivide( B, A, L ) is a particular solution of the inhomogeneous (one sided) linear system of equations XA+YL=B
in case it is solvable (for some Y which is forgotten). Otherwise fail is returned. The name RightDivide suggests "X=BA−1 modulo L".

```jldoctest
julia> A = HomalgMatrix(1:9, 3, 3, ZZ)
[1   2   3]
[4   5   6]
[7   8   9]

julia> B = HomalgMatrix([3, 5, 7, 13, 16, 19, 29, 33, 37], 3, 3, ZZ)
[ 3    5    7]
[13   16   19]
[29   33   37]

julia> L = HomalgMatrix(2:10, 3, 3, ZZ)
[2   3    4]
[5   6    7]
[8   9   10]

julia> X = RightDivide(B, A, L)
[0   0   -2]
[0   0   -1]
[0   0    0]

julia> B = HomalgMatrix([3, 5, 7, 0, 16, 19, 0, 33, 37], 3, 3, ZZ)
[3    5    7]
[0   16   19]
[0   33   37]

julia> RightDivide(B, A, L)
"fail"
```
"""
function RightDivide(B, A, L)::Union{TypeOfMatrixForHomalg, String}
    try
        return SafeRightDivide(B, A, L)
    catch
        return "fail"
    end
end

"""
    SafeRightDivide(B, A)

Returns: a homalg matrix

Same as RightDivide, but asserts that the result is not fail.

```jldoctest
julia> A = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> B = HomalgMatrix(3:8, 3, 2, ZZ)
[3   4]
[5   6]
[7   8]

julia> SafeRightDivide(B, A)
[0    1   0]
[0    0   1]
[0   -1   2]

julia> SafeRightDivide(A, B)
[0   3   -2]
[0   2   -1]
[0   1    0]

julia> SafeRightDivide(B, B)
[0   2   -1]
[0   1    0]
[0   0    1]

julia> B = HomalgMatrix(4:9, 3, 2, ZZ)
[4   5]
[6   7]
[8   9]

julia> SafeRightDivide(A, B)
ERROR: Unable to solve linear system
```
"""
function SafeRightDivide(B, A)::TypeOfMatrixForHomalg
    return SafeRightDivide(B, A, HomalgZeroMatrix(0, NumberColumns(A), HomalgRing(A)))
end

# function SafeRightDivide(mat2, mat1)::Union{TypeOfMatrixForHomalg}
#     return Nemo.solve_left(mat1, mat2)
# end

"""
    UniqueRightDivide(B, A)

Returns: a homalg matrix

Same as SafeRightDivide, but asserts that the solution is unique.

```jldoctest
julia> A = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> B = HomalgMatrix(3:8, 3, 2, ZZ)
[3   4]
[5   6]
[7   8]

julia> RightDivide(B, A)
[0    1   0]
[0    0   1]
[0   -1   2]

julia> UniqueRightDivide(B, A)
ERROR: The inhomogeneous linear system of equations XA=B has no unique solution

julia> mat = HomalgIdentityMatrix(3, ZZ)
[1   0   0]
[0   1   0]
[0   0   1]

julia> UniqueRightDivide(mat, mat)
[1   0   0]
[0   1   0]
[0   0   1]
```
"""
function UniqueRightDivide(B, A)::Union{TypeOfMatrixForHomalg}

    if NumberRows(BasisOfRows(A)) != NumberRows(A)
        error("The inhomogeneous linear system of equations XA=B has no unique solution")
    end

    return SafeRightDivide(B, A)
end

"""
    RightDivide(B, A)

Returns: a homalg matrix or fail

Let B and A be matrices having the same number of columns and defined over the same ring.
The matrix RightDivide(B, A) is a particular solution of the inhomogeneous (one sided) linear system of equations XA=B in case it is solvable.
Otherwise fail is returned. The name RightDivide suggests "X=BA^-1".

```jldoctest
julia> A = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> B = HomalgMatrix(3:8, 3, 2, ZZ)
[3   4]
[5   6]
[7   8]

julia> RightDivide(B, A)
[0    1   0]
[0    0   1]
[0   -1   2]

julia> RightDivide(A, B)
[0   3   -2]
[0   2   -1]
[0   1    0]

julia> C = HomalgMatrix(4:9, 3, 2, ZZ)
[4   5]
[6   7]
[8   9]

julia> RightDivide(A, C)
"fail"

julia> RightDivide(C, A)
"fail"
```
"""
function RightDivide(B, A)::Union{TypeOfMatrixForHomalg, String}
    try
        return SafeRightDivide(B, A)
    catch
        return "fail"
    end
end

"""
    SafeLeftDivide(A, B, L)

Returns: a homalg matrix

Same as LeftDivide, but asserts that the result is not fail.

```jldoctest
julia> A = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> L = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> B = HomalgMatrix([5, 5, 11, 9, 17, 13], 3, 2, ZZ)
[ 5    5]
[11    9]
[17   13]

julia> SafeLeftDivide(A, B, L)
[0   0]
[0   0]

julia> B = HomalgMatrix([3, 5, 0, 9, 11, 13], 3, 2, ZZ)
[ 3    5]
[ 0    9]
[11   13]

julia> SafeLeftDivide(A, B, L)
ERROR: Unable to solve linear system
```
"""
function SafeLeftDivide(A, B, L)::TypeOfMatrixForHomalg
    return TransposedMatrix(SafeRightDivide(TransposedMatrix(B), TransposedMatrix(A), TransposedMatrix(L)))
end

"""
    LeftDivide(A, B, L)

Returns: a homalg matrix or fail

Let A, B and L be matrices having the same number of columns and defined over the same ring.
The matrix LeftDivide( A, B, L ) is a particular solution of the inhomogeneous (one sided) linear system of equations AX+LY=B
in case it is solvable (for some Y which is forgotten). Otherwise fail is returned. The name LeftDivide suggests "X=A−1B modulo L"

```jldoctest
julia> A = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> L = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> B = HomalgMatrix([5, 5, 11, 9, 17, 13], 3, 2, ZZ)
[ 5    5]
[11    9]
[17   13]

julia> LeftDivide(A, B, L)
[0   0]
[0   0]

julia> B = HomalgMatrix([3, 5, 0, 9, 11, 13], 3, 2, ZZ)
[ 3    5]
[ 0    9]
[11   13]

julia> LeftDivide(A, B, L)
"fail"
```
"""
function LeftDivide(A, B, L)::Union{TypeOfMatrixForHomalg, String}
    try
        return SafeLeftDivide(A, B, L)
    catch y
        return "fail"
    end
end

"""
    SafeLeftDivide(A, B)

Returns: a homalg matrix

Same as LeftDivide, but asserts that the result is not fail.

```jldoctest
julia> A = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> B = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> C = HomalgMatrix([1,0,0,0,0,0], 3, 2, ZZ)
[1   0]
[0   0]
[0   0]

julia> SafeLeftDivide(B, B)
[1   0]
[0   1]

julia> SafeLeftDivide(A, B)
[0   -1]
[1    2]

julia> SafeLeftDivide(C, B)
ERROR: Unable to solve linear system
```
"""
function SafeLeftDivide(A, B)::TypeOfMatrixForHomalg
    return SafeLeftDivide(A, B, HomalgZeroMatrix(NumberRows(A), 0, HomalgRing(A)))
end

#= function SafeLeftDivide(A, B)::Union{TypeOfMatrixForHomalg}
    return Nemo.solve(A, B)
end =#

"""
    UniqueLeftDivide(A, B)

Returns: a homalg matrix

Same as SafeLeftDivide, but asserts that the solution is unique.

```jldoctest
julia> A = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> B = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> UniqueLeftDivide(B, B)
[1   0]
[0   1]

julia> A = HomalgMatrix([1,2,3,0,5,6,0,0,0], 3, 3, ZZ)
[1   2   3]
[0   5   6]
[0   0   0]

julia> B = HomalgMatrix([1,2,0], 3, 1, ZZ)
[1]
[2]
[0]

julia> UniqueLeftDivide(A, B)
ERROR: The inhomogeneous linear system of equations AX=B has no unique solution
```
"""
function UniqueLeftDivide(A, B)::TypeOfMatrixForHomalg
    
    if NumberColumns(BasisOfColumns(A)) != NumberColumns(A)
        error("The inhomogeneous linear system of equations AX=B has no unique solution")
    end

    return SafeLeftDivide(A, B)
end

"""
    LeftDivide(A, B)

Returns: a homalg matrix or fail

Let A and B be matrices having the same number of rows and defined over the same ring.
The matrix LeftDivide(A, B) is a particular solution of the inhomogeneous (one sided) linear system of equations AX=B in case it is solvable.
Otherwise fail is returned. The name LeftDivide suggests "X=A^{-1}B".

```jldoctest
julia> A = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> B = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> C = HomalgMatrix([1,0,0,0,0,0], 3, 2, ZZ)
[1   0]
[0   0]
[0   0]

julia> LeftDivide(B, B)
[1   0]
[0   1]

julia> LeftDivide(A, B)
[0   -1]
[1    2]

julia> LeftDivide(C, B)
"fail"
```
"""
function LeftDivide(A, B)::Union{TypeOfMatrixForHomalg, String}
    try
        return SafeLeftDivide(A, B)
    catch y
        return "fail"
    end
end

"""
    DecideZeroRows(B, A)

    Returns: a homalg matrix

Let B and A be matrices having the same number of columns and defined over the same ring R and S be the row span of A,
i.e. the R-submodule of the free left module R(1×NrColumns(A)) spanned by the rows of A.
The result is a matrix B′ having the same shape as B, for which the i-th row B′^i is equivalent to the i-th row B^i of B modulo S,
i.e. B′^i−B^i is an element of the row span S of A.
Moreover, the row B′^i is zero, if and only if the row B^i is an element of S. So DecideZeroRows decides which rows of B are zero modulo the rows of A.

```jldoctest
julia> A = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> B = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> A′ = DecideZeroRows(B, A)
[0   1]
[0   1]
[0   1]

julia> RightDivide(B, A)
"fail"

julia> RightDivide(B - A′, A)
[0   -1   1]
[0   -2   2]
[0   -3   3]

julia> DecideZeroRows(A, A)
[0   0]
[0   0]
[0   0]

julia> C = HomalgMatrix([4, 6, 2, 2], 2, 2, ZZ)
[4   6]
[2   2]

julia> DecideZeroRows(C, A)
[0   0]
[0   0]

julia> DecideZeroRows(A, C)
[1   0]
[1   0]
[1   0]
```
"""
function DecideZeroRows(B,A)::TypeOfMatrixForHomalg
    #A,B are defined over the same ring
    ring = HomalgRing(B)

    nr_rows_A = NumberRows(A)
    nr_rows_B = NumberRows(B)

    #A,B have the same number of columns
    nr_cols = NumberColumns(B)

    null_mat_a = HomalgZeroMatrix(nr_rows_A, nr_rows_B, ring)
    ident_mat_b = HomalgIdentityMatrix(nr_rows_B, ring)

    list_of_matrices = [UnionOfRows(ring, nr_rows_B, [ident_mat_b, null_mat_a]), UnionOfRows(ring, nr_cols,[B, A])]
    union_mat = UnionOfColumns(ring, nr_rows_B + nr_rows_A,list_of_matrices)
    temp_mat = BasisOfRows(union_mat)

    return temp_mat[1:nr_rows_B, nr_rows_B+1:nr_rows_B+nr_cols]
end

"""
    DecideZeroColumns(B, A)

    Returns: a homalg matrix

Let B and A be matrices having the same number of rows and defined over the same ring R and S be the column span of A,
i.e. the R-submodule of the free right module R(NrRows(A)×1) spanned by the columns of A.
The result is a matrix B′ having the same shape as B, for which the i-th column B′_i is equivalent to the i-th column B_i of B modulo S,
i.e. B′_i−B_i is an element of the column span S of A. Moreover, the column B′_i is zero, if and only if the column B_i is an element of S.
So DecideZeroColumns decides which columns of B are zero modulo the columns of A.

```jldoctest
julia> A = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> B = HomalgMatrix([3, 1, 7, 1, 11, 1], 3, 2, ZZ)
[ 3   1]
[ 7   1]
[11   1]

julia> DecideZeroColumns(B, A)
[0   0]
[0   0]
[0   0]
```
"""
function DecideZeroColumns(B, A)::TypeOfMatrixForHomalg
    return TransposedMatrix(DecideZeroRows(TransposedMatrix(B), TransposedMatrix(A)))
end

"""
    SyzygiesOfRows(A, N)

    Returns: a homalg matrix

Let R be the ring over which A is defined. The matrix of relative row syzygies SyzygiesOfRows(A, N) is a matrix whose rows span the left kernel of A modulo N,
i.e. the R-submodule of the free left module R(1×NrRows(A)) consisting of all rows K satisfying KA+LN=0 for some row L∈R(1×NrRows(N)).

```jldoctest
julia> A = HomalgMatrix(1:9, 3, 3, ZZ)
[1   2   3]
[4   5   6]
[7   8   9]

julia> N = HomalgMatrix(2:10, 3, 3, ZZ)
[2   3    4]
[5   6    7]
[8   9   10]

julia> SyzygiesOfRows(A, N)
[1   0   2]
[0   1   2]
[0   0   3]

julia> SyzygiesOfRows(A,A)
[1   0   0]
[0   1   0]
[0   0   1]
```
"""
function SyzygiesOfRows(A, N)
    ring = HomalgRing(A)
    nr_cols = NumberColumns(A)
    union_a_n = UnionOfRows(ring, nr_cols, [A, N])
    ident_a = HomalgIdentityMatrix(NumberRows(A), ring)
    zero_n = HomalgZeroMatrix(NumberRows(N), NumberRows(A), ring)
    union_mat = UnionOfColumns(ring, NumberRows(A) + NumberRows(N), [union_a_n, UnionOfRows(ring, NumberRows(A), [ident_a, zero_n])])
    temp_mat = RowReducedEchelonForm(union_mat)[1]

    first_zero_row = FirstZeroRow(temp_mat[:, 1:nr_cols])
    first_zero_row_2 = FirstZeroRow(temp_mat)
    return temp_mat[first_zero_row : (first_zero_row_2-1), (nr_cols +1) : (nr_cols + NumberRows(A))]
end

"""
    SyzygiesOfColumns(A, N)

    Returns: a homalg matrix

Let R be the ring over which M is defined. The matrix of relative column syzygies SyzygiesOfColumns(A,N) is a matrix whose columns span the right kernel of A modulo N,
i.e. the R-submodule of the free right module R(NrColumns(A)×1)consisting of all columns K satisfying AK+NL=0 for some column L∈R(NrColumns(N)×1).

```jldoctest
julia> A = HomalgMatrix(1:6, 3, 2, ZZ)
[1   2]
[3   4]
[5   6]

julia> N = HomalgMatrix(2:7, 3, 2, ZZ)
[2   3]
[4   5]
[6   7]

julia> SyzygiesOfColumns(A, N)
[1   0]
[0   1]
```
"""
function SyzygiesOfColumns(A, N)
    return TransposedMatrix(SyzygiesOfRows(TransposedMatrix(A), TransposedMatrix(N)))
end

## ReducedSyzygiesOfRows / ReducedSyzygiesOfColumns

"""
    ReducedSyzygiesOfRows(M)

Return the reduced row syzygies of M, i.e. a matrix S with
reduced rows satisfying S * M = 0.

```jldoctest
julia> A = HomalgMatrix(4:9, 3, 2, ZZ)
[4   5]
[6   7]
[8   9]

julia> S = ReducedSyzygiesOfRows(A)
[1   -2   1]

julia> S * A
[0   0]
```
"""
function ReducedSyzygiesOfRows(M)
    return BasisOfRows(SyzygiesOfRows(M))
end

"""
    ReducedSyzygiesOfRows(M, M2)

Return the reduced relative row syzygies of M modulo M2,
i.e. a matrix K with reduced rows satisfying K * M + L * M2 = 0 for some L.
"""
function ReducedSyzygiesOfRows(M, M2)
    return BasisOfRows(SyzygiesOfRows(M, M2))
end

"""
    ReducedSyzygiesOfColumns(M)

Return the reduced column syzygies of M, i.e. a matrix S with
reduced columns satisfying M * S = 0.

```jldoctest
julia> A = TransposedMatrix(HomalgMatrix(4:9, 3, 2, ZZ))
[4   6   8]
[5   7   9]

julia> S = ReducedSyzygiesOfColumns(A)
[ 1]
[-2]
[ 1]

julia> A * S
[0]
[0]
```
"""
function ReducedSyzygiesOfColumns(M)
    return TransposedMatrix(ReducedSyzygiesOfRows(TransposedMatrix(M)))
end

"""
    ReducedSyzygiesOfColumns(M, M2)

Return the reduced relative column syzygies of M modulo M2,
i.e. a matrix K with reduced columns satisfying M * K + M2 * L = 0 for some L.
"""
function ReducedSyzygiesOfColumns(M, M2)
    return TransposedMatrix(ReducedSyzygiesOfRows(TransposedMatrix(M), TransposedMatrix(M2)))
end

export UnionOfRows, UnionOfColumns, KroneckerMat, CertainColumns, CertainRows,
    SafeRightDivide, UniqueRightDivide, RightDivide, SafeLeftDivide, UniqueLeftDivide, LeftDivide,
    DecideZeroRows, DecideZeroColumns

end
