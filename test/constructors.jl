@testset "constructors" begin
    mat = HomalgMatrix(1:6, 2, 3, ZZ)
    @test NumberRows(mat) == 2
    @test NumberColumns(mat) == 3

    idmat = HomalgIdentityMatrix(3, ZZ)
    @test NumberRows(idmat) == 3
    @test NumberColumns(idmat) == 3

    zeromat = HomalgZeroMatrix(3, 2, ZZ)
    @test NumberRows(zeromat) == 3
    @test NumberColumns(zeromat) == 2

    rowvector = HomalgRowVector(1:5, 5, ZZ)
    @test NumberRows(rowvector) == 1
    @test NumberColumns(rowvector) == 5

    rowvector2 = HomalgRowVector(1:5, ZZ)
    @test NumberRows(rowvector2) == 1
    @test NumberColumns(rowvector2) == 5

    columnvector = HomalgColumnVector(1:5, 5, ZZ)
    @test NumberRows(columnvector) == 5
    @test NumberColumns(columnvector) == 1

    columnvector2 = HomalgColumnVector(1:5, ZZ)
    @test NumberRows(columnvector2) == 5
    @test NumberColumns(columnvector2) == 1

    diagonalmat = HomalgDiagonalMatrix(1:5, ZZ)
    @test NumberRows(diagonalmat) == 5
    @test NumberColumns(diagonalmat) == 5

    m_zz = RandomMatrix(2, 3, ZZ)
    @test NumberRows(m_zz) == 2
    @test NumberColumns(m_zz) == 3

    m_qq = RandomMatrix(2, 3, QQ)
    @test NumberRows(m_qq) == 2
    @test NumberColumns(m_qq) == 3

end
