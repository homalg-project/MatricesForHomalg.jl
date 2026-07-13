@testset "properties" begin
    mat = HomalgMatrix(1:6, 2, 3, ZZ)
    idmat = HomalgIdentityMatrix(3, ZZ)
    zeromat = HomalgZeroMatrix(2, 3, ZZ)
    emptymat = HomalgMatrix([], 0, 0, ZZ)
    emptyrows = HomalgMatrix([], 0, 3, ZZ)
    emptycols = HomalgMatrix([], 2, 0, ZZ)
    symmetricmat = HomalgMatrix([1,2,3,2,4,5,3,5,6], 3, 3, ZZ)
    uppertriangularmat = HomalgMatrix([1,2,3,0,4,5,0,0,6], 3, 3, ZZ)

    @test IsOne(mat) == false
    @test IsOne(idmat) == true
    @test IsOne(emptymat) == true
    @test IsOne(emptyrows) == false
    @test IsOne(emptycols) == false


    @test IsZero(mat) == false
    @test IsZero(zeromat) == true
    @test IsZero(emptymat) == true
    @test IsZero(emptyrows) == true
    @test IsZero(emptycols) == true

    @test IsEmptyMatrix(mat) == false
    @test IsEmptyMatrix(emptymat) == true
    @test IsEmptyMatrix(emptyrows) == true
    @test IsEmptyMatrix(emptycols) == true

    @test emptymat * emptyrows == emptyrows
    @test emptycols * emptymat == emptycols
    @test emptycols * emptyrows == zeromat

    @test IsSymmetricMatrix(mat) == false
    @test IsSymmetricMatrix(symmetricmat) == true

    @test HasHasInvariantBasisProperty(ZZ) == true
    @test HasHasInvariantBasisProperty(QQ) == true
    @test HasInvariantBasisProperty(ZZ) == true
    @test HasInvariantBasisProperty(QQ) == true
    @test HasHasInvariantBasisProperty(Nemo.residue_ring(ZZ, 5)[1]) == false

end