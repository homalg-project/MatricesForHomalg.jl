@testset "ViewStringOfHomalgMatrix" begin
    mat = HomalgMatrix(1:6, 2, 3, ZZ)
    idmat = HomalgIdentityMatrix(3, ZZ)
    zeromat = HomalgZeroMatrix(3, 2, ZZ)
    diagmat = HomalgDiagonalMatrix(1:5, ZZ)
    emptymat = HomalgMatrix([], 0, 0, ZZ)
    emptyrows = HomalgMatrix([], 0, 3, ZZ)

    @test ViewStringOfHomalgMatrix(mat) == "<A 2 x 3 matrix over a ring>"
    @test ViewStringOfHomalgMatrix(idmat) == "<A 3 x 3 identity matrix over a ring>"
    @test ViewStringOfHomalgMatrix(zeromat) == "<A 3 x 2 zero matrix over a ring>"
    @test ViewStringOfHomalgMatrix(diagmat) == "<A diagonal 5 x 5 matrix over a ring>"
    @test ViewStringOfHomalgMatrix(emptymat) == "<A 0 x 0 identity matrix over a ring>"
    @test ViewStringOfHomalgMatrix(emptyrows) == "<A 0 x 3 zero matrix over a ring>"

    qmat = HomalgMatrix(1:6, 2, 3, QQ)
    @test ViewStringOfHomalgMatrix(qmat) == "<A 2 x 3 matrix over a ring>"
end

@testset "DisplayStringOfHomalgMatrix" begin
    mat = HomalgMatrix([1,-22,333,4444,5,-6], 2, 3, ZZ)
    @test DisplayStringOfHomalgMatrix(mat) == "[ [     1,   -22,   333 ],\n  [  4444,     5,    -6 ] ]\n"

    rowvec = HomalgRowVector(1:5, ZZ)
    @test DisplayStringOfHomalgMatrix(rowvec) == "[ [  1,  2,  3,  4,  5 ] ]\n"

    onebyone = HomalgMatrix([42], 1, 1, ZZ)
    @test DisplayStringOfHomalgMatrix(onebyone) == "[ [  42 ] ]\n"

    qmat = HomalgMatrix([QQ(1,2), QQ(0), QQ(0), QQ(1,2)], 2, 2, QQ)
    @test DisplayStringOfHomalgMatrix(qmat) == "[ [  1/2,    0 ],\n  [    0,  1/2 ] ]\n"

    emptymat = HomalgMatrix([], 0, 0, ZZ)
    @test DisplayStringOfHomalgMatrix(emptymat) == "(an empty 0 x 0 matrix)\n"

    emptycols = HomalgMatrix([], 2, 0, ZZ)
    @test DisplayStringOfHomalgMatrix(emptycols) == "(an empty 2 x 0 matrix)\n"
end

@testset "DisplayHomalgMatrix" begin
    mat = HomalgMatrix([1,2,3,4], 2, 2, ZZ)
    pipe = Pipe()
    redirect_stdout(pipe) do
        DisplayHomalgMatrix(mat)
    end
    close(pipe.in)
    @test read(pipe.out, String) == DisplayStringOfHomalgMatrix(mat)
end
