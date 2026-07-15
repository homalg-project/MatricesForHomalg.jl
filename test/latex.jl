@testset "LaTeXOutputOfHomalgRingElement" begin
    @test LaTeXOutputOfHomalgRingElement(ZZ(0))    == "0"
    @test LaTeXOutputOfHomalgRingElement(ZZ(-3))   == "-3"
    @test LaTeXOutputOfHomalgRingElement(ZZ(42))   == "42"
    @test LaTeXOutputOfHomalgRingElement(QQ(0))    == "0"
    @test LaTeXOutputOfHomalgRingElement(QQ(2))    == "2"
    @test LaTeXOutputOfHomalgRingElement(QQ(1,3))  == "\\frac{1}{3}"
    @test LaTeXOutputOfHomalgRingElement(QQ(-1,4)) == "-\\frac{1}{4}"
end

@testset "LaTeXOutputOfHomalgRing" begin
    @test LaTeXOutputOfHomalgRing(ZZ) == "\\mathbb{Z}"
    @test LaTeXOutputOfHomalgRing(QQ) == "\\mathbb{Q}"
end

@testset "LaTeXOutputOfHomalgMatrix" begin
    # integer matrix: zero becomes \cdot
    mat = HomalgMatrix([0, 1, 2, 0], 2, 2, ZZ)
    @test LaTeXOutputOfHomalgMatrix(mat) == "\\left( \\begin{array}{rr}\n \\cdot & 1 \\\\\n 2 & \\cdot \n\\end{array} \\right)"

    # rational matrix: fraction and zero
    mat = HomalgMatrix([QQ(0), QQ(1, 3)], 1, 2, QQ)
    @test LaTeXOutputOfHomalgMatrix(mat) == "\\left( \\begin{array}{rr}\n \\cdot & \\frac{1}{3} \n\\end{array} \\right)"

    # rational matrix: negative fraction
    mat = HomalgMatrix([QQ(-1, 4)], 1, 1, QQ)
    @test LaTeXOutputOfHomalgMatrix(mat) == "\\left( \\begin{array}{r}\n -\\frac{1}{4} \n\\end{array} \\right)"

    # identity matrix
    mat = HomalgIdentityMatrix(2, ZZ)
    @test LaTeXOutputOfHomalgMatrix(mat) == "\\left( \\begin{array}{rr}\n 1 & \\cdot \\\\\n \\cdot & 1 \n\\end{array} \\right)"
end
