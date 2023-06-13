using MatricesForHomalg
using Test
using Documenter

DocMeta.setdocmeta!(MatricesForHomalg, :DocTestSetup, :(using MatricesForHomalg); recursive = true)

include("constructors.jl")
include("properties.jl")
include("attributes.jl")
include("testmanual.jl")
