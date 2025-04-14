module SauterSchwabQuadrature1d

# -------- exportet parts
# types
export SauterSchwabStrategy1d
export CommonEdge, CommonVertex

# functions
export sauterschwab_parameterized1d, _JoshuasRules, reorder

# -------- included files
include("doublesauterschwabint.jl")

end