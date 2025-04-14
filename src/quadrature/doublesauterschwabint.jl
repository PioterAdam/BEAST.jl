# -------- used packages

using LinearAlgebra
using StaticArrays

# -------- included files
include("gqlog.jl")


abstract type SauterSchwabStrategy1d end

struct CommonEdge{A} <: SauterSchwabStrategy1d
    qps::A
end
struct CommonVertex{A} <: SauterSchwabStrategy1d
    qps::A
end



"""
	(::CommonEdge)(f, ξ, η)

Regularizing coordinate transform for parametrization on the unit line: [0,1] ↦ Γ.

Common face case.
"""


#=
# Dutch version
function (::CommonEdge)(f, u, v)

    return (1 - v) *
            (
            f( v + (1 - v) * u, (1 - v) * u )  +  
            f( (1 - v) * u, v + (1 - v) * u )
            )
end
=#

# our version
function (::CommonEdge)(f, w, z)

    return (1-z) *
            (
            f( (1 - w) * (1 - z), (1 - w) * (1 - z) + z )  +
            f( 1 - (1 - w) * (1 - z),  1 - (1 - w) * (1 - z) - z)
            )
end



"""
	(::CommonVertex)(f, ξ, η)

Regularizing coordinate transform for parametrization on the unit line: [0,1] ↦ Γ.

Common vertex case.
"""

#=
#Dutch version
function (::CommonVertex)(f, u, v)

    return v * (
        f( 1 - ( 1 - u) * v, u * v) +
        f( u * v, 1 - (1 - u) * v)
        )
end
=#


# our version
function (::CommonVertex)(f, w, z)
    return z * (
        f((1 - w)*z, 1 - w*z) +
        f(1 - (1 - w)*z, w*z)
    )
end


#=
function (::CommonVertex)(f, ξ, η)

    return (f(ξ, η))
end
=#


function _JoshuasRules(n,a,b)

    x, w = generalizedquadrature(n)
    return collect(zip(x,w))
end


function sauterschwab_parameterized1d(integrand, strategy::SauterSchwabStrategy1d)

    qps = strategy.qps
    return sum(w1 * w2 * strategy(integrand, ξ, η) for (ξ, w1) in qps, (η, w2) in qps)
end

