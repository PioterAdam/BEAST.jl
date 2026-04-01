# minimal_singleshot_efie_tm.jl
println("RUNNING FILE: minimal_singleshot_efie_tm.jl")
using LinearAlgebra, Random, Statistics, Printf
using BEAST, CompScienceMeshes
using AdaptiveCrossApproximation
using OhMyThreads
using ParallelKMeans
using H2Trees
using SpecialFunctions
#using H2Trees: BoundingBallTree

import BEAST: kernelvals, KernelValsHelmholtz2D
import SpecialFunctions: hankelh2
import CompScienceMeshes: normal
import LinearAlgebra: norm, dot

# ------------------------------------------------------------
# Complex-k patch
# ------------------------------------------------------------
function kernelvals(biop::BEAST.HelmholtzOperator2D{T, K}, tgeo, bgeo) where {T, K <: Complex}
    γ = biop.gamma
    k = iszero(real(γ)) ? imag(γ) : -im * γ

    r = tgeo.cart - bgeo.cart
    R = norm(r)

    kr = k * R
    H0 = hankelh2(0, kr)
    H1 = hankelh2(1, kr)

    green     = -im / 4 * H0
    gradgreen =  k * im / 4 * H1 * (r / R)

    txty = dot(normal(tgeo), normal(bgeo))
    return KernelValsHelmholtz2D(γ, r, R, green, gradgreen, txty)
end

# ------------------------------------------------------------
# Block tree
# ------------------------------------------------------------
function build_blocktree(X; minvalues=50)
#function build_blocktree(X; minvalues=50)
    #testtree  = KMeansTree(X.pos, 2; minvalues=minvalues)
    #trialtree = KMeansTree(X.pos, 2; minvalues=minvalues)
    testtree  = BoundingBallTree(X.pos; minvalues=minvalues)
    trialtree = BoundingBallTree(X.pos; minvalues=minvalues)    
    
    # niedzialaja
    #testtree  = SimpleHybridTree(X.pos, 2.0; minvalues=minvalues)
    #trialtree = SimpleHybridTree(X.pos, 2.0; minvalues=minvalues)
    #testtree  = TwoNTree(X.pos, 2.0; minvalues=minvalues)
    #trialtree = TwoNTree(X.pos, 2.0; minvalues=minvalues)
    #testtree  = BoundingBallTree(X.pos, 2; minvalues=minvalues)
    #trialtree = BoundingBallTree(X.pos, 2; minvalues=minvalues)
    println("BUILD_BLOCKTREE = KMEANS VERSION")
    return BlockTree(testtree, trialtree)
end

# ------------------------------------------------------------
# ACA/HMatrix assembly for SINGLE-LAYER
# ------------------------------------------------------------
#function hassemble_with_tree(op, X, tree;
#                             nearquadstrat = BEAST.DoubleNumSauterQstrat(7,7,0,4,9,9),
#                             farquadstrat  = BEAST.DoubleNumQStrat(3,4),
#                             #nearquadstrat = BEAST.DoubleNumSauterQstrat(30,30,0,4,25,25),
#                             #farquadstrat  = BEAST.DoubleNumQStrat(25,30),
#                             scheduler     = DynamicScheduler())
#    return AdaptiveCrossApproximation.HMatrix(
#        op, X, X, tree;
#        nearquadstrat = nearquadstrat,
#        farquadstrat  = farquadstrat,
#        scheduler     = scheduler
#    )
#end

function hassemble_with_tree(op, X, tree;
    η = 0.25,
    nearquadstrat = BEAST.DoubleNumSauterQstrat(30,30,0,4,25,25),
    farquadstrat  = BEAST.DoubleNumQStrat(25,30),
    #nearquadstrat = BEAST.DoubleNumSauterQstrat(7,7,0,4,9,9),
    #farquadstrat  = BEAST.DoubleNumQStrat(3,4),
    #scheduler     = DynamicScheduler()
    #scheduler     = StaticScheduler()
    scheduler     = SerialScheduler() 
    )
return AdaptiveCrossApproximation.HMatrix(
op, X, X, tree;
#isnear        = AdaptiveCrossApproximation.isnear(; η = η),
isnear = AdaptiveCrossApproximation.isnear(; η = η),
nearquadstrat = nearquadstrat,
farquadstrat  = farquadstrat,
scheduler     = scheduler
)
end

# ------------------------------------------------------------
# Context for EFIE
# ------------------------------------------------------------
mutable struct CircleEFIECtx
    X
    tree
    Scache::Dict{ComplexF64, Any}
end

function CircleEFIECtx(X; minvalues=50)
    tree = build_blocktree(X; minvalues=minvalues)
    return CircleEFIECtx(X, tree, Dict{ComplexF64,Any}())
end

# ------------------------------------------------------------
# Cache ACA single-layer operator S(k)
# ------------------------------------------------------------
function getS!(ctx::CircleEFIECtx, k::ComplexF64)
    get!(ctx.Scache, k) do
        𝒮 = Helmholtz2D.singlelayer(; wavenumber = k)
        @time hassemble_with_tree(𝒮, ctx.X, ctx.tree)
    end
end

# ------------------------------------------------------------
# Dense single-layer operator S(k)
# ------------------------------------------------------------
function build_dense_S(X, k::ComplexF64; quadstrat)
    𝒮 = Helmholtz2D.singlelayer(; wavenumber = k)
    return assemble(𝒮, X, X; quadstrat = quadstrat)
end

# ------------------------------------------------------------
# Single-shot comparison:
#   ||(Sfmm - Sd_hi)v|| / ||Sd_hi v||
#   ||(Sd_lo - Sd_hi)v|| / ||Sd_hi v||
# ------------------------------------------------------------
function singleshot_efie(ctx, X, khat::ComplexF64;
        q_lo = BEAST.DoubleNumSauterQstrat(3,4,0,4,5,6),
        q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
        #q_lo = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
        #q_hi = BEAST.DoubleNumSauterQstrat(30,30,0,4,25,25),
        nvec::Int = 5,
        seed::Int = 1)

    Sfmm = getS!(ctx, khat)

    Sd_lo = build_dense_S(X, khat; quadstrat = q_lo)
    Sd_hi = build_dense_S(X, khat; quadstrat = q_hi)

    D = size(Sd_hi, 2)
    Random.seed!(seed)

    ef = zeros(Float64, nvec)
    ed = zeros(Float64, nvec)

    for j in 1:nvec
        v = randn(ComplexF64, D)
        hiv = Sd_hi * v
        ef[j] = norm(Sfmm * v - Sd_hi * v) / (norm(hiv) + eps())
        ed[j] = norm(Sd_lo * v - Sd_hi * v) / (norm(hiv) + eps())
    end

    println("mean(fmm vs dense_hi) = ", mean(ef), "   max = ", maximum(ef))
    println("mean(dense_lo vs hi)  = ", mean(ed), "   max = ", maximum(ed))

    return ef, ed
end

# ------------------------------------------------------------
# Random operator mismatch test
# ------------------------------------------------------------
function op_mismatch_random_efie(ctx, X; k::ComplexF64, ntests::Int=5, seed::Int=1,
                                 quad_dense = BEAST.DoubleNumSauterQstrat(8,9,0,4,9,10))
    rng = MersenneTwister(seed)

    Sdense = build_dense_S(X, k; quadstrat = quad_dense)
    Sfmm   = getS!(ctx, k)

    D = size(Sdense, 2)
    errs = zeros(Float64, ntests)

    for t in 1:ntests
        v = randn(rng, ComplexF64, D)
        Tv_dense = Sdense * v
        Tv_fmm   = Sfmm * v
        errs[t]  = norm(Tv_fmm - Tv_dense) / max(norm(Tv_dense), 1e-300)
        @printf("test %d: err = %.3e\n", t, errs[t])
    end

    @printf("err mean = %.3e, max = %.3e  (k = %.10f%+.10fi)\n",
            mean(errs), maximum(errs), real(k), imag(k))

    return errs
end




# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
a       = 1.0
p_order = 2
h       = 2π * a / (2*246)

circle = CompScienceMeshes.meshcircle(a, h)
X = lagrangecx(circle, order = p_order)

#test
#using H2Trees
#testtree = BoundingBallTree(X.pos; minvalues=20)
#println(testtree)
#bt = BlockTree(testtree, testtree)
#println(bt)

#bt = build_blocktree(X)
#println(bt)

#testtree = KMeansTree(X.pos, 2; minvalues=20)
#@show typeof(testtree)
#@show typeof(data(testtree, root(testtree)))
#@which KMeansTree(X.pos, 2; minvalues=20)


ctx = CircleEFIECtx(X; minvalues=100)

khat = ComplexF64(2.404825557695773 + 0im)

println("=== single-shot EFIE/TM ===")
ef, ed = singleshot_efie(ctx, X, khat; nvec=7, seed=1)

println("\n=== operator mismatch EFIE/TM ===")
errs = op_mismatch_random_efie(ctx, X; k=khat, ntests=5, seed=1)

# test kolumny

function column_test_efie(ctx, X; k::ComplexF64,
                          cols = nothing,
                          quad_dense = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15))

    Sdense = build_dense_S(X, k; quadstrat=quad_dense)
    Sfmm   = getS!(ctx, k)

    n = size(Sdense, 2)

    if cols === nothing
        cols = [1, div(n,2), n]
    end

    errs = Float64[]

    println("\n=== column test EFIE/TM ===")
    @printf("k = %.10f%+.10fi\n", real(k), imag(k))

    for i in cols
        e = zeros(ComplexF64, n)
        e[i] = 1 + 0im

        col_dense = Sdense * e
        col_fmm   = Sfmm * e

        err = norm(col_fmm - col_dense) / max(norm(col_dense), 1e-300)
        push!(errs, err)

        @printf("col %4d: relerr = %.3e   ||dense col|| = %.3e   ||fmm col|| = %.3e\n",
                i, err, norm(col_dense), norm(col_fmm))
    end

    @printf("column err mean = %.3e, max = %.3e\n", mean(errs), maximum(errs))
    return errs
end

function inspect_one_column_efie(ctx, X; k::ComplexF64, col::Int,
                                 quad_dense = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
                                 topn::Int = 20)

    Sdense = build_dense_S(X, k; quadstrat=quad_dense)
    Sfmm   = getS!(ctx, k)

    n = size(Sdense, 2)
    e = zeros(ComplexF64, n)
    e[col] = 1 + 0im

    cd = Sdense * e
    cf = Sfmm * e
    diff = cf - cd

    relerr = norm(diff) / max(norm(cd), 1e-300)
    @printf("\n=== inspect one column EFIE/TM ===\n")
    @printf("col = %d, relerr = %.3e\n", col, relerr)

    idx = sortperm(abs.(diff), rev=true)[1:min(topn, n)]

    println("largest entrywise differences:")
    for j in idx
        @printf("row %4d: |diff|=%.3e   dense=%.3e%+.3ei   fmm=%.3e%+.3ei\n",
                j,
                abs(diff[j]),
                real(cd[j]), imag(cd[j]),
                real(cf[j]), imag(cf[j]))
    end

    return cd, cf, diff, idx
end


n = size(build_dense_S(X, khat; quadstrat=BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15)), 2)
colerrs = column_test_efie(ctx, X; k=khat, cols=[1, 10, div(n,4), div(n,2), n-10, n])

cd, cf, diff, idx = inspect_one_column_efie(ctx, X; k=khat, col=1)

function linearity_test_efie(ctx, X; k::ComplexF64, i::Int=1, j::Int=10)
    Sfmm = getS!(ctx, k)
    n = length(X)

    e1 = zeros(ComplexF64, n); e1[i] = 1
    e2 = zeros(ComplexF64, n); e2[j] = 1

    lhs = Sfmm * (e1 + e2)
    rhs = Sfmm * e1 + Sfmm * e2

    err = norm(lhs - rhs) / max(norm(rhs), 1e-300)

    @printf("linearity test: i=%d, j=%d, relerr = %.3e\n", i, j, err)
    return err
end

err_lin = linearity_test_efie(ctx, X; k=khat, i=1, j=10)

function column_sweep_efie(ctx, X; k::ComplexF64,
                           step::Int = 100,
                           quad_dense = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15))

    Sdense = build_dense_S(X, k; quadstrat=quad_dense)
    Sfmm   = getS!(ctx, k)

    n = size(Sdense, 2)
    cols = collect(1:step:n)
    errs = Float64[]

    println("\n=== column sweep EFIE/TM ===")
    @printf("k = %.10f%+.10fi\n", real(k), imag(k))

    for i in cols
        e = zeros(ComplexF64, n)
        e[i] = 1 + 0im

        cd = Sdense * e
        cf = Sfmm * e

        err = norm(cf - cd) / max(norm(cd), 1e-300)
        push!(errs, err)

        @printf("col %4d: relerr = %.3e\n", i, err)
    end

    @printf("column sweep mean = %.3e, max = %.3e\n", mean(errs), maximum(errs))
    return cols, errs
end
cols, errs = column_sweep_efie(ctx, X; k=khat, step=100)

cols1, errs1 = inspect_one_column_efie(ctx, X; k=khat, col=901)

cols2, errs2 = inspect_one_column_efie(ctx, X; k=khat, col=401)







# dodatkowe
#testtree = TwoNTree(X.pos, 2.0; minvalues=20)
#trialtree = TwoNTree(X.pos, 2.0; minvalues=20)
#
#bt = BlockTree(testtree, trialtree)
#
#@show testtree
#@show trialtree
#@show bt
#
#@show length(testtree)
#@show length(trialtree)
#
#display(testtree)
#display(trialtree)
#display(bt)
#
#testtree  = KMeansTree(X.pos, 2; minvalues=20)
#trialtree = KMeansTree(X.pos, 2; minvalues=20)
#
#@show testtree
#@show trialtree
#
#bt = BlockTree(testtree, trialtree)
#@show bt
#
#θ = atan.(getindex.(X.pos,2), getindex.(X.pos,1))
#θu = copy(θ)
#
#for i in 2:length(θu)
#    while θu[i] < θu[i-1] - π
#        θu[i] += 2π
#    end
#end
#
#issorted(θu)