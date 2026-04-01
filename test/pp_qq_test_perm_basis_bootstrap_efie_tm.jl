# qq_test_perm_basis_bootstrap_efie_tm.jl
println("RUNNING FILE: qq_test_perm_basis_bootstrap_efie_tm.jl")

using LinearAlgebra, Random, Statistics, Printf
using BEAST, CompScienceMeshes
using AdaptiveCrossApproximation
using OhMyThreads
using Krylov
using H2Trees
using SpecialFunctions

import BEAST: kernelvals, KernelValsHelmholtz2D
import SpecialFunctions: hankelh2
import CompScienceMeshes: normal
import LinearAlgebra: norm, dot
import AdaptiveCrossApproximation: permute

# ============================================================
# Permute adapter for BEAST basis
# ============================================================

function permute(space::BEAST.LagrangeBasis, p::AbstractVector{<:Integer})
    pI = collect(Int, p)
    mesh = getfield(space, 1)
    fns  = getfield(space, 2)
    pos  = getfield(space, 3)
    return typeof(space)(mesh, fns[pI], pos[pI])
end

# ============================================================
# Complex-k patch for Helmholtz 2D kernel
# ============================================================

function kernelvals(biop::BEAST.HelmholtzOperator2D{T,K}, tgeo, bgeo) where {T,K<:Complex}
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

# ============================================================
# Contexts
# ============================================================

struct CircleEFIECtx
    X
    tree
end

struct CircleEFIEPermCtx
    Xorig
    Xtest_p
    Xtrial_p
    tree
end

function build_blocktree(X; minvalues=100)
    testtree  = BoundingBallTree(X.pos; minvalues=minvalues)
    trialtree = BoundingBallTree(X.pos; minvalues=minvalues)
    return BlockTree(testtree, trialtree)
end

function make_ctx_circle(; radius::Float64=1.0, Nel::Int=512, minvalues::Int=100)
    h = 2π * radius / Nel
    Γ = CompScienceMeshes.meshcircle(radius, h)
    X = BEAST.lagrangecx(Γ, order=2)
    tree = build_blocktree(X; minvalues=minvalues)
    return CircleEFIECtx(X, tree)
end

"""
Bootstrap path:
- start from original basis
- make fresh copies
- call HMatrix(...; perm=true) ONCE
- let ACA reorder the copies in-place
- afterward use those permuted copies with perm=false
"""
function make_ctx_circle_perm_bootstrap(; radius::Float64=1.0,
                                         Nel::Int=512,
                                         minvalues::Int=100,
                                         kref::ComplexF64 = ComplexF64(2.404825557695773, 0.0),
                                         tol::Float64 = 1e-4,
                                         η::Float64 = 0.25,
                                         scheduler = SerialScheduler(),
                                         nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                                         farquadstrat  = BEAST.DoubleNumQStrat(5,4))

    ctx = make_ctx_circle(; radius=radius, Nel=Nel, minvalues=minvalues)

    Xtest_p  = deepcopy(ctx.X)
    Xtrial_p = deepcopy(ctx.X)

    Sop = Helmholtz2D.singlelayer(; wavenumber = kref)

    _ = AdaptiveCrossApproximation.HMatrix(
        Sop, Xtest_p, Xtrial_p, ctx.tree;
        compressor    = ACA(; convergence = FNormEstimator(tol)),
        isnear        = AdaptiveCrossApproximation.isnear(; η = η),
        nearquadstrat = nearquadstrat,
        farquadstrat  = farquadstrat,
        scheduler     = scheduler,
        perm          = true,
    )

    changed_test  = any(Xtest_p.pos[i]  != ctx.X.pos[i] for i in eachindex(ctx.X.pos))
    changed_trial = any(Xtrial_p.pos[i] != ctx.X.pos[i] for i in eachindex(ctx.X.pos))

    println("perm bootstrap changed test order  = ", changed_test)
    println("perm bootstrap changed trial order = ", changed_trial)

    return CircleEFIEPermCtx(ctx.X, Xtest_p, Xtrial_p, ctx.tree)
end

# ============================================================
# Builders: original basis
# ============================================================

function build_dense_S(ctx::CircleEFIECtx, k::ComplexF64;
                       quadstrat = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15))
    Sop = Helmholtz2D.singlelayer(; wavenumber = k)
    S   = assemble(Sop, ctx.X, ctx.X; quadstrat = quadstrat)
    return Matrix(S)
end

function build_fmm_S(ctx::CircleEFIECtx, k::ComplexF64;
                     tol::Float64 = 1e-4,
                     η::Float64   = 0.25,
                     scheduler     = SerialScheduler(),
                     nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                     farquadstrat  = BEAST.DoubleNumQStrat(5,4))
    Sop = Helmholtz2D.singlelayer(; wavenumber = k)
    Sfmm = AdaptiveCrossApproximation.HMatrix(
        Sop, ctx.X, ctx.X, ctx.tree;
        compressor    = ACA(; convergence = FNormEstimator(tol)),
        isnear        = AdaptiveCrossApproximation.isnear(; η = η),
        nearquadstrat = nearquadstrat,
        farquadstrat  = farquadstrat,
        scheduler     = scheduler,
        perm          = false,
    )
    return Sfmm
end

# ============================================================
# Builders: permuted-basis bootstrap path
# ============================================================

function build_dense_S(ctx::CircleEFIEPermCtx, k::ComplexF64;
                       quadstrat = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15))
    Sop = Helmholtz2D.singlelayer(; wavenumber = k)
    S   = assemble(Sop, ctx.Xtest_p, ctx.Xtrial_p; quadstrat = quadstrat)
    return Matrix(S)
end

function build_fmm_S(ctx::CircleEFIEPermCtx, k::ComplexF64;
                     tol::Float64 = 1e-4,
                     η::Float64   = 0.25,
                     scheduler     = SerialScheduler(),
                     nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                     farquadstrat  = BEAST.DoubleNumQStrat(5,4))
    Sop = Helmholtz2D.singlelayer(; wavenumber = k)
    Sfmm = AdaptiveCrossApproximation.HMatrix(
        Sop, ctx.Xtest_p, ctx.Xtrial_p, ctx.tree;
        compressor    = ACA(; convergence = FNormEstimator(tol)),
        isnear        = AdaptiveCrossApproximation.isnear(; η = η),
        nearquadstrat = nearquadstrat,
        farquadstrat  = farquadstrat,
        scheduler     = scheduler,
        perm          = false,
    )
    return Sfmm
end

# ============================================================
# Helpers
# ============================================================

function relerrs_apply(A, B; nvec::Int=7, seed::Int=1)
    Random.seed!(seed)
    n = size(B, 2)
    errs = Vector{Float64}(undef, nvec)

    for i in 1:nvec
        v = randn(ComplexF64, n)
        Av = A * v
        Bv = B * v
        errs[i] = norm(Av - Bv) / norm(Bv)
    end
    return errs
end

# ============================================================
# Operator comparison
# ============================================================

function singleshot_compare(ctx, k::ComplexF64;
                            q_lo = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                            q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
                            tol = 1e-4,
                            η   = 0.25,
                            scheduler = SerialScheduler(),
                            nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                            farquadstrat  = BEAST.DoubleNumQStrat(5,4),
                            nvec::Int=7,
                            seed::Int=1,
                            label::String="EFIE")

    @printf("\n=== singleshot_compare (%s) ===\n", label)
    @printf("k = %+.12f%+.12fi\n", real(k), imag(k))
    @printf("scheduler = %s, ACA tol = %.1e, η = %.3f\n",
            string(typeof(scheduler)), tol, η)

    Td_lo = build_dense_S(ctx, k; quadstrat=q_lo)
    Td_hi = build_dense_S(ctx, k; quadstrat=q_hi)

    TF = build_fmm_S(ctx, k;
        tol=tol, η=η, scheduler=scheduler,
        nearquadstrat=nearquadstrat, farquadstrat=farquadstrat)

    e_fmm   = relerrs_apply(TF,    Td_hi; nvec=nvec, seed=seed)
    e_dense = relerrs_apply(Td_lo, Td_hi; nvec=nvec, seed=seed)

    @printf("mean(fmm   vs dense_hi) = %.12g   max = %.12g\n",
            mean(e_fmm), maximum(e_fmm))
    @printf("mean(dense_lo vs hi)    = %.12g   max = %.12g\n",
            mean(e_dense), maximum(e_dense))

    return e_fmm, e_dense
end

# ============================================================
# Solve comparison on permuted-basis bootstrap path
# ============================================================

function compare_solve(ctx::CircleEFIEPermCtx, k::ComplexF64;
                       tol::Float64 = 1e-4,
                       η::Float64   = 0.25,
                       gmres_rtol::Float64 = 1e-10,
                       gmres_itmax::Int = 1500,
                       scheduler = SerialScheduler(),
                       q_dense = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
                       nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                       farquadstrat  = BEAST.DoubleNumQStrat(5,4),
                       seed::Int = 123)

    Random.seed!(seed)

    Sdense = build_dense_S(ctx, k; quadstrat=q_dense)
    Sfmm   = build_fmm_S(ctx, k;
                         tol=tol, η=η, scheduler=scheduler,
                         nearquadstrat=nearquadstrat,
                         farquadstrat=farquadstrat)

    n = size(Sdense, 1)
    b = randn(ComplexF64, n)

    xd = Sdense \ b
    xf, st = Krylov.gmres(Sfmm, b; rtol=gmres_rtol, itmax=gmres_itmax, verbose=0)

    rd_dense = norm(Sdense * xd - b) / norm(b)
    rf_fmmop = norm(Sfmm   * xf - b) / norm(b)
    rf_dense = norm(Sdense * xf - b) / norm(b)

    @printf("\n=== compare_solve (bootstrap permuted basis) ===\n")
    @printf("k = %+.12f%+.12fi\n", real(k), imag(k))
    @printf("GMRES solved = %s, status = %s, niter = %d\n",
            string(st.solved), string(st.status), st.niter)
    @printf("dense solve residual in dense op = %.3e\n", rd_dense)
    @printf("ACA   solve residual in ACA op   = %.3e\n", rf_fmmop)
    @printf("ACA   solve residual in dense op = %.3e\n", rf_dense)
    @printf("relative difference xf vs xd     = %.3e\n", norm(xf - xd) / norm(xd))

    return xd, xf, st
end

# ============================================================
# Main
# ============================================================

function run_all()
    k0 = ComplexF64(2.404825557695773, 0.0)

    println("\n================ ORIGINAL BASIS =================\n")
    ctx0 = make_ctx_circle(; radius=1.0, Nel=512, minvalues=100)

    singleshot_compare(ctx0, k0;
        q_lo = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
        tol = 1e-4,
        η   = 0.25,
        scheduler = SerialScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        farquadstrat  = BEAST.DoubleNumQStrat(5,4),
        nvec = 5,
        seed = 1,
        label = "ACA original basis"
    )

    println("\n================ PERMUTED-BASIS BOOTSTRAP =================\n")
    ctxp = make_ctx_circle_perm_bootstrap(; radius=1.0, Nel=512, minvalues=100,
        kref = k0,
        tol = 1e-4,
        η   = 0.25,
        scheduler = SerialScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        farquadstrat  = BEAST.DoubleNumQStrat(5,4)
    )

    singleshot_compare(ctxp, k0;
        q_lo = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
        tol = 1e-4,
        η   = 0.25,
        scheduler = SerialScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        farquadstrat  = BEAST.DoubleNumQStrat(5,4),
        nvec = 5,
        seed = 1,
        label = "ACA bootstrap permuted basis"
    )

    compare_solve(ctxp, k0 + 0.03im;
        tol = 1e-4,
        η   = 0.25,
        gmres_rtol = 1e-10,
        gmres_itmax = 1500,
        scheduler = SerialScheduler(),
        q_dense = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        farquadstrat  = BEAST.DoubleNumQStrat(5,4),
        seed = 123
    )

    return nothing
end

run_all()