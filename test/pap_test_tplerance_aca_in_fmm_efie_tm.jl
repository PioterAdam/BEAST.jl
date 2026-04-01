# pap_test_aca_sweep_efie.jl
#
# ACA/FMM diagnostic sweep for 2D Helmholtz single-layer operator (EFIE)
# on a circle boundary.
#
# Goal:
#   compare ACA/HMatrix operator against dense reference
#   before using Beyn / contour solvers.
#
println("RUNNING FILE: pap_test_tplerance_aca_in_fmm_efie_tm.jl")

using LinearAlgebra, SparseArrays, Random, Statistics, Printf
using BEAST, CompScienceMeshes
using AdaptiveCrossApproximation
using OhMyThreads
using Krylov
using ParallelKMeans
using H2Trees
using SpecialFunctions
using LinearMaps

import BEAST: kernelvals, KernelValsHelmholtz2D
import SpecialFunctions: hankelh2
import CompScienceMeshes: normal
import LinearAlgebra: norm, dot

import AdaptiveCrossApproximation: permute

function permute(space::BEAST.LagrangeBasis, p::AbstractVector{<:Integer})
    pI = collect(Int, p)

    mesh = getfield(space, 1)
    fns  = getfield(space, 2)
    pos  = getfield(space, 3)

    return typeof(space)(mesh, fns[pI], pos[pI])
end

#import AdaptiveCrossApproximation: PermutedHMatrix

#AdaptiveCrossApproximation.PermutedHMatrix(
#    perms::Tuple,
#    H::AdaptiveCrossApproximation.HMatrix,
#) = AdaptiveCrossApproximation.PermutedHMatrix(H, perms)

# ============================================================
# Complex-k patch for Helmholtz 2D kernel
# ============================================================

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

# ============================================================
# Context: geometry + tree + cache
# ============================================================

struct CircleEFIECtx
    X
    tree
    cache::Dict{ComplexF64,Any}
end

function build_blocktree(X; minvalues=50)
    #testtree  = KMeansTree(X.pos, 2; minvalues=minvalues)
    #trialtree = KMeansTree(X.pos, 2; minvalues=minvalues)
    testtree  = BoundingBallTree(X.pos; minvalues=minvalues)
    trialtree = BoundingBallTree(X.pos; minvalues=minvalues) 
    return BlockTree(testtree, trialtree)
end

"""
Create circle mesh with Nel segments and scalar basis X.
"""
function make_ctx_circle(; radius::Float64=1.0, Nel::Int=512, minvalues::Int=100)
    h = 2π * radius / Nel
    Γ = CompScienceMeshes.meshcircle(radius, h)

    # keep the older/simple basis choice from your harness
    #X = BEAST.lagrangec0d1(Γ)
    p_order = 2
    X = lagrangecx(Γ, order = p_order)

    tree = build_blocktree(X; minvalues=minvalues)

    return CircleEFIECtx(X, tree, Dict{ComplexF64,Any}())
end

# ============================================================
# Assemblers for EFIE: T(k) = S(k)
# ============================================================

"""
Dense operator S(k), assembled with a chosen quadrature rule.
"""
function build_dense_T(ctx::CircleEFIECtx, k::ComplexF64; quadstrat)
    Sop = Helmholtz2D.singlelayer(; wavenumber = k)
    S   = assemble(Sop, ctx.X, ctx.X; quadstrat = quadstrat)
    return S
end

"""
ACA/HMatrix approximation of the same single-layer operator S(k).
"""
function build_fmm_T(ctx::CircleEFIECtx, k::ComplexF64;
                     tol::Float64 = 1e-4,
                     η::Float64   = 1.0,
                     #scheduler    = DynamicScheduler(),
                     scheduler    = SerialScheduler(),
                     nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                     farquadstrat  = BEAST.DoubleNumQStrat(4,4))

    Sop = Helmholtz2D.singlelayer(; wavenumber = k)

    Sfmm = AdaptiveCrossApproximation.HMatrix(
        Sop, ctx.X, ctx.X, ctx.tree;
        #compressor    = ACA(; tol = tol),
        compressor = ACA(; convergence = FNormEstimator(tol)),
        isnear        = AdaptiveCrossApproximation.isnear(; η = η),
        nearquadstrat = nearquadstrat,
        farquadstrat  = farquadstrat,
        scheduler     = scheduler,
        perm          = false,
    )

    return Sfmm
end

# ============================================================
# Robust relative error on random vectors:
#   ||(A v - B v)|| / ||B v||
# ============================================================

function relerrs_apply(A, B; nvec::Int=7, seed::Int=1)
    Random.seed!(seed)
    n = size(B, 1)
    errs = Vector{Float64}(undef, nvec)

    for i in 1:nvec
        v = randn(ComplexF64, n)
        Av = A * v
        Bv = B * v
        errs[i] = norm(Av - Bv) / norm(Bv)
    end

    return errs
end

function relerrs_apply_(A, B; nvec::Int=7, seed::Int=1)
    Random.seed!(seed)
    n = size(B, 1)
    errs = Vector{Float64}(undef, nvec)

    for i in 1:nvec
        v = randn(ComplexF64, n)
        Av = A * v
        Bv = B * v

        if !all(isfinite, Av)
            @printf("relerrs_apply: non-finite Av at vec %d\n", i)
        end
        if !all(isfinite, Bv)
            @printf("relerrs_apply: non-finite Bv at vec %d\n", i)
        end

        nB = norm(Bv)
        nD = norm(Av - Bv)

        @printf("vec %d: norm(Av-Bv)=%.6e  norm(Bv)=%.6e\n", i, nD, nB)

        errs[i] = nD / nB
    end

    return errs
end

# ============================================================
# Single-shot comparison
# ============================================================

function singleshot_compare(ctx::CircleEFIECtx, k::ComplexF64;
                            q_lo = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                            q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
                            tol = 1e-4,
                            η   = 0.25,
                            scheduler = SerialScheduler(), # albo StaticScheduler()
                            nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                            farquadstrat  = BEAST.DoubleNumQStrat(4,4),
                            nvec::Int=7,
                            seed::Int=1)

    @printf("\n=== singleshot_compare (EFIE / single-layer) ===\n")
    @printf("k = %+.12f%+.12fi\n", real(k), imag(k))
    @printf("scheduler = %s, ACA tol = %.1e, η = %.3f\n",
            string(typeof(scheduler)), tol, η)

    Td_lo = build_dense_T(ctx, k; quadstrat = q_lo)
    Td_hi = build_dense_T(ctx, k; quadstrat = q_hi)

    TF = build_fmm_T(ctx, k;
        tol           = tol,
        η             = η,
        scheduler     = scheduler,
        nearquadstrat = nearquadstrat,
        farquadstrat  = farquadstrat,
    )

    e_fmm   = relerrs_apply(TF,    Td_hi; nvec=nvec, seed=seed)
    e_dense = relerrs_apply(Td_lo, Td_hi; nvec=nvec, seed=seed)

    @printf("mean(fmm   vs dense_hi) = %.12g   max = %.12g\n",
            mean(e_fmm), maximum(e_fmm))
    @printf("mean(dense_lo vs hi)    = %.12g   max = %.12g\n",
            mean(e_dense), maximum(e_dense))

    return e_fmm, e_dense
end

function singleshot_compare_(ctx::CircleEFIECtx, k::ComplexF64;
                            q_lo = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                            q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
                            tol = 1e-4,
                            η   = 0.25,
                            scheduler = SerialScheduler(),#StaticScheduler(),
                            nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                            farquadstrat  = BEAST.DoubleNumQStrat(4,4),
                            nvec::Int=7,
                            seed::Int=1)

    @printf("\n=== singleshot_compare (EFIE / single-layer) ===\n")
    @printf("k = %+.12f%+.12fi\n", real(k), imag(k))
    @printf("scheduler = %s, ACA tol = %.1e, η = %.3f\n", string(typeof(scheduler)), tol, η)

    Td_lo = build_dense_T(ctx, k; quadstrat = q_lo)
    Td_hi = build_dense_T(ctx, k; quadstrat = q_hi)

    TF = build_fmm_T(ctx, k;
        tol           = tol,
        η             = η,
        scheduler     = scheduler,
        nearquadstrat = nearquadstrat,
        farquadstrat  = farquadstrat,
    )

    vtest = randn(ComplexF64, size(Td_hi, 2))

    @show @which(TF * vtest)
    @show @which(mul!(zeros(ComplexF64, size(Td_hi,1)), TF, vtest))
    @show @which(LinearMaps._unsafe_mul!(zeros(ComplexF64, size(Td_hi,1)), TF, vtest))

    @show typeof(TF)
    if hasproperty(TF, :permutation)
        @show length(TF.permutation[1]), length(TF.permutation[2])
        @show TF.permutation[1][1:10]
        @show TF.permutation[2][1:10]
    end

    if hasproperty(TF, :permutation)
        p1, p2 = TF.permutation
        @show p1 == collect(1:length(p1))
        @show p2 == collect(1:length(p2))
    end

        # ------------------------------------------------------------
    # Internal consistency test:
    # compare internal HMatrix TF.M against dense operator
    # permuted by the same row/col permutations
    # ------------------------------------------------------------
    if hasproperty(TF, :permutation) && hasproperty(TF, :M)
        ptest, ptrial = TF.permutation

        Td_perm = Td_hi[ptest, ptrial]

        Random.seed!(seed)
        n_int = size(Td_perm, 2)
        e_int = Vector{Float64}(undef, nvec)

        for i in 1:nvec
            xp = randn(ComplexF64, n_int)
            Hv = TF.M * xp
            Dv = Td_perm * xp
            e_int[i] = norm(Hv - Dv) / norm(Dv)
        end

        @printf("mean(TF.M vs Td_hi[ptest,ptrial]) = %.12g   max = %.12g\n",
                mean(e_int), maximum(e_int))
            # probe test: internal HMatrix vs dense-permuted on one fixed vector
        ptest, ptrial = TF.permutation
        vprobe = randn(ComplexF64, size(Td_hi, 2))

        xp_probe = vprobe[ptrial]
        yp_h = TF.M * xp_probe
        yp_d = Td_hi[ptest, ptrial] * xp_probe

        @printf("probe: TF.M*xp vs denseperm*xp = %.12g\n",
                norm(yp_h - yp_d) / norm(yp_d))

        y_tf = TF * vprobe
        y_manual = zeros(ComplexF64, size(Td_hi, 1))
        y_manual[ptest] .= yp_h

        @printf("probe: TF*v vs manual-from-TF.M = %.12g\n",
                norm(y_tf - y_manual) / norm(y_manual))
    end

        # ------------------------------------------------------------
    # Wrapper consistency test:
    # compare TF * v against manual gather/scatter built from TF.M
    # ------------------------------------------------------------
    if hasproperty(TF, :permutation) && hasproperty(TF, :M)
        ptest, ptrial = TF.permutation

        Random.seed!(seed)
        e_wrap = Vector{Float64}(undef, nvec)

        for i in 1:nvec
            v = randn(ComplexF64, size(Td_hi, 2))

            y_tf = TF * v

            xp = v[ptrial]
            yp = TF.M * xp

            y_manual = zeros(ComplexF64, size(Td_hi, 1))
            y_manual[ptest] .= yp

            e_wrap[i] = norm(y_tf - y_manual) / norm(y_manual)
        end

        @printf("mean(TF vs manual wrapper) = %.12g   max = %.12g\n",
                mean(e_wrap), maximum(e_wrap))
    end

        # ------------------------------------------------------------
    # Dense permutation-wrapper consistency test
    # ------------------------------------------------------------
    if hasproperty(TF, :permutation)
        ptest, ptrial = TF.permutation

        Random.seed!(seed)
        e_densewrap = Vector{Float64}(undef, nvec)

        for i in 1:nvec
            v = randn(ComplexF64, size(Td_hi, 2))

            yp = Td_hi[ptest, ptrial] * v[ptrial]

            y_manual_dense = zeros(ComplexF64, size(Td_hi, 1))
            y_manual_dense[ptest] .= yp

            y_true = Td_hi * v

            e_densewrap[i] = norm(y_manual_dense - y_true) / norm(y_true)
        end

        @printf("mean(dense perm-wrapper vs dense_hi) = %.12g   max = %.12g\n",
                mean(e_densewrap), maximum(e_densewrap))
    end

    e_fmm   = relerrs_apply(TF,    Td_hi; nvec=nvec, seed=seed)
    e_dense = relerrs_apply(Td_lo, Td_hi; nvec=nvec, seed=seed)

    @printf("mean(fmm   vs dense_hi) = %.12g   max = %.12g\n", mean(e_fmm), maximum(e_fmm))
    @printf("mean(dense_lo vs hi)    = %.12g   max = %.12g\n", mean(e_dense), maximum(e_dense))

    return e_fmm, e_dense
end

# ============================================================
# Sweep ACA tol / scheduler / eta
# ============================================================

function sweep_aca(ctx::CircleEFIECtx, k::ComplexF64;
                   tols = [1e-3, 1e-4, 1e-5, 1e-6],
                   etas = [1.0],
                   schedulers = [StaticScheduler(), DynamicScheduler()],
                   q_lo = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                   q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
                   nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                   farquadstrat  = BEAST.DoubleNumQStrat(4,4),
                   nvec::Int=7,
                   seed::Int=1)

    @printf("\n============================================\n")
    @printf("ACA sweep for EFIE single-layer (fixed X/tree)\n")
    @printf("k = %+.12f%+.12fi\n", real(k), imag(k))
    @printf("============================================\n")

    Td_lo = build_dense_T(ctx, k; quadstrat = q_lo)
    Td_hi = build_dense_T(ctx, k; quadstrat = q_hi)

    e_dense = relerrs_apply(Td_lo, Td_hi; nvec=nvec, seed=seed)
    @printf("Dense self-check: mean=%.3e  max=%.3e\n", mean(e_dense), maximum(e_dense))

    for sch in schedulers
        @printf("\n=== scheduler = %s ===\n", string(typeof(sch)))
        for η in etas
            @printf("--- η = %.3f ---\n", η)
            #for tol in tols
            #    try
            #        TF = build_fmm_T(ctx, k;
            #            tol           = tol,
            #            η             = η,
            #            scheduler     = sch,
            #            nearquadstrat = nearquadstrat,
            #            farquadstrat  = farquadstrat,
            #        )
            #        e_fmm = relerrs_apply(TF, Td_hi; nvec=nvec, seed=seed)
            #        @printf("tol=%-8.1e  mean(FMM vs dense_hi)=%.3e  max=%.3e\n",
            #                tol, mean(e_fmm), maximum(e_fmm))
            #    catch err
            #        @printf("tol=%-8.1e  ERROR: %s\n", tol, sprint(showerror, err))
            #    end
            #end
            for tol in tols
                try
                    TF = build_fmm_T(ctx, k;
                        tol           = tol,
                        η             = η,
                        scheduler     = sch,
                        nearquadstrat = nearquadstrat,
                        farquadstrat  = farquadstrat,
                    )
            
                    e_fmm = relerrs_apply(TF, Td_hi; nvec=nvec, seed=seed)
            
                    if any(x -> !isfinite(x), e_fmm)
                        @printf("tol=%-8.1e  NON-FINITE result (likely far-assembly / hankelh2 issue)\n", tol)
                    else
                        @printf("tol=%-8.1e  mean(FMM vs dense_hi)=%.3e  max=%.3e\n",
                                tol, mean(e_fmm), maximum(e_fmm))
                    end
            
                catch err
                    msg = sprint(showerror, err)
                    if occursin("AmosException", msg) || occursin("hankelh2", msg)
                        @printf("tol=%-8.1e  ERROR: far-assembly/kernel sampling issue (%s)\n", tol, msg)
                    else
                        @printf("tol=%-8.1e  ERROR: %s\n", tol, msg)
                    end
                end
            end
        end
    end

    return nothing
end

# ============================================================
# Optional solve comparison at one k
# ============================================================

function compare_solve(ctx::CircleEFIECtx, k::ComplexF64;
                       tol::Float64 = 1e-6,
                       η::Float64   = 1.0,
                       gmres_rtol::Float64 = 1e-6,
                       gmres_itmax::Int = 800,
                       scheduler = DynamicScheduler(),
                       q_dense = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
                       nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                       farquadstrat  = BEAST.DoubleNumQStrat(4,4),
                       seed::Int = 123)

    Random.seed!(seed)

    Sdense = build_dense_T(ctx, k; quadstrat=q_dense)
    Sfmm   = build_fmm_T(ctx, k;
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

    @printf("\n=== compare_solve (EFIE / single-layer) ===\n")
    @printf("k = %+.12f%+.12fi\n", real(k), imag(k))
    @printf("GMRES solved = %s, status = %s, niter = %d\n", st.solved, st.status, st.niter)
    @printf("dense solve residual in dense op = %.3e\n", rd_dense)
    @printf("FMM solve residual in FMM op     = %.3e\n", rf_fmmop)
    @printf("FMM solve residual in dense op   = %.3e\n", rf_dense)
    @printf("relative difference xf vs xd     = %.3e\n", norm(xf - xd) / norm(xd))

    return xd, xf, st
end

# ============================================================
# MAIN
# ============================================================

function main1()
    Nel = 512
    ctx = make_ctx_circle(; radius=1.0, Nel=Nel, minvalues=100)

    k0 = ComplexF64(2.404825557695773, 0.0)

    singleshot_compare(ctx, k0;
        q_lo = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
        tol = 1e-4,
        η   = 0.25,
        scheduler = StaticScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        farquadstrat  = BEAST.DoubleNumQStrat(4,4),
        nvec = 7,
        seed = 123
    )
end


function main()
    Nel = 512
    ctx = make_ctx_circle(; radius=1.0, Nel=Nel, minvalues=100)

    # first Dirichlet eigenvalue of disk, useful diagnostic target
    k0 = ComplexF64(2.404825557695773, 0.0)

    # 1) single shot
    singleshot_compare(ctx, k0;
        q_lo = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
        tol = 1e-4,
        η   = 0.25,
        scheduler = StaticScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        farquadstrat  = BEAST.DoubleNumQStrat(4,4),
        nvec = 7,
        seed = 123
    )

    # 2) sweep
    sweep_aca(ctx, k0;
        tols = [1e-4],#[1e-3, 1e-4, 1e-5, 1e-6],
        etas = [0.25], #[0.125 0.25 0.5 0.75 1.0], #, 2.0, 10.0
        schedulers = [StaticScheduler(), DynamicScheduler()],
        q_lo = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        farquadstrat  = BEAST.DoubleNumQStrat(4,4),
        nvec = 7,
        seed = 123
    )

    # 3) optional one-k solve test
    compare_solve(ctx, k0 + 0.03im;
        tol = 1e-6,
        η   = 0.25,
        gmres_rtol = 1e-6,
        gmres_itmax = 800,
        scheduler = DynamicScheduler(),
        q_dense = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        farquadstrat  = BEAST.DoubleNumQStrat(4,4),
        seed = 123
    )
end

ctx = make_ctx_circle(; radius=1.0, Nel=512, minvalues=100)
k0 = ComplexF64(2.404825557695773, 0.0)

singleshot_compare(ctx, k0;
    q_lo = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
    q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
    tol = 1e-4,
    η   = 0.25, #0.25--5
    scheduler = StaticScheduler(),
    nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
    farquadstrat  = BEAST.DoubleNumQStrat(4,4),
    nvec = 7,
    seed = 123
)

main1()


main()

#using BEAST, CompScienceMeshes, AdaptiveCrossApproximation, ParallelKMeans, H2Trees, OhMyThreads

#import AdaptiveCrossApproximation: permute

#function permute(space::BEAST.LagrangeBasis, p::AbstractVector{<:Integer})
#    pI = collect(Int, p)
#
#    mesh = getfield(space, 1)
#    fns  = getfield(space, 2)
#    pos  = getfield(space, 3)
#
#    return typeof(space)(mesh, fns[pI], pos[pI])
#end
#
#Γ = CompScienceMeshes.meshcircle(1.0, 2π/128)
#p_order = 2
#X = lagrangecx(Γ, order = p_order)
#X = BEAST.lagrangec0d1(Γ)

#function build_blocktree(X; minvalues=100)
#    testtree  = KMeansTree(X.pos, 2; minvalues=minvalues)
#    trialtree = KMeansTree(X.pos, 2; minvalues=minvalues)
#    return BlockTree(testtree, trialtree)
##end
#
#tree = build_blocktree(X)
#
#k = ComplexF64(2.404825557695773, 0.0)
#Sop = Helmholtz2D.singlelayer(; wavenumber = k)

#H = AdaptiveCrossApproximation.HMatrix(
#    Sop, X, X, tree;
#    compressor    = ACA(; tol = 1e-4),
#    isnear        = AdaptiveCrossApproximation.isnear(; η = 1.0),
#    nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
#    farquadstrat  = BEAST.DoubleNumQStrat(4,4),
#    scheduler     = StaticScheduler(),
#    perm          = false,
#)