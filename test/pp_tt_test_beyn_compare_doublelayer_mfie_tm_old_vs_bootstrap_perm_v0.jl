println("RUNNING FILE: tt_test_beyn_compare_doublelayer_mfie_tm_old_vs_bootstrap_perm.jl")

using LinearAlgebra, Random, Statistics, Printf
using BEAST, CompScienceMeshes
using AdaptiveCrossApproximation
using OhMyThreads
using Krylov
using H2Trees
using SpecialFunctions
using FunctionZeros
using LinearMaps

import BEAST: kernelvals, KernelValsHelmholtz2D
import SpecialFunctions: hankelh2
import CompScienceMeshes: normal
import LinearAlgebra: norm, dot, mul!
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
# Generic Beyn solver
# ============================================================

function BeynSolve_auto(z0, Rx, Ry, UserBeynFunction, ApplyT, D;
                        L::Int = D,
                        N::Int = 64,
                        Kmax::Int = 3,
                        tol_svd::Float64 = 1e-8,
                        tol_res::Float64 = 1e-8,
                        seed::Int = 1,
                        verbose::Bool = true)

    Random.seed!(seed)
    VHat = randn(ComplexF64, D, L)

    nMom_max = 2 * Kmax
    Ap = [zeros(ComplexF64, D, L) for _ in 1:nMom_max]

    Δθ = 2π / N
    MInvVHat = zeros(ComplexF64, D, L)

    for nn in 0:(N-1)
        θ  = nn * Δθ
        cθ = cos(θ)
        sθ = sin(θ)

        w  = Rx * cθ + im * Ry * sθ
        z  = z0 + w
        dz = (-Rx * sθ + im * Ry * cθ) * Δθ

        UserBeynFunction(z, VHat, MInvVHat)

        wp = one(w)
        @inbounds for p in 1:nMom_max
            Ap[p] .+= wp * dz .* MInvVHat
            wp *= w
        end
    end

    function Beyn_from_moments(K::Int)
        m  = D
        ℓ  = L
        Km = K * m
        Kl = K * ℓ

        B0 = zeros(ComplexF64, Km, Kl)
        B1 = zeros(ComplexF64, Km, Kl)

        for i in 1:K
            rows = (i-1)*m+1 : i*m
            for j in 1:K
                cols = (j-1)*ℓ+1 : j*ℓ
                p0 = i + j - 2
                p1 = i + j - 1
                B0[rows, cols] .= Ap[p0+1]
                B1[rows, cols] .= Ap[p1+1]
            end
        end

        U, σ, W = svd(B0)
        k = count(>(tol_svd), abs.(σ))
        verbose && @printf("K=%d: rank(B0) ≈ %d\n", K, k)

        if k == 0
            return ComplexF64[], zeros(ComplexF64, D, 0), Int[], Float64[]
        end

        U0 = U[:, 1:k]
        σk = σ[1:k]
        W0 = W[:, 1:k]

        Dred = U0' * B1 * W0 * Diagonal(1.0 ./ σk)
        F = eigen(Dred)

        μ = F.values
        S = F.vectors

        Vblock1 = U0[1:m, :]
        Vfull   = Vblock1 * S
        λ       = μ .+ z0

        r = zeros(Float64, length(λ))
        good_idxs = Int[]

        verbose && @printf("  Residuals (K=%d):\n", K)
        for j in eachindex(λ)
            vj = Vfull[:, j]
            Tv = ApplyT(λ[j], vj)
            rj = norm(Tv) / norm(vj)
            r[j] = rj
            verbose && @printf("    res[%2d] = %.3e\n", j, rj)
            if rj <= tol_res
                push!(good_idxs, j)
            end
        end

        return λ, Vfull, good_idxs, r
    end

    λ_best    = ComplexF64[]
    V_best    = zeros(ComplexF64, D, 0)
    good_best = 0
    good_prev = 0
    K_used    = 0

    for K in 1:Kmax
        λ, V, good_idxs, _ = Beyn_from_moments(K)
        ngood = length(good_idxs)

        verbose && @printf("  K=%d: %d good eigenpairs (res ≤ %.1e)\n", K, ngood, tol_res)

        if ngood == 0
            continue
        end

        λg = λ[good_idxs]
        Vg = V[:, good_idxs]

        if good_best == 0 || ngood > good_prev
            good_best = ngood
            λ_best    = λg
            V_best    = Vg
            K_used    = K
        elseif ngood == good_prev
            verbose && @printf("  #good eigenpairs stabilized at %d for K=%d.\n", ngood, K)
            λ_best = λg
            V_best = Vg
            K_used = K
            break
        end

        good_prev = ngood
    end

    if good_best == 0 || isempty(λ_best)
        verbose && println("No acceptable eigenpairs found.")
        return ComplexF64[], zeros(ComplexF64, D, 0), false
    end

    p = sortperm(λ_best, by = x -> (real(x), imag(x)))
    λ_best = λ_best[p]
    V_best = V_best[:, p]

    verbose && @printf("Using K=%d with %d good eigenpairs.\n", K_used, good_best)

    return λ_best, V_best, true
end

function BeynSolve_auto(z0, R, UserBeynFunction, ApplyT, D; kwargs...)
    return BeynSolve_auto(z0, R, R, UserBeynFunction, ApplyT, D; kwargs...)
end

# ============================================================
# Contexts
# ============================================================

struct CircleDLKCtx
    X
    tree
end

struct CircleDLKPermCtx
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
    return CircleDLKCtx(X, tree)
end

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

    Kop = Helmholtz2D.doublelayer(; wavenumber = kref)

    _ = AdaptiveCrossApproximation.HMatrix(
        Kop, Xtest_p, Xtrial_p, ctx.tree;
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

    return CircleDLKPermCtx(ctx.X, Xtest_p, Xtrial_p, ctx.tree)
end

# ============================================================
# Dense builders: D and Iop
# ============================================================

function build_dense_D(ctx::CircleDLKCtx, k::ComplexF64;
                       quadstrat = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15))
    Kop = Helmholtz2D.doublelayer(; wavenumber = k)
    D = assemble(Kop, ctx.X, ctx.X; quadstrat = quadstrat)
    return Matrix(D)
end

function build_dense_D(ctx::CircleDLKPermCtx, k::ComplexF64;
                       quadstrat = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15))
    Kop = Helmholtz2D.doublelayer(; wavenumber = k)
    D = assemble(Kop, ctx.Xtest_p, ctx.Xtrial_p; quadstrat = quadstrat)
    return Matrix(D)
end

function build_dense_Iop(ctx::CircleDLKCtx)
    Iop = assemble(BEAST.Identity(), ctx.X, ctx.X)
    return Matrix(Iop)
end

function build_dense_Iop(ctx::CircleDLKPermCtx)
    Iop = assemble(BEAST.Identity(), ctx.Xtest_p, ctx.Xtrial_p)
    return Matrix(Iop)
end

# ============================================================
# ACA builder: D only
# ============================================================

function build_fmm_D(ctx::CircleDLKCtx, k::ComplexF64;
                     tol::Float64 = 1e-4,
                     η::Float64   = 0.25,
                     scheduler = SerialScheduler(),
                     nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                     farquadstrat  = BEAST.DoubleNumQStrat(5,4))
    Kop = Helmholtz2D.doublelayer(; wavenumber = k)
    return AdaptiveCrossApproximation.HMatrix(
        Kop, ctx.X, ctx.X, ctx.tree;
        compressor    = ACA(; convergence = FNormEstimator(tol)),
        isnear        = AdaptiveCrossApproximation.isnear(; η = η),
        nearquadstrat = nearquadstrat,
        farquadstrat  = farquadstrat,
        scheduler     = scheduler,
        perm          = false,
    )
end

function build_fmm_D(ctx::CircleDLKPermCtx, k::ComplexF64;
                     tol::Float64 = 1e-4,
                     η::Float64   = 0.25,
                     scheduler = SerialScheduler(),
                     nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                     farquadstrat  = BEAST.DoubleNumQStrat(5,4))
    Kop = Helmholtz2D.doublelayer(; wavenumber = k)
    return AdaptiveCrossApproximation.HMatrix(
        Kop, ctx.Xtest_p, ctx.Xtrial_p, ctx.tree;
        compressor    = ACA(; convergence = FNormEstimator(tol)),
        isnear        = AdaptiveCrossApproximation.isnear(; η = η),
        nearquadstrat = nearquadstrat,
        farquadstrat  = farquadstrat,
        scheduler     = scheduler,
        perm          = false,
    )
end

# ============================================================
# T = D + jump_coeff * Iop
# Here jump_coeff = -0.5 matches your old dense file
# ============================================================

function apply_T_dense(Dmat::AbstractMatrix, Iop::AbstractMatrix, v;
                       jump_coeff::Float64 = -0.5)
    return Dmat * v .+ jump_coeff .* (Iop * v)
end

function dense_matrix_from_parts(Dmat::AbstractMatrix, Iop::AbstractMatrix;
                                 jump_coeff::Float64 = -0.5)
    return Dmat .+ jump_coeff .* Iop
end

function make_shifted_map(Dop, Iop::AbstractMatrix, D::Int;
                          jump_coeff::Float64 = -0.5)
    return LinearMap{ComplexF64}(
        x -> (Dop * x .+ jump_coeff .* (Iop * x)),
        D, D;
        ismutating = false
    )
end

# ============================================================
# ACA solve/apply for Beyn
# ============================================================

function Beyn_circle_D_mfie_tm_ACA(k::ComplexF64,
                                   VHat::AbstractMatrix{ComplexF64},
                                   MInvVHat::AbstractMatrix{ComplexF64},
                                   ctx;
                                   tol::Float64 = 1e-4,
                                   η::Float64   = 0.25,
                                   scheduler = SerialScheduler(),
                                   nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                                   farquadstrat  = BEAST.DoubleNumQStrat(5,4),
                                   jump_coeff::Float64 = -0.5,
                                   gmres_rtol::Float64 = 1e-10,
                                   gmres_itmax::Int = 1500,
                                   verbose::Bool = false)

    Dfmm = build_fmm_D(ctx, k;
        tol=tol, η=η, scheduler=scheduler,
        nearquadstrat=nearquadstrat, farquadstrat=farquadstrat)

    Iop = build_dense_Iop(ctx)
    Ddim = size(VHat, 1)
    Top = make_shifted_map(Dfmm, Iop, Ddim; jump_coeff=jump_coeff)

    nrhs = size(VHat, 2)
    for j in 1:nrhs
        b = view(VHat, :, j)
        x, st = Krylov.gmres(Top, b; rtol=gmres_rtol, itmax=gmres_itmax, verbose=0)

        if verbose && !st.solved
            @printf("GMRES warning at k=%+.6e%+.6ei, rhs=%d, status=%s, niter=%d\n",
                    real(k), imag(k), j, string(st.status), st.niter)
        end

        MInvVHat[:, j] .= x
    end

    return MInvVHat
end

function ApplyT_circle_D_mfie_tm_ACA(k::ComplexF64,
                                     v::AbstractVecOrMat{ComplexF64},
                                     ctx;
                                     tol::Float64 = 1e-4,
                                     η::Float64   = 0.25,
                                     scheduler = SerialScheduler(),
                                     nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                                     farquadstrat  = BEAST.DoubleNumQStrat(5,4),
                                     jump_coeff::Float64 = -0.5)
    Dfmm = build_fmm_D(ctx, k;
        tol=tol, η=η, scheduler=scheduler,
        nearquadstrat=nearquadstrat, farquadstrat=farquadstrat)

    Iop = build_dense_Iop(ctx)
    return Dfmm * v .+ jump_coeff .* (Iop * v)
end

# ============================================================
# Dense Beyn runner
# ============================================================

function run_beyn_dense_doublelayer_mfie_tm(;
        radius::Float64 = 1.0,
        Nel::Int = 512,
        z0::ComplexF64 = ComplexF64(2.4067799554075626, 0.0),
        Rx::Float64 = 0.05,
        Ry::Float64 = 0.02,
        jump_coeff::Float64 = -0.5,
        L::Int = 8,
        N::Int = 32,
        Kmax::Int = 2,
        tol_svd::Float64 = 1e-8,
        tol_res::Float64 = 1e-8,
        seed::Int = 1,
        verbose::Bool = true,
        quadstrat = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15))

    println("=== DENSE Beyn / double-layer / mfie_tm ===")
    @printf("jump_coeff = %.3f\n", jump_coeff)

    ctx = make_ctx_circle(; radius=radius, Nel=Nel, minvalues=100)

    Ddim = length(ctx.X)
    k_true = besselj_zero(0, 1) / radius

    @printf("D = %d\n", Ddim)
    @printf("Analytic k_(0,1) ≈ %.12f\n", k_true)

    UserFun = (k, Vhat, MIV) -> begin
        Dmat = build_dense_D(ctx, k; quadstrat=quadstrat)
        Iop  = build_dense_Iop(ctx)
        Tmat = dense_matrix_from_parts(Dmat, Iop; jump_coeff=jump_coeff)
        MIV .= Tmat \ Vhat
        MIV
    end

    ApplyT = (k, v) -> begin
        Dmat = build_dense_D(ctx, k; quadstrat=quadstrat)
        Iop  = build_dense_Iop(ctx)
        apply_T_dense(Dmat, Iop, v; jump_coeff=jump_coeff)
    end

    @time λ, V, ok = BeynSolve_auto(z0, Rx, Ry, UserFun, ApplyT, Ddim;
        L=L, N=N, Kmax=Kmax,
        tol_svd=tol_svd, tol_res=tol_res,
        seed=seed, verbose=verbose)

    p = sortperm(λ, by = x -> abs(x - k_true))
    λ = λ[p]

    println("ok = ", ok)
    println("Eigenvalues inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e %+.12ei   |k-k_true| ≈ %.3e\n",
                j, real(val), imag(val), abs(val - k_true))
    end

    return λ, V, ok
end

# ============================================================
# ACA original-basis Beyn runner
# ============================================================

function run_beyn_aca_original_doublelayer_mfie_tm(;
        radius::Float64 = 1.0,
        Nel::Int = 512,
        minvalues::Int = 100,
        tol::Float64 = 1e-4,
        η::Float64 = 0.25,
        scheduler = SerialScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        farquadstrat  = BEAST.DoubleNumQStrat(5,4),
        z0::ComplexF64 = ComplexF64(2.4067799554075626, 0.0),
        Rx::Float64 = 0.05,
        Ry::Float64 = 0.02,
        jump_coeff::Float64 = -0.5,
        L::Int = 8,
        N::Int = 32,
        Kmax::Int = 2,
        tol_svd::Float64 = 1e-8,
        tol_res::Float64 = 1e-5,
        gmres_rtol::Float64 = 1e-10,
        gmres_itmax::Int = 1500,
        seed::Int = 1,
        verbose::Bool = true)

    println("=== ACA Beyn / original basis / double-layer / mfie_tm ===")
    @printf("jump_coeff = %.3f\n", jump_coeff)

    ctx = make_ctx_circle(; radius=radius, Nel=Nel, minvalues=minvalues)

    Ddim = length(ctx.X)
    k_true = besselj_zero(0, 1) / radius

    @printf("D = %d\n", Ddim)
    @printf("Analytic k_(0,1) ≈ %.12f\n", k_true)
    @printf("ACA params: tol = %.1e, η = %.3f, scheduler = %s\n",
            tol, η, string(typeof(scheduler)))

    UserFun = (k, Vhat, MIV) -> Beyn_circle_D_mfie_tm_ACA(k, Vhat, MIV, ctx;
        tol=tol, η=η, scheduler=scheduler,
        nearquadstrat=nearquadstrat, farquadstrat=farquadstrat,
        jump_coeff=jump_coeff,
        gmres_rtol=gmres_rtol, gmres_itmax=gmres_itmax, verbose=false)

    ApplyT = (k, v) -> ApplyT_circle_D_mfie_tm_ACA(k, v, ctx;
        tol=tol, η=η, scheduler=scheduler,
        nearquadstrat=nearquadstrat, farquadstrat=farquadstrat,
        jump_coeff=jump_coeff)

    @time λ, V, ok = BeynSolve_auto(z0, Rx, Ry, UserFun, ApplyT, Ddim;
        L=L, N=N, Kmax=Kmax,
        tol_svd=tol_svd, tol_res=tol_res,
        seed=seed, verbose=verbose)

    p = sortperm(λ, by = x -> abs(x - k_true))
    λ = λ[p]

    println("ok = ", ok)
    println("Eigenvalues inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e %+.12ei   |k-k_true| ≈ %.3e\n",
                j, real(val), imag(val), abs(val - k_true))
    end

    return λ, V, ok
end

# ============================================================
# ACA bootstrap-permuted-basis Beyn runner
# ============================================================

function run_beyn_aca_bootstrap_perm_doublelayer_mfie_tm(;
        radius::Float64 = 1.0,
        Nel::Int = 512,
        minvalues::Int = 100,
        tol::Float64 = 1e-4,
        η::Float64 = 0.25,
        scheduler = SerialScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        farquadstrat  = BEAST.DoubleNumQStrat(5,4),
        z0::ComplexF64 = ComplexF64(2.4067799554075626, 0.0),
        Rx::Float64 = 0.05,
        Ry::Float64 = 0.02,
        jump_coeff::Float64 = -0.5,
        L::Int = 8,
        N::Int = 32,
        Kmax::Int = 2,
        tol_svd::Float64 = 1e-8,
        tol_res::Float64 = 1e-5,
        gmres_rtol::Float64 = 1e-10,
        gmres_itmax::Int = 1500,
        seed::Int = 1,
        verbose::Bool = true)

    println("=== ACA Beyn / bootstrap-permuted basis / double-layer / mfie_tm ===")
    @printf("jump_coeff = %.3f\n", jump_coeff)

    ctxp = make_ctx_circle_perm_bootstrap(; radius=radius, Nel=Nel, minvalues=minvalues,
        kref=z0, tol=tol, η=η, scheduler=scheduler,
        nearquadstrat=nearquadstrat, farquadstrat=farquadstrat)

    Ddim = length(ctxp.Xtest_p)
    k_true = besselj_zero(0, 1) / radius

    @printf("D = %d\n", Ddim)
    @printf("Analytic k_(0,1) ≈ %.12f\n", k_true)
    @printf("ACA bootstrap-perm params: tol = %.1e, η = %.3f, scheduler = %s\n",
            tol, η, string(typeof(scheduler)))

    UserFun = (k, Vhat, MIV) -> Beyn_circle_D_mfie_tm_ACA(k, Vhat, MIV, ctxp;
        tol=tol, η=η, scheduler=scheduler,
        nearquadstrat=nearquadstrat, farquadstrat=farquadstrat,
        jump_coeff=jump_coeff,
        gmres_rtol=gmres_rtol, gmres_itmax=gmres_itmax, verbose=false)

    ApplyT = (k, v) -> ApplyT_circle_D_mfie_tm_ACA(k, v, ctxp;
        tol=tol, η=η, scheduler=scheduler,
        nearquadstrat=nearquadstrat, farquadstrat=farquadstrat,
        jump_coeff=jump_coeff)

    @time λ, V, ok = BeynSolve_auto(z0, Rx, Ry, UserFun, ApplyT, Ddim;
        L=L, N=N, Kmax=Kmax,
        tol_svd=tol_svd, tol_res=tol_res,
        seed=seed, verbose=verbose)

    p = sortperm(λ, by = x -> abs(x - k_true))
    λ = λ[p]

    println("ok = ", ok)
    println("Eigenvalues inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e %+.12ei   |k-k_true| ≈ %.3e\n",
                j, real(val), imag(val), abs(val - k_true))
    end

    return λ, V, ok
end

# ============================================================
# Comparison helper
# ============================================================

function print_compare_table(λd, λa0, λap; kref)
    println("\n=== COMPARISON TO REFERENCE ===")
    nd  = length(λd)
    na0 = length(λa0)
    nap = length(λap)

    println("dense count                = ", nd)
    println("ACA original basis count   = ", na0)
    println("ACA bootstrap-perm count   = ", nap)

    println("\nClosest eigenvalues to reference:")
    if !isempty(λd)
        j = argmin(abs.(λd .- kref))
        @printf("dense               : %+.12e %+.12ei   |·-kref|=%.3e\n",
                real(λd[j]), imag(λd[j]), abs(λd[j]-kref))
    end
    if !isempty(λa0)
        j = argmin(abs.(λa0 .- kref))
        @printf("ACA original basis  : %+.12e %+.12ei   |·-kref|=%.3e\n",
                real(λa0[j]), imag(λa0[j]), abs(λa0[j]-kref))
    end
    if !isempty(λap)
        j = argmin(abs.(λap .- kref))
        @printf("ACA bootstrap-perm  : %+.12e %+.12ei   |·-kref|=%.3e\n",
                real(λap[j]), imag(λap[j]), abs(λap[j]-kref))
    end
end

# ============================================================
# Main run block
# ============================================================

function run_compare_doublelayer_mfie_tm_beyn()
    radius = 1.0
    Nel = 512

    z0 = ComplexF64(2.4067799554075626, 0.0)
    Rx = 0.05
    Ry = 0.02

    tol = 1e-4
    η   = 0.25
    jump_coeff = -0.5

    L = 8
    N = 32
    Kmax = 2

    k_true = besselj_zero(0,1) / radius

    println("\n================ DENSE =================\n")
    λd, Vd, okd = run_beyn_dense_doublelayer_mfie_tm(
        radius = radius,
        Nel = Nel,
        z0 = z0,
        Rx = Rx,
        Ry = Ry,
        jump_coeff = jump_coeff,
        L = L,
        N = N,
        Kmax = Kmax,
        tol_svd = 1e-8,
        tol_res = 1e-8,
        seed = 1,
        verbose = true,
        quadstrat = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
    )

    println("\n================ ACA ORIGINAL BASIS =================\n")
    λa0, Va0, oka0 = run_beyn_aca_original_doublelayer_mfie_tm(
        radius = radius,
        Nel = Nel,
        minvalues = 100,
        tol = tol,
        η = η,
        scheduler = SerialScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        farquadstrat  = BEAST.DoubleNumQStrat(5,4),
        z0 = z0,
        Rx = Rx,
        Ry = Ry,
        jump_coeff = jump_coeff,
        L = L,
        N = N,
        Kmax = Kmax,
        tol_svd = 1e-8,
        tol_res = 1e-5,
        gmres_rtol = 1e-10,
        gmres_itmax = 1500,
        seed = 1,
        verbose = true,
    )

    println("\n================ ACA BOOTSTRAP-PERM BASIS =================\n")
    λap, Vap, okap = run_beyn_aca_bootstrap_perm_doublelayer_mfie_tm(
        radius = radius,
        Nel = Nel,
        minvalues = 100,
        tol = tol,
        η = η,
        scheduler = SerialScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        farquadstrat  = BEAST.DoubleNumQStrat(5,4),
        z0 = z0,
        Rx = Rx,
        Ry = Ry,
        jump_coeff = jump_coeff,
        L = L,
        N = N,
        Kmax = Kmax,
        tol_svd = 1e-8,
        tol_res = 1e-5,
        gmres_rtol = 1e-10,
        gmres_itmax = 1500,
        seed = 1,
        verbose = true,
    )

    print_compare_table(λd, λa0, λap; kref=k_true)

    return λd, λa0, λap
end

λd, λa0, λap = run_compare_doublelayer_mfie_tm_beyn()