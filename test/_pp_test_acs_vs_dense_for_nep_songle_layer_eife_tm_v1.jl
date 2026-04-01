# pp_test_beyn_aca_single_layer_efie_tm.jl
println("RUNNING FILE: pp_test_beyn_aca_single_layer_efie_tm.jl")

using LinearAlgebra, Random, Statistics, Printf
using BEAST, CompScienceMeshes
using AdaptiveCrossApproximation
using OhMyThreads
using Krylov
using H2Trees
using SpecialFunctions
using LinearMaps
using FunctionZeros
using Base.Threads

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
# Geometry / tree / context
# ============================================================

struct CircleEFIECtx
    X
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

    p_order = 2
    X = lagrangecx(Γ, order = p_order)

    tree = build_blocktree(X; minvalues=minvalues)
    return CircleEFIECtx(X, tree)
end

# ============================================================
# Dense single-layer operator (optional reference / residual)
# ============================================================

function build_dense_S(ctx::CircleEFIECtx, k::ComplexF64;
                       quadstrat = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15))
    Sop = Helmholtz2D.singlelayer(; wavenumber = k)
    S   = assemble(Sop, ctx.X, ctx.X; quadstrat = quadstrat)
    return Matrix(S)
end

# ============================================================
# ACA/HMatrix single-layer operator
# ============================================================

function build_fmm_S(ctx::CircleEFIECtx, k::ComplexF64;
                     tol::Float64 = 1e-4,
                     η::Float64   = 0.25,
                     scheduler     = SerialScheduler() ,
                     nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                     farquadstrat  = BEAST.DoubleNumQStrat(4,4))

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
# ACA linear solve S(k) Y = VHat
# solve each RHS with GMRES
# ============================================================

function Beyn_circle_D_ACA(k::ComplexF64,
                           VHat::AbstractMatrix{ComplexF64},
                           MInvVHat::AbstractMatrix{ComplexF64},
                           ctx::CircleEFIECtx;
                           tol::Float64 = 1e-4,
                           η::Float64   = 0.25,
                           scheduler     = SerialScheduler() ,
                           nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                           farquadstrat  = BEAST.DoubleNumQStrat(4,4),
                           gmres_rtol::Float64 = 1e-10,
                           gmres_itmax::Int = 1500,
                           verbose::Bool = false)

    S = build_fmm_S(ctx, k;
                    tol = tol,
                    η = η,
                    scheduler = scheduler,
                    nearquadstrat = nearquadstrat,
                    farquadstrat = farquadstrat)

    nrhs = size(VHat, 2)

    for j in 1:nrhs
        b = view(VHat, :, j)
        x, st = Krylov.gmres(S, b; rtol=gmres_rtol, itmax=gmres_itmax, verbose=0)

        if verbose && !st.solved
            @printf("GMRES warning at k=%+.6e%+.6ei, rhs=%d, status=%s, niter=%d\n",
                    real(k), imag(k), j, string(st.status), st.niter)
        end

        MInvVHat[:, j] .= x
    end

    return MInvVHat
end

function ApplyT_circle_D_ACA(k::ComplexF64,
                             v::AbstractVecOrMat{ComplexF64},
                             ctx::CircleEFIECtx;
                             tol::Float64 = 1e-4,
                             η::Float64   = 0.25,
                             scheduler     = SerialScheduler() ,
                             nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                             farquadstrat  = BEAST.DoubleNumQStrat(4,4))
    S = build_fmm_S(ctx, k;
                    tol = tol,
                    η = η,
                    scheduler = scheduler,
                    nearquadstrat = nearquadstrat,
                    farquadstrat = farquadstrat)
    return S * v
end

# ============================================================
# Beyn solver
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

#   nT = nthreads()
#   Ap_thread   = [[zeros(ComplexF64, D, L) for _ in 1:nMom_max] for _ in 1:nT]
#   MInv_thread = [zeros(ComplexF64, D, L) for _ in 1:nT]
#
#    Δθ = 2π / N
#
#    @threads for nn in 0:(N-1)
#        tid = threadid()
#        Ap_loc   = Ap_thread[tid]
#        MInvVHat = MInv_thread[tid]
#
#        θ  = nn * Δθ
#        cθ = cos(θ)
#        sθ = sin(θ)
#
#        w  = Rx * cθ + im * Ry * sθ
#        z  = z0 + w
#        dz = (-Rx * sθ + im * Ry * cθ) * Δθ
#
#        UserBeynFunction(z, VHat, MInvVHat)
#
#        wp = one(w)
#        @inbounds for p in 1:nMom_max
#            Ap_loc[p] .+= wp * dz .* MInvVHat
#            wp *= w
#        end
#    end
#
#    for t in 1:nT
#        Ap_loc = Ap_thread[t]
#        @inbounds for p in 1:nMom_max
#            Ap[p] .+= Ap_loc[p]
#        end
#    end

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
            rows = (i-1)*m + 1 : i*m
            for j in 1:K
                cols = (j-1)*ℓ + 1 : j*ℓ
                p0 = i + j - 2
                p1 = i + j - 1
                B0[rows, cols] .= Ap[p0 + 1]
                B1[rows, cols] .= Ap[p1 + 1]
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

        if ApplyT !== nothing
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
        else
            good_idxs = collect(eachindex(λ))
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
# ACA Beyn test
# ============================================================

function test_circle_cavity_Dirichlet_Beyn_ACA(;
        radius::Float64 = 1.0,
        Nel::Int = 512,
        minvalues::Int = 100,
        tol::Float64 = 1e-4,
        η::Float64 = 0.25,
        scheduler = StaticScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        farquadstrat  = BEAST.DoubleNumQStrat(4,4),
        z0 = nothing,
        R::Float64 = 0.4,
        L::Int = 8,
        N::Int = 32,
        Kmax::Int = 2,
        tol_svd::Float64 = 1e-8,
        tol_res::Float64 = 1e-8,
        gmres_rtol::Float64 = 1e-10,
        gmres_itmax::Int = 1500,
        seed::Int = 1,
        verbose::Bool = true)

    println("=== Beyn eigenvalue test: circular cavity (Dirichlet, interior, EFIE, ACA) ===")

    ctx = make_ctx_circle(; radius=radius, Nel=Nel, minvalues=minvalues)
    D = length(ctx.X)

    k_true = besselj_zero(0, 1) / radius
    z0 === nothing && (z0 = complex(k_true))

    @printf("D = %d\n", D)
    @printf("Analytic k_(0,1) ≈ %.10f\n", k_true)
    @printf("ACA params: tol = %.1e, η = %.3f, scheduler = %s\n",
            tol, η, string(typeof(scheduler)))

    UserFun = (k, Vhat, MIV) -> Beyn_circle_D_ACA(k, Vhat, MIV, ctx;
        tol=tol,
        η=η,
        scheduler=scheduler,
        nearquadstrat=nearquadstrat,
        farquadstrat=farquadstrat,
        gmres_rtol=gmres_rtol,
        gmres_itmax=gmres_itmax,
        verbose=false)

    ApplyT = (k, v) -> ApplyT_circle_D_ACA(k, v, ctx;
        tol=tol,
        η=η,
        scheduler=scheduler,
        nearquadstrat=nearquadstrat,
        farquadstrat=farquadstrat)

    @time λ, V, ok = BeynSolve_auto(z0, R, UserFun, ApplyT, D;
        L = L,
        N = N,
        Kmax = Kmax,
        tol_svd = tol_svd,
        tol_res = tol_res,
        seed = seed,
        verbose = verbose)

    println("ok = ", ok)

    p = sortperm(λ, by = x -> (real(x), imag(x)))
    λ = λ[p]

    println("Eigenvalues inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e %+.12ei   |k-k_true| ≈ %.3e\n",
                j, real(val), imag(val), abs(val - k_true))
    end

    return λ, V, ok
end

# ============================================================
# Run
# ============================================================

λa, Va, oka = test_circle_cavity_Dirichlet_Beyn_ACA(
    radius = 1.0,
    Nel = 512,
    minvalues = 100,
    tol = 1e-4,
    η = 0.25,
    scheduler = SerialScheduler(),
    nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
    farquadstrat  = BEAST.DoubleNumQStrat(4,4),
    R = 0.4,
    L = 8,
    N = 32,
    Kmax = 2,
    tol_svd = 1e-8,
    tol_res = 1e-8,
    gmres_rtol = 1e-10,
    gmres_itmax = 1500,
    seed = 1,
    verbose = true,
)


# λa, Va, oka = test_circle_cavity_Dirichlet_Beyn_ACA(
#     radius = 1.0,
#     Nel = 512,
#     minvalues = 100,
#     tol = 1e-4,
#     η = 0.25,
#     scheduler = StaticScheduler(),
#     nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
#     farquadstrat  = BEAST.DoubleNumQStrat(4,4),
#     R = 0.4,
#     L = 8,
#     N = 32,
#     Kmax = 2,
#     tol_svd = 1e-8,
#     tol_res = 1e-8,
#     gmres_rtol = 1e-10,
#     gmres_itmax = 1500,
#     seed = 1,
#     verbose = true,
# )