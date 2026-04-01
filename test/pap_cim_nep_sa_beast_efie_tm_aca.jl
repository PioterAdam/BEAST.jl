using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test
using SpecialFunctions
using Gmsh
using Printf
using Base.Threads
using FunctionZeros

using ParallelKMeans
using H2Trees
using AdaptiveCrossApproximation
using OhMyThreads
using Krylov

const Shape = BEAST.Shape
const LagrangeBasis = BEAST.LagrangeBasis

import LinearAlgebra: norm, dot
import BEAST: kernelvals, KernelValsHelmholtz2D
import SpecialFunctions: hankelh2
import CompScienceMeshes: normal

const j01 = besselj_zero(0, 1)

# ============================================================
# Complex-k patch for Helmholtz 2D kernels
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
# Dense EFIE operator for 2D TM interior Dirichlet
# ============================================================

"""
    assemble_T_circle_Dirichlet(k, X)

Dense EFIE / single-layer operator for the 2D TM interior Dirichlet cavity:
    T(k) = S(k)
"""
function assemble_T_circle_Dirichlet(k::Number, X)
    𝒮 = Helmholtz2D.singlelayer(; wavenumber = k)
    S  = assemble(𝒮, X, X)
    return Matrix(S)
end

function Beyn_circle_D_dense(k, VHat, MInvVHat, X)
    T = assemble_T_circle_Dirichlet(k, X)
    MInvVHat .= T \ VHat
    return MInvVHat
end

function ApplyT_circle_D_dense(k::ComplexF64,
                               v::AbstractVecOrMat{ComplexF64},
                               X)
    T = assemble_T_circle_Dirichlet(k, X)
    return T * v
end

# ============================================================
# Beyn solver
# ============================================================

function BeynSolve_auto(z0, Rx, Ry, UserBeynFunction, ApplyT, D;
                        L::Int       = D,
                        N::Int       = 64,
                        Kmax::Int    = 3,
                        tol_svd::Float64 = 1e-8,
                        tol_res::Float64 = 1e-8,
                        verbose::Bool = true)

    # Fixed probing matrix
    VHat = randn(D, L) .+ im * randn(D, L)

    # Sequential version for robustness/debugging
    nMom_max = 2*Kmax
    Ap       = [zeros(ComplexF64, D, L) for _ in 1:nMom_max]
    MInvVHat = zeros(ComplexF64, D, L)

    Δθ = 2π / N
    for n = 0:N-1
        θ  = n * Δθ
        cθ = cos(θ)
        sθ = sin(θ)

        w  = Rx * cθ + im * Ry * sθ
        z  = z0 + w
        dz = (-Rx * sθ + im * Ry * cθ) * Δθ

        MInvVHat = UserBeynFunction(z, VHat, MInvVHat)

        wp = one(w)
        @inbounds for p = 1:nMom_max
            Ap[p] .+= wp * dz .* MInvVHat
            wp *= w
        end
    end

    function Beyn_from_moments(K::Int)
        m    = D
        ℓ    = L
        Km   = K*m
        Kl   = K*ℓ

        B0 = zeros(ComplexF64, Km, Kl)
        B1 = zeros(ComplexF64, Km, Kl)

        for i in 1:K
            rows = (i-1)*m+1 : i*m
            for j in 1:K
                cols = (j-1)*ℓ+1 : j*ℓ
                p0 = i + j - 2
                p1 = i + j - 1
                @inbounds begin
                    B0[rows, cols] .= Ap[p0+1]
                    B1[rows, cols] .= Ap[p1+1]
                end
            end
        end

        U, σ, W = svd(B0)
        k = count(>(tol_svd), abs.(σ))
        if verbose
            @printf("K=%d: rank(B0) ≈ k = %d\n", K, k)
        end
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

        r  = zeros(Float64, length(λ))
        good_idxs = Int[]

        if ApplyT === nothing
            return λ, Vfull, collect(1:length(λ)), r
        else
            if verbose
                @printf("  Residuals (K=%d):\n", K)
            end
            for j in 1:length(λ)
                vj = Vfull[:, j]
                Tv = ApplyT(λ[j], vj)
                rj = norm(Tv) / norm(vj)
                r[j] = rj
                if verbose
                    @printf("    res[%2d] = %.3e\n", j, rj)
                end
                if rj <= tol_res
                    push!(good_idxs, j)
                end
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
        λ, V, good_idxs, r = Beyn_from_moments(K)
        ngood = length(good_idxs)

        if verbose
            @printf("  K=%d: %d good eigenpairs (res ≤ %.1e)\n", K, ngood, tol_res)
        end

        if ngood == 0
            continue
        end

        λg = λ[good_idxs]
        Vg = V[:, good_idxs]

        if good_best == 0
            good_best = ngood
            λ_best    = λg
            V_best    = Vg
            K_used    = K
        end

        if ngood > good_prev
            good_best = ngood
            λ_best    = λg
            V_best    = Vg
            K_used    = K
        elseif ngood == good_prev && ngood > 0
            if verbose
                @printf("  #good eigenpairs stabilized at %d for K=%d.\n", ngood, K)
            end
            λ_best = λg
            V_best = Vg
            K_used = K
            break
        end

        good_prev = ngood
    end

    if good_best == 0 || isempty(λ_best)
        if verbose
            @printf("No acceptable eigenpairs found.\n")
        end
        return ComplexF64[], zeros(ComplexF64, D, 0), false
    end

    p = sortperm(λ_best, by = x -> (real(x), imag(x)))
    λ_best = λ_best[p]
    V_best = V_best[:, p]

    if verbose
        @printf("Using K=%d with %d good eigenpairs.\n", K_used, good_best)
    end

    return λ_best, V_best, true
end

function BeynSolve_auto(z0, R, UserBeynFunction, ApplyT, D; kwargs...)
    return BeynSolve_auto(z0, R, R, UserBeynFunction, ApplyT, D; kwargs...)
end

# ============================================================
# Asakura / block Sakurai-Sugiura
# ============================================================

function AsakuraSolve_circle_auto(z0::Complex, R::Float64,
                                  UserFun, ApplyT, D::Int;
                                  L::Int         = 8,
                                  Kblocks::Int   = 4,
                                  N::Int         = 64,
                                  δ_svd::Float64 = 1e-12,
                                  tol_res::Float64 = 1e-8,
                                  verbose::Bool  = true)

    Kblocks ≥ 1 || error("Kblocks must be ≥ 1")
    L       ≥ 1 || error("L must be ≥ 1")
    N       ≥ 4 || error("N must be ≥ 4")

    nMom = 2*Kblocks

    verbose && println("=== Asakura / block-SS contour solver (circle) ===")
    verbose && @printf("Center z0 = %+.6e%+.6ei, R = %.3f, D = %d, L = %d, Kblocks = %d, N = %d\n",
                       real(z0), imag(z0), R, D, L, Kblocks, N)

    V        = randn(ComplexF64, D, L)
    S_blocks = [zeros(ComplexF64, D, L) for _ in 1:nMom]
    MInv     = zeros(ComplexF64, D, L)

    Δθ = 2π / N

    for n = 0:N-1
        θ  = (n + 0.5) * Δθ
        ζ  = exp(im * θ)
        z  = z0 + R * ζ
        dz = im * R * ζ * Δθ

        w_dθ = dz / (2π * im)
        w = z - z0

        G = UserFun(z, V, MInv)

        wpow = one(w)
        for k = 0:(nMom-1)
            S_blocks[k+1] .+= (wpow * w_dθ) .* G
            wpow *= w
        end
    end

    rows_H = Kblocks * D
    cols_H = Kblocks * L

    H   = zeros(ComplexF64, rows_H, cols_H)
    Hlt = zeros(ComplexF64, rows_H, cols_H)

    for i in 1:Kblocks
        for j in 1:Kblocks
            Sk   = S_blocks[i + j - 1]
            Sk1  = S_blocks[i + j]
            rows = (i-1)*D+1 : i*D
            cols = (j-1)*L+1 : j*L
            @inbounds begin
                H[rows, cols]   .= Sk
                Hlt[rows, cols] .= Sk1
            end
        end
    end

    if verbose
        σvals = svdvals(H)
        σmax  = maximum(σvals)
        @printf("SVD(H): σ_max = %.3e, σ_min = %.3e\n", σmax, minimum(σvals))
    end

    UH, σ, VH = svd(H)
    σmax = maximum(σ)
    keep = findall(σ .>= δ_svd * σmax)
    mprime = length(keep)

    if mprime == 0
        verbose && @warn "All singular values below threshold; no eigenvalues returned."
        return ComplexF64[], zeros(ComplexF64, D, 0), false
    end

    UH1 = UH[:, keep]
    VH1 = VH[:, keep]
    Σ1  = Diagonal(σ[keep])

    Hlt_tilde = UH1' * Hlt * VH1
    G         = Σ1 \ Hlt_tilde

    eigG = eigen(G)
    μ    = eigG.values
    Y    = eigG.vectors

    λ_hat = z0 .+ μ

    S_matrix = hcat(S_blocks[1:Kblocks]...)
    W_orig   = VH1 * Y

    X_hat = Matrix{ComplexF64}(undef, D, mprime)
    for j in 1:mprime
        xj = S_matrix * W_orig[:, j]
        n  = norm(xj)
        X_hat[:, j] .= n > 0 ? xj ./ n : xj
    end

    good    = Bool[]
    resvals = Float64[]

    for j in 1:mprime
        v = X_hat[:, j]
        r = ApplyT(λ_hat[j], v)
        rnorm = norm(r)
        vnorm = norm(v)
        ρj = vnorm == 0 ? Inf : rnorm / vnorm
        push!(resvals, ρj)
        push!(good, ρj ≤ tol_res)
    end

    if verbose
        println("Residual check (tol_res = $(@sprintf("%.3e", tol_res))):")
        for j in 1:mprime
            @printf("  j=%2d: |λ̂|=% .6e, res = %.3e%s\n",
                    j, abs(λ_hat[j]), resvals[j],
                    good[j] ? "  (good)" : "")
        end
    end

    idx_good = findall(good)
    if isempty(idx_good)
        verbose && @warn "No eigenpairs passed the residual tolerance."
        return λ_hat, X_hat, false
    end

    λ_good = λ_hat[idx_good]
    X_good = X_hat[:, idx_good]

    return λ_good, X_good, true
end

# ============================================================
# ACA / FMM helper layer for EFIE
# ============================================================

function build_blocktree(X; minvalues=100)
    testtree  = KMeansTree(X.pos, 2; minvalues=minvalues)
    trialtree = KMeansTree(X.pos, 2; minvalues=minvalues)
    return BlockTree(testtree, trialtree)
end

"""
Simple ACA/HMatrix assembly for EFIE.
Start conservatively with eta ≈ 1.
"""
function hassemble_with_tree(op, X, tree;
                             aca_tol::Float64 = 1e-6,
                             η::Float64 = 1.0,
                             nearquadstrat = BEAST.DoubleNumSauterQstrat(30,30,0,4,25,25),
                             farquadstrat  = BEAST.DoubleNumQStrat(7,8),
                             scheduler     = DynamicScheduler())

    return AdaptiveCrossApproximation.HMatrix(
        op, X, X, tree;
        compressor    = ACA(; tol = aca_tol),
        isnear        = AdaptiveCrossApproximation.isnear(; η = η),
        nearquadstrat = nearquadstrat,
        farquadstrat  = farquadstrat,
        scheduler     = scheduler,
    )
end

mutable struct CircleEFIECtx
    X
    tree
    Scache::Dict{ComplexF64, Any}
end

function CircleEFIECtx(X; minvalues=100)
    tree = build_blocktree(X; minvalues=minvalues)
    return CircleEFIECtx(X, tree, Dict{ComplexF64,Any}())
end

function getS!(ctx::CircleEFIECtx, k::ComplexF64;
               aca_tol::Float64 = 1e-6,
               η::Float64 = 1.0,
               nearquadstrat = BEAST.DoubleNumSauterQstrat(30,30,0,4,25,25),
               farquadstrat  = BEAST.DoubleNumQStrat(7,8),
               scheduler     = DynamicScheduler())

    get!(ctx.Scache, k) do
        𝒮 = Helmholtz2D.singlelayer(; wavenumber = k)
        @time hassemble_with_tree(
            𝒮, ctx.X, ctx.tree;
            aca_tol       = aca_tol,
            η             = η,
            nearquadstrat = nearquadstrat,
            farquadstrat  = farquadstrat,
            scheduler     = scheduler,
        )
    end
end

function Beyn_circle_D(k, VHat, MInvVHat, ctx::CircleEFIECtx;
                       rtol=1e-6, itmax=800, verbose=0,
                       aca_tol=1e-6, η=1.0,
                       nearquadstrat = BEAST.DoubleNumSauterQstrat(30,30,0,4,25,25),
                       farquadstrat  = BEAST.DoubleNumQStrat(7,8),
                       scheduler     = DynamicScheduler())

    k = ComplexF64(k)
    Sk = getS!(ctx, k;
               aca_tol       = aca_tol,
               η             = η,
               nearquadstrat = nearquadstrat,
               farquadstrat  = farquadstrat,
               scheduler     = scheduler)

    n, L = size(VHat)
    @assert size(MInvVHat) == (n, L)

    for j in 1:L
        x, st = Krylov.gmres(Sk, view(VHat, :, j); rtol=rtol, itmax=itmax, verbose=verbose)
        st.solved || @warn "GMRES not solved" k=k j=j status=st.status niter=st.niter
        MInvVHat[:, j] .= x
    end
    return MInvVHat
end

function ApplyT_num(k::ComplexF64, V, ctx::CircleEFIECtx;
                    aca_tol=1e-6, η=1.0,
                    nearquadstrat = BEAST.DoubleNumSauterQstrat(30,30,0,4,25,25),
                    farquadstrat  = BEAST.DoubleNumQStrat(7,8),
                    scheduler     = DynamicScheduler())

    Sk = getS!(ctx, k;
               aca_tol       = aca_tol,
               η             = η,
               nearquadstrat = nearquadstrat,
               farquadstrat  = farquadstrat,
               scheduler     = scheduler)

    if V isa AbstractVector
        y = similar(V, ComplexF64)
        mul!(y, Sk, V)
        return y
    else
        @assert size(V, 1) == size(ctx.X.pos, 2) || true
        Y = similar(V, ComplexF64)
        for j in axes(V, 2)
            vj = view(V, :, j)
            yj = view(Y, :, j)
            mul!(yj, Sk, vj)
        end
        return Y
    end
end

# ============================================================
# Dense EFIE tests
# ============================================================

function test_circle_cavity_Dirichlet_Beyn_dense_EFIE()
    println("=== Beyn eigenvalue test: circular cavity (Dirichlet, interior, dense EFIE) ===")

    a       = 1.0
    p_order = 2
    h       = 2π * a / 64

    circle = CompScienceMeshes.meshcircle(a, h)
    X = lagrangecx(circle, order=p_order)

    T0 = assemble_T_circle_Dirichlet(2.0 + 0im, X)
    D  = size(T0, 1)
    @printf("D = %d\n", D)

    k_true = j01 / a
    @printf("Analytic k_{0,1} ≈ %.10f\n", k_true)

    z0 = k_true + 0im
    R  = 0.4

    L       = 8
    N       = 32
    Kmax    = 2
    tol_svd = 1e-8
    tol_res = 1e-8
    verbose = true

    UserFun = (k, Vhat, MIV) -> Beyn_circle_D_dense(k, Vhat, MIV, X)
    ApplyT  = (k, v)         -> ApplyT_circle_D_dense(ComplexF64(k), v, X)

    @time λ, V, ok = BeynSolve_auto(
        z0, R,
        UserFun, ApplyT, D;
        L       = L,
        N       = N,
        Kmax    = Kmax,
        tol_svd = tol_svd,
        tol_res = tol_res,
        verbose = verbose
    )

    println("ok = ", ok)

    p = sortperm(λ, by = x -> (real(x), imag(x)))
    λ = λ[p]
    V = V[:, p]

    println("Eigenvalues k inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e, %+.12e   |k-ktrue|=%.3e\n",
                j, real(val), imag(val), abs(val-k_true))
    end

    return λ, V, ok
end

function test_circle_cavity_Dirichlet_Asakura_dense_EFIE(;
        a::Float64      = 1.0,
        p_order::Int    = 2,
        h::Float64      = 2π * a / 64,
        L::Int          = 8,
        Kblocks::Int    = 2,
        N::Int          = 64,
        δ_svd::Float64  = 1e-12,
        tol_res::Float64 = 1e-8,
        verbose::Bool   = true)

    println("=== Asakura eigenvalue test: circular cavity (Dirichlet, interior, dense EFIE) ===")

    circle = CompScienceMeshes.meshcircle(a, h)
    X = lagrangecx(circle, order=p_order)

    T0 = assemble_T_circle_Dirichlet(2.0 + 0im, X)
    D  = size(T0, 1)
    @printf("D = %d\n", D)

    k_true = j01 / a
    @printf("Analytic k_{0,1} ≈ %.10f\n", k_true)

    z0 = k_true + 0im
    R  = 0.4

    UserFun = (k, Vhat, MInvVHat) -> Beyn_circle_D_dense(k, Vhat, MInvVHat, X)
    ApplyT  = (k, v)              -> ApplyT_circle_D_dense(ComplexF64(k), v, X)

    @time λ, V, ok = AsakuraSolve_circle_auto(
        z0, R,
        UserFun, ApplyT, D;
        L        = L,
        Kblocks  = Kblocks,
        N        = N,
        δ_svd    = δ_svd,
        tol_res  = tol_res,
        verbose  = verbose,
    )

    println("ok = ", ok)

    p = sortperm(λ, by = x -> (real(x), imag(x)))
    λ = λ[p]
    V = V[:, p]

    println("Eigenvalues k inside contour (Asakura):")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e, %+.12e   |k-ktrue|=%.3e\n",
                j, real(val), imag(val), abs(val-k_true))
    end

    return λ, V, ok
end

# ============================================================
# ACA / mFMM EFIE tests
# ============================================================

function test_circle_cavity_Dirichlet_Beyn_mFMM_EFIE()
    println("=== Beyn eigenvalue test: circular cavity (Dirichlet, interior, EFIE, mFMM) ===")

    a       = 1.0
    p_order = 2
    h       = 2π * a / (2 * 246)

    circle = CompScienceMeshes.meshcircle(a, h)
    X      = lagrangecx(circle, order=p_order)

    ctx = CircleEFIECtx(X; minvalues=100)

    D = length(X.fns)
    @printf("D = %d\n", D)

    k_true = j01 / a
    @printf("Analytic k_{0,1} ≈ %.10f\n", k_true)

    z0 = k_true + 0im
    R  = 0.10

    L       = 12
    N       = 64
    Kmax    = 2
    tol_svd = 1e-8
    tol_res = 1e-5
    verbose = true

    aca_tol = 1e-6
    η       = 1.0
    nearqs  = BEAST.DoubleNumSauterQstrat(30,30,0,4,25,25)
    farqs   = BEAST.DoubleNumQStrat(7,8)
    sched   = DynamicScheduler()

    UserFun = (k, Vhat, MIV) -> Beyn_circle_D(
        k, Vhat, MIV, ctx;
        rtol          = 1e-6,
        itmax         = 800,
        verbose       = 0,
        aca_tol       = aca_tol,
        η             = η,
        nearquadstrat = nearqs,
        farquadstrat  = farqs,
        scheduler     = sched
    )

    ApplyT = (k, V) -> ApplyT_num(
        ComplexF64(k), V, ctx;
        aca_tol       = aca_tol,
        η             = η,
        nearquadstrat = nearqs,
        farquadstrat  = farqs,
        scheduler     = sched
    )

    # consistency check
    kcheck = ComplexF64(z0 + R * exp(im * 0.37))
    Ltest  = 3
    Vhat   = randn(ComplexF64, D, Ltest)
    MIV    = similar(Vhat)

    UserFun(kcheck, Vhat, MIV)
    Tmiv = ApplyT(kcheck, MIV)
    err  = norm(Tmiv - Vhat) / norm(Vhat)
    @printf("Consistency check: ||T*MIV - Vhat||/||Vhat|| = %.3e\n", err)

    @time λ, V, ok = BeynSolve_auto(
        z0, R,
        UserFun, ApplyT, D;
        L       = L,
        N       = N,
        Kmax    = Kmax,
        tol_svd = tol_svd,
        tol_res = tol_res,
        verbose = verbose
    )

    println("ok = ", ok)
    if !ok
        println("Warning: BeynSolve_auto did not find a fully verified set of eigenpairs.")
    end

    p = sortperm(λ, by = x -> (real(x), imag(x)))
    λ = λ[p]
    V = V[:, p]

    println("Eigenvalues k inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e, %+.12e   |k-ktrue|=%.3e\n",
                j, real(val), imag(val), abs(val-k_true))
    end

    if !isempty(λ)
        jbest = argmin(abs.(λ .- k_true))
        khat  = λ[jbest]
        vhat  = V[:, jbest]

        Sfmm = getS!(ctx, ComplexF64(khat);
                     aca_tol       = aca_tol,
                     η             = η,
                     nearquadstrat = nearqs,
                     farquadstrat  = farqs,
                     scheduler     = sched)
        rfmm = norm(ApplyT(ComplexF64(khat), vhat)) / norm(vhat)

        Sop = Helmholtz2D.singlelayer(; wavenumber = ComplexF64(khat))
        nearquadstrat_dense = BEAST.DoubleNumSauterQstrat(8, 9, 0, 4, 9, 10)
        Sdense = assemble(Sop, X, X; quadstrat = nearquadstrat_dense)
        rdense = norm(Sdense * vhat) / norm(vhat)

        @printf("Residual(best): khat=%+.6f%+.6fi |khat-ktrue|=%.3e rfmm=%.3e rdense=%.3e\n",
                real(khat), imag(khat), abs(khat-k_true), rfmm, rdense)

        ratio = rfmm / max(rdense, 1e-300)
        @printf("  residual ratio rfmm/rdense = %.3e\n", ratio)

        Tv_fmm   = ApplyT(ComplexF64(khat), vhat)
        Tv_dense = Sdense * vhat
        op_err   = norm(Tv_fmm - Tv_dense) / max(norm(Tv_dense), 1e-300)
        @printf("  operator mismatch on vhat = %.3e\n", op_err)
    end

    return λ, V, ok
end

function test_circle_cavity_Dirichlet_Asakura_mFMM_EFIE(;
        a::Float64      = 1.0,
        p_order::Int    = 2,
        h::Float64      = 2π * a / (2 * 246),
        L::Int          = 8,
        Kblocks::Int    = 2,
        N::Int          = 64,
        δ_svd::Float64  = 1e-12,
        tol_res::Float64 = 1e-5,
        verbose::Bool   = true)

    println("=== Asakura eigenvalue test: circular cavity (Dirichlet, interior, EFIE, mFMM) ===")

    circle = CompScienceMeshes.meshcircle(a, h)
    X      = lagrangecx(circle, order=p_order)

    ctx = CircleEFIECtx(X; minvalues=100)

    D = length(X.fns)
    @printf("D = %d\n", D)

    k_true = j01 / a
    @printf("Analytic k_{0,1} ≈ %.10f\n", k_true)

    z0 = k_true + 0im
    R  = 0.10

    aca_tol = 1e-6
    η       = 1.0
    nearqs  = BEAST.DoubleNumSauterQstrat(30,30,0,4,25,25)
    farqs   = BEAST.DoubleNumQStrat(7,8)
    sched   = DynamicScheduler()

    UserFun = (k, Vhat, MIV) -> Beyn_circle_D(
        k, Vhat, MIV, ctx;
        rtol          = 1e-6,
        itmax         = 800,
        verbose       = 0,
        aca_tol       = aca_tol,
        η             = η,
        nearquadstrat = nearqs,
        farquadstrat  = farqs,
        scheduler     = sched
    )

    ApplyT = (k, V) -> ApplyT_num(
        ComplexF64(k), V, ctx;
        aca_tol       = aca_tol,
        η             = η,
        nearquadstrat = nearqs,
        farquadstrat  = farqs,
        scheduler     = sched
    )

    @time λ, V, ok = AsakuraSolve_circle_auto(
        z0, R,
        UserFun, ApplyT, D;
        L        = L,
        Kblocks  = Kblocks,
        N        = N,
        δ_svd    = δ_svd,
        tol_res  = tol_res,
        verbose  = verbose
    )

    println("ok = ", ok)

    p = sortperm(λ, by = x -> (real(x), imag(x)))
    λ = λ[p]
    V = V[:, p]

    println("Eigenvalues k inside contour (Asakura, mFMM):")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e, %+.12e   |k-ktrue|=%.3e\n",
                j, real(val), imag(val), abs(val-k_true))
    end

    if !isempty(λ)
        jbest = argmin(abs.(λ .- k_true))
        khat  = λ[jbest]
        vhat  = V[:, jbest]

        ApplyT_best = (k, V) -> ApplyT_num(
            ComplexF64(k), V, ctx;
            aca_tol       = aca_tol,
            η             = η,
            nearquadstrat = nearqs,
            farquadstrat  = farqs,
            scheduler     = sched
        )

        rfmm = norm(ApplyT_best(ComplexF64(khat), vhat)) / norm(vhat)

        Sop = Helmholtz2D.singlelayer(; wavenumber = ComplexF64(khat))
        nearquadstrat_dense = BEAST.DoubleNumSauterQstrat(8, 9, 0, 4, 9, 10)
        Sdense = assemble(Sop, X, X; quadstrat = nearquadstrat_dense)
        rdense = norm(Sdense * vhat) / norm(vhat)

        @printf("Residual(best): khat=%+.6f%+.6fi |khat-ktrue|=%.3e rfmm=%.3e rdense=%.3e\n",
                real(khat), imag(khat), abs(khat-k_true), rfmm, rdense)

        ratio = rfmm / max(rdense, 1e-300)
        @printf("  residual ratio rfmm/rdense = %.3e\n", ratio)

        Tv_fmm   = ApplyT_best(ComplexF64(khat), vhat)
        Tv_dense = Sdense * vhat
        op_err   = norm(Tv_fmm - Tv_dense) / max(norm(Tv_dense), 1e-300)
        @printf("  operator mismatch on vhat = %.3e\n", op_err)
    end

    return λ, V, ok
end

# ============================================================
# Optional runs
# ============================================================

# Dense reference:
# λ_dense, V_dense, ok_dense = test_circle_cavity_Dirichlet_Beyn_dense_EFIE()
# λA_dense, XA_dense, okA_dense = test_circle_cavity_Dirichlet_Asakura_dense_EFIE()

# ACA / mFMM reference:
λ_m, V_m, ok_m = test_circle_cavity_Dirichlet_Beyn_mFMM_EFIE()
# λA_m, XA_m, okA_m = test_circle_cavity_Dirichlet_Asakura_mFMM_EFIE()