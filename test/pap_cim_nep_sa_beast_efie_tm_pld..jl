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

const Shape = BEAST.Shape
const LagrangeBasis = BEAST.LagrangeBasis

import LinearAlgebra: norm, dot
import BEAST: kernelvals, KernelValsHelmholtz2D
import SpecialFunctions: hankelh2
import CompScienceMeshes: normal

j01 = besselj_zero(0, 1)

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

    VHat = randn(D, L) .+ im * randn(D, L)

    nMom_max = 2*Kmax
    Ap = [zeros(ComplexF64, D, L) for _ in 1:nMom_max]

    nT          = nthreads()
    Ap_thread   = [ [zeros(ComplexF64, D, L) for _ in 1:nMom_max] for _ in 1:nT ]
    MInv_thread = [ zeros(ComplexF64, D, L) for _ in 1:nT ]

    Δθ = 2π / N

    @threads for nn in 0:(N-1)
        tid = threadid()
        Ap_loc   = Ap_thread[tid]
        MInvVHat = MInv_thread[tid]

        θ  = nn * Δθ
        cθ = cos(θ)
        sθ = sin(θ)

        w  = Rx * cθ + im * Ry * sθ
        z  = z0 + w
        dz = (-Rx * sθ + im * Ry * cθ) * Δθ

        MInvVHat = UserBeynFunction(z, VHat, MInvVHat)

        wp = one(w)
        @inbounds for p = 1:nMom_max
            Ap_loc[p] .+= wp * dz .* MInvVHat
            wp *= w
        end
    end

    for t in 1:nT
        Ap_loc = Ap_thread[t]
        @inbounds for p in 1:nMom_max
            Ap[p] .+= Ap_loc[p]
        end
    end

    function Beyn_from_moments(K::Int)
        nMom = 2*K
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
# Asakura / block-SS
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

    KL = Kblocks * L
    nMom = 2*Kblocks

    verbose && println("=== Asakura / block-SS contour solver (circle) ===")
    verbose && @printf("Center z0 = %+.6e%+.6ei, R = %.3f, D = %d, L = %d, Kblocks = %d, N = %d\n",
                       real(z0), imag(z0), R, D, L, Kblocks, N)

    V   = randn(ComplexF64, D, L)
    S_blocks = [zeros(ComplexF64, D, L) for _ in 1:nMom]
    MInv = zeros(ComplexF64, D, L)

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
    G        = Σ1 \ Hlt_tilde

    eigG = eigen(G)
    μ    = eigG.values
    Y    = eigG.vectors

    λ_hat = z0 .+ μ

    S_matrix = hcat(S_blocks[1:Kblocks]...)
    W_orig = VH1 * Y

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
    ok = true

    return λ_good, X_good, ok
end

# ============================================================
# Patch for complex k in 2D Helmholtz kernel
# ============================================================

function kernelvals(biop::BEAST.HelmholtzOperator2D{T, K}, tgeo, bgeo) where {T, K <: Complex}
    γ = biop.gamma
    k = iszero(real(γ)) ? imag(γ) : -im * γ

    r = tgeo.cart - bgeo.cart
    R = norm(r)

    kr = k * R
    H0 = hankelh2(0, kr)
    H1 = hankelh2(1, kr)

    green     = -im/4 * H0
    gradgreen =  k*im/4 * H1 * (r / R)

    txty = dot(normal(tgeo), normal(bgeo))
    return KernelValsHelmholtz2D(γ, r, R, green, gradgreen, txty)
end

# ============================================================
# EFIE / single-layer operator for TM interior Dirichlet
# ============================================================

"""
    assemble_T_circle_Dirichlet(k, X)

Dense EFIE operator for 2D TM interior Dirichlet cavity:

    T(k) = S(k)

where S is the Helmholtz single-layer boundary operator.
"""
function assemble_T_circle_Dirichlet(k::Number, X)
    𝒮 = Helmholtz2D.singlelayer(; wavenumber = k)
    S = assemble(𝒮, X, X)
    return Matrix(S)
end

"""
    Beyn_circle_D(k, VHat, MInvVHat, X)

User function for contour solvers:
solve S(k) * Y = VHat.
"""
function Beyn_circle_D(k, VHat, MInvVHat, X)
    T = assemble_T_circle_Dirichlet(k, X)
    MInvVHat .= T \ VHat
    return MInvVHat
end

"""
    ApplyT_circle_D(k, v, X)

Apply S(k) to vector or block vector v.
"""
function ApplyT_circle_D(k::ComplexF64,
                         v::AbstractVecOrMat{ComplexF64},
                         X)
    T = assemble_T_circle_Dirichlet(k, X)
    return T * v
end

# ============================================================
# Optional diagnostic: sigma_min scan
# ============================================================

function scan_sigma_min_circle_EFIE(; a=1.0, p_order=2, h=2π*a/64,
                                     kmin=2.0, kmax=2.8, nk=60)

    println("=== sigma_min scan: circular cavity, dense EFIE ===")

    circle = CompScienceMeshes.meshcircle(a, h)
    X = lagrangecx(circle, order=p_order)

    ks = collect(range(kmin, kmax, length=nk))
    σmins = zeros(Float64, nk)

    for (i, k) in enumerate(ks)
        T = assemble_T_circle_Dirichlet(complex(k), X)
        σmins[i] = minimum(svdvals(T))
        @printf("k = %.8f   sigma_min = %.6e\n", k, σmins[i])
    end

    imin = argmin(σmins)
    @printf("Best k from sigma_min scan: %.10f   sigma_min = %.6e\n", ks[imin], σmins[imin])

    return ks, σmins
end

# ============================================================
# Beyn test
# ============================================================

function test_circle_cavity_Dirichlet_Beyn()
    println("=== Beyn eigenvalue test: circular cavity (Dirichlet, interior, EFIE) ===")

    a       = 1.0
    p_order = 2
    h       = 2π * a / 64

    circle = CompScienceMeshes.meshcircle(a, h)
    X = lagrangecx(circle, order=p_order)

    T0 = assemble_T_circle_Dirichlet(2.0 + 0im, X)
    D  = size(T0, 1)
    @printf("D = %d\n", D)

    j01 = besselj_zero(0, 1)
    k_true = j01 / a
    @printf("Analytic k_{0,1} ≈ %.10f\n", k_true)

    z0 = k_true + 0.0im
    R  = 0.4

    L       = 8
    N       = 32
    Kmax    = 2
    tol_svd = 1e-8
    tol_res = 1e-8
    verbose = true

    UserFun = (k, Vhat, MIV) -> Beyn_circle_D(k, Vhat, MIV, X)
    ApplyT  = (k, v)         -> ApplyT_circle_D(k, v, X)

    @time λ, V, ok = BeynSolve_auto(z0, R,
                                    UserFun, ApplyT, D;
                                    L       = L,
                                    N       = N,
                                    Kmax    = Kmax,
                                    tol_svd = tol_svd,
                                    tol_res = tol_res,
                                    verbose = verbose)

    println("ok = ", ok)
    if !ok
        println("Warning: BeynSolve_auto did not find a fully verified set of eigenpairs.")
    end

    p = sortperm(λ, by = x -> (real(x), imag(x)))
    λ = λ[p]

    println("Eigenvalues k inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e, %+.12e   |k - k_true| ≈ %.3e\n",
                j, real(val), imag(val), abs(val - k_true))
    end

    return λ, V, ok
end

# ============================================================
# Asakura test
# ============================================================

function test_circle_cavity_Dirichlet_Asakura(;
        a::Float64      = 1.0,
        p_order::Int    = 2,
        h::Float64      = 2π * a / 64,
        L::Int          = 8,
        Kblocks::Int    = 2,
        N::Int          = 64,
        δ_svd::Float64  = 1e-12,
        tol_res::Float64 = 1e-8,
        verbose::Bool   = true)

    println("=== Asakura eigenvalue test: circular cavity (Dirichlet, interior, EFIE) ===")

    circle = CompScienceMeshes.meshcircle(a, h)
    X = lagrangecx(circle, order=p_order)

    T0 = assemble_T_circle_Dirichlet(2.0 + 0im, X)
    D  = size(T0, 1)
    @printf("D = %d\n", D)

    j01    = besselj_zero(0, 1)
    k_true = j01 / a
    @printf("Analytic k_{0,1} ≈ %.10f\n", k_true)

    z0 = k_true + 0.0im
    R  = 0.4

    UserFun = (k, Vhat, MInvVHat) -> Beyn_circle_D(k, Vhat, MInvVHat, X)
    ApplyT  = (k, v)              -> ApplyT_circle_D(k, v, X)

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

    println("Eigenvalues k inside contour (Asakura):")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e, %+.12e   |k - k_true| ≈ %.3e\n",
                j, real(val), imag(val), abs(val - k_true))
    end

    return λ, V, ok
end

# ============================================================
# Run
# ============================================================

# 1) optional: first do sigma_min scan
# ks, σmins = scan_sigma_min_circle_EFIE()

# 2) Beyn
λ, V, ok = test_circle_cavity_Dirichlet_Beyn()

# 3) Asakura
λA, XA, okA = test_circle_cavity_Dirichlet_Asakura()