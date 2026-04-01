using LinearAlgebra, Random, Printf
using BEAST, CompScienceMeshes
using SpecialFunctions
using FunctionZeros

import BEAST: kernelvals, KernelValsHelmholtz2D
import SpecialFunctions: hankelh2
import CompScienceMeshes: normal
import LinearAlgebra: norm, dot

# ============================================================
# Complex-k patch
# ============================================================

function kernelvals(biop::BEAST.HelmholtzOperator2D{T,K}, tgeo, bgeo) where {T,K<:Complex}
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

    nMom_max = 2*Kmax
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
        Km = K*m
        Kl = K*ℓ

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
# Dense EFIE single-layer
# ============================================================

function build_dense_S(X, k::ComplexF64;
                       quadstrat = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15))
    Sop = Helmholtz2D.singlelayer(; wavenumber = k)
    S   = assemble(Sop, X, X; quadstrat = quadstrat)
    return Matrix(S)
end

function Beyn_circle_D_dense(k::ComplexF64,
                             VHat::AbstractMatrix{ComplexF64},
                             MInvVHat::AbstractMatrix{ComplexF64},
                             X;
                             quadstrat = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15))
    S = build_dense_S(X, k; quadstrat=quadstrat)
    MInvVHat .= S \ VHat
    return MInvVHat
end

function ApplyT_circle_D_dense(k::ComplexF64,
                               v::AbstractVecOrMat{ComplexF64},
                               X;
                               quadstrat = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15))
    S = build_dense_S(X, k; quadstrat=quadstrat)
    return S * v
end

# ============================================================
# Dense test
# ============================================================

function test_circle_cavity_Dirichlet_Beyn_dense(;
        radius::Float64 = 1.0,
        Nel::Int = 512,
        z0::ComplexF64 = ComplexF64(2.404825557695773, 0.0),
        R::Float64 = 0.4,
        L::Int = 8,
        N::Int = 32,
        Kmax::Int = 2,
        tol_svd::Float64 = 1e-8,
        tol_res::Float64 = 1e-8,
        seed::Int = 1,
        verbose::Bool = true,
        quadstrat = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15))

    println("=== Beyn eigenvalue test: circular cavity (Dirichlet, interior, EFIE, dense) ===")

    h = 2π * radius / Nel
    Γ = CompScienceMeshes.meshcircle(radius, h)
    X = lagrangecx(Γ, order=2)

    D = length(X)
    k_true = besselj_zero(0, 1) / radius

    @printf("D = %d\n", D)
    @printf("Analytic k_(0,1) ≈ %.10f\n", k_true)
    @printf("Dense contour center z0 = %+.12f%+.12fi, R = %.4f\n", real(z0), imag(z0), R)

    UserFun = (k, Vhat, MIV) -> Beyn_circle_D_dense(k, Vhat, MIV, X; quadstrat=quadstrat)
    ApplyT  = (k, v)         -> ApplyT_circle_D_dense(k, v, X; quadstrat=quadstrat)

    @time λ, V, ok = BeynSolve_auto(z0, R, UserFun, ApplyT, D;
        L = L,
        N = N,
        Kmax = Kmax,
        tol_svd = tol_svd,
        tol_res = tol_res,
        seed = seed,
        verbose = verbose)

    println("ok = ", ok)

    p = sortperm(λ, by = x -> abs(x - k_true))
    λ = λ[p]

    println("Eigenvalues inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e %+.12ei   |k-k_true| ≈ %.3e\n",
                j, real(val), imag(val), abs(val - k_true))
    end

    return λ, V, ok
end