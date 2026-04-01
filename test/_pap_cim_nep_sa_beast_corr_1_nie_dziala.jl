using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test
using SpecialFunctions
using Gmsh

import LinearAlgebra: norm
const Shape = BEAST.Shape
const LagrangeBasis = BEAST.LagrangeBasis

using Printf
using Base.Threads
using FunctionZeros
j01 = besselj_zero(0, 1)

# ============================================================
# BeynSolve_auto (bez zmian)
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
        cθ = cos(θ); sθ = sin(θ)

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
        verbose && @printf("K=%d: rank(B0) ≈ k = %d\n", K, k)

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
            verbose && @printf("  Residuals (K=%d):\n", K)
            for j in 1:length(λ)
                vj = Vfull[:, j]
                Tv = ApplyT(λ[j], vj)
                rj = norm(Tv) / norm(vj)
                r[j] = rj
                verbose && @printf("    res[%2d] = %.3e\n", j, rj)
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

        verbose && @printf("  K=%d: %d good eigenpairs (res ≤ %.1e)\n",
                           K, ngood, tol_res)

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
            verbose && @printf("  #good eigenpairs stabilized at %d for K=%d.\n", ngood, K)
            λ_best = λg
            V_best = Vg
            K_used = K
            break
        end

        good_prev = ngood
    end

    if good_best == 0 || isempty(λ_best)
        verbose && @printf("No acceptable eigenpairs found.\n")
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
# kernelvals override (bez zmian)
# ============================================================
import BEAST: kernelvals, KernelValsHelmholtz2D
import SpecialFunctions: hankelh2
import CompScienceMeshes: normal
import LinearAlgebra: norm, dot

function kernelvals(biop::BEAST.HelmholtzOperator2D{T, K}, tgeo, bgeo) where {T, K <: Complex}
    γ = biop.gamma
    k = iszero(real(γ)) ? imag(γ) : -im*γ

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
# DENSE: SINGLE-LAYER operator  T(k) = S(k)
# (tu była jedyna zmiana: singlelayer + brak -0.5 I)
# ============================================================
"""
    assemble_T_circle_Dirichlet(k, X)

TU TERAZ: T(k) = S(k) (single-layer).
Zostawiam nazwę funkcji, żeby nic innego nie ruszać.
"""
function assemble_T_circle_Dirichlet(k::Number, X)
    𝒮 = Helmholtz2D.singlelayer(; wavenumber = k)  # <<< CHANGED
    S  = assemble(𝒮, X, X)
    I  = assemble(BEAST.Identity(), X, X)          # zostawiam (0*I) żeby reszta kodu się nie zmieniała
    return Matrix(0.0 .* I .+ S)                   # <<< CHANGED (α=0)
end

function Beyn_circle_D(k, VHat, MInvVHat, X)
    T = assemble_T_circle_Dirichlet(k, X)
    MInvVHat .= T \ VHat
    return MInvVHat
end

function ApplyT_circle_D(k::ComplexF64, v::AbstractVecOrMat{ComplexF64}, X)
    T = assemble_T_circle_Dirichlet(k, X)
    return T * v
end

function test_circle_cavity_Dirichlet_Beyn()
    println("=== Beyn eigenvalue test: circular cavity (SINGLE-LAYER) ===")

    a       = 1.0
    p_order = 2
    h       = 2π * a / 64

    circle = CompScienceMeshes.meshcircle(a, h)
    X = lagrangecx(circle, order=p_order)

    T0 = assemble_T_circle_Dirichlet(2.0, X)
    D  = size(T0, 1)
    @printf("D = %d\n", D)

    j01 = besselj_zero(0, 1)
    k_true = j01 / a
    @printf("Reference j01/a ≈ %.10f\n", k_true)

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

    @time λ, V, ok = BeynSolve_auto(z0, R, UserFun, ApplyT, D;
                                   L=L, N=N, Kmax=Kmax,
                                   tol_svd=tol_svd, tol_res=tol_res,
                                   verbose=verbose)

    println("ok = ", ok)
    p = sortperm(λ, by = x -> (real(x), imag(x)))
    λ = λ[p]

    println("Eigenvalues k inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e, %+.12e   |k-kref| ≈ %.3e\n",
                j, real(val), imag(val), abs(val - k_true))
    end
    return λ, V, ok
end

# ============================================================
# mFMM / ACA / H2
# ============================================================
using ParallelKMeans
using H2Trees
using AdaptiveCrossApproximation
using OhMyThreads
using Krylov

import LinearAlgebra: mul!

struct DirichletTOp{K,I,T}
    K::K
    I::I
    α::T
end
Base.size(A::DirichletTOp) = size(A.I)

function mul!(y::AbstractVector, A::DirichletTOp, x::AbstractVector)
    mul!(y, A.K, x)
    mul!(y, A.I, x, A.α, one(A.α))
    return y
end

# --- BlockTree (zostawiam jak u Ciebie: KMeansTree(X.pos,...) ) ---
function build_blocktree(X; minvalues=100)
    testtree  = KMeansTree(X.pos, 2; minvalues=minvalues)
    trialtree = KMeansTree(X.pos, 2; minvalues=minvalues)
    return BlockTree(testtree, trialtree)
end

function hassemble_with_tree(op, X, tree;
    aca_tol::Float64 = 1e-4,
    η::Float64 = 1.0,
    nearquadstrat = BEAST.DoubleNumSauterQstrat(8,8,0,4,7,7),
    farquadstrat  = BEAST.DoubleNumQStrat(4,5),
    scheduler     = DynamicScheduler(),
)
    @show aca_tol η
    return AdaptiveCrossApproximation.HMatrix(
        op, X, X, tree;
        compressor    = ACA(; tol = aca_tol),
        isnear        = AdaptiveCrossApproximation.isnear(; η = η),
        nearquadstrat = nearquadstrat,
        farquadstrat  = farquadstrat,
        scheduler     = scheduler,
    )
end

mutable struct CircleDirichletCtx
    X
    tree
    I
    Kcache::Dict{ComplexF64, Any}
end

function CircleDirichletCtx(X; minvalues=100)
    tree = build_blocktree(X; minvalues=minvalues)
    I    = assemble(BEAST.Identity(), X, X)    # zostaje (ale α=0)
    return CircleDirichletCtx(X, tree, I, Dict{ComplexF64,Any}())
end

# Get/build "K(k)" — TU TERAZ jest SINGLELAYER (S(k))
function getK!(ctx::CircleDirichletCtx, k::ComplexF64)
    get!(ctx.Kcache, k) do
        𝒮 = Helmholtz2D.singlelayer(; wavenumber=k)  # <<< CHANGED
        @time hassemble_with_tree(𝒮, ctx.X, ctx.tree; aca_tol=1e-6, η=10.0)
    end
end

function Beyn_circle_D(k, VHat, MInvVHat, ctx; rtol=1e-12, itmax=8000, verbose=0)
    k  = ComplexF64(k)
    Kk = getK!(ctx, k)

    α0 = ComplexF64(0.0)                         # <<< CHANGED (było -0.5)
    Top = DirichletTOp(Kk, ctx.I, α0)

    n, L = size(VHat)
    @assert size(MInvVHat) == (n, L)

    for j in 1:L
        x, st = Krylov.gmres(Top, view(VHat, :, j); rtol=rtol, itmax=itmax, verbose=verbose)
        st.solved || @warn "GMRES not solved" k=k j=j status=st.status niter=st.niter
        MInvVHat[:, j] .= x
    end
    return MInvVHat
end

function test_circle_cavity_Dirichlet_Beyn_mFMM()
    println("=== Beyn eigenvalue test: circular cavity (SINGLE-LAYER) ===_mFMM")

    a       = 1.0
    p_order = 2
    h       = 2π * a / (2*246)

    circle = CompScienceMeshes.meshcircle(a, h)
    X = lagrangecx(circle, order=p_order)

    ctx = CircleDirichletCtx(X; minvalues=100)
    D = size(ctx.I, 1)
    @printf("D = %d\n", D)

    j01 = besselj_zero(0, 1)
    k_true = j01 / a
    @printf("Reference j01/a ≈ %.10f\n", k_true)

    z0 = k_true + 0.0im
    R  = 0.10

    L       = 12
    N       = 64
    Kmax    = 2
    tol_svd = 1e-8
    tol_res = 1e-3
    verbose = true

    UserFun = (k, Vhat, MIV) -> Beyn_circle_D(k, Vhat, MIV, ctx; rtol=1e-8, itmax=400)

    function ApplyT_num(k::ComplexF64, V)
        Kk = getK!(ctx, k)
        I  = ctx.I
        α  = ComplexF64(0.0)         # <<< CHANGED (było -0.5)
        β  = ComplexF64(1.0)

        if V isa AbstractVector
            y = similar(V, ComplexF64)
            mul!(y, Kk, V)
            mul!(y, I, V, α, β)      # dodaje 0*I*V (czyli nic)
            return y
        else
            @assert size(V,1) == size(I,1)
            Y = similar(V, ComplexF64)
            for j in axes(V,2)
                vj = view(V, :, j)
                yj = view(Y, :, j)
                mul!(yj, Kk, vj)
                mul!(yj, I, vj, α, β)
            end
            return Y
        end
    end
    ApplyT = (k, V) -> ApplyT_num(ComplexF64(k), V)

    # Consistency check
    kcheck = ComplexF64(z0 + R*exp(im*0.37))
    Ltest  = 3
    Vhat   = randn(ComplexF64, D, Ltest)
    MIV    = similar(Vhat)
    UserFun(kcheck, Vhat, MIV)
    Tmiv = ApplyT(kcheck, MIV)
    err  = norm(Tmiv - Vhat) / norm(Vhat)
    @printf("Consistency check: ||T*MIV - Vhat||/||Vhat|| = %.3e\n", err)

    @time λ, V, ok = BeynSolve_auto(z0, R, UserFun, ApplyT, D;
                                   L=L, N=N, Kmax=Kmax,
                                   tol_svd=tol_svd, tol_res=tol_res,
                                   verbose=verbose)

    println("ok = ", ok)

    p = sortperm(λ, by = x -> (real(x), imag(x)))
    λ = λ[p]; V = V[:, p]

    println("Eigenvalues k inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e, %+.12e   |k-kref|=%.3e\n",
                j, real(val), imag(val), abs(val-k_true))
    end

    if !isempty(λ)
        jbest = argmin(abs.(λ .- k_true))
        khat  = λ[jbest]
        vhat  = V[:, jbest]

        Tfmm = DirichletTOp(getK!(ctx, ComplexF64(khat)), ctx.I, ComplexF64(0.0)) # <<< CHANGED α=0
        rfmm = norm(ApplyT(ComplexF64(khat), vhat)) / norm(vhat)

        Dop = Helmholtz2D.singlelayer(; wavenumber = ComplexF64(khat))           # <<< CHANGED
        nearquadstrat_dense = BEAST.DoubleNumSauterQstrat(8, 9, 0, 4, 9, 10)
        Kdense = assemble(Dop, X, X; quadstrat = nearquadstrat_dense)
        Tdense = Matrix(Kdense)                                                  # <<< CHANGED (bez -0.5 I)
        rdense = norm(Tdense * vhat) / norm(vhat)

        @printf("Residual(best): khat=%+.6f%+.6fi |khat-kref|=%.3e rfmm=%.3e rdense=%.3e\n",
                real(khat), imag(khat), abs(khat-k_true), rfmm, rdense)

        Tv_fmm   = ApplyT(ComplexF64(khat), vhat)
        Tv_dense = Tdense * vhat
        op_err   = norm(Tv_fmm - Tv_dense) / max(norm(Tv_dense), 1e-300)
        @printf("  operator mismatch on vhat = %.3e\n", op_err)
    end

    return λ, V, ok
end

λ_m, V_m, ok_m = test_circle_cavity_Dirichlet_Beyn_mFMM()

# (reszta, np. Asakura, możesz zostawić/wyciąć — nie ruszałem tutaj, bo prosiłeś “nic więcej”)