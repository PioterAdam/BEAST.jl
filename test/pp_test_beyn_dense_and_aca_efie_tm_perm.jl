
# pp_test_beyn_dense_and_aca_efie_tm_perm.jl
println("RUNNING FILE: pp_test_beyn_dense_and_aca_efie_tm_perm.jl")

using LinearAlgebra, Random, Statistics, Printf
using BEAST, CompScienceMeshes
using AdaptiveCrossApproximation
using OhMyThreads
using Krylov
using H2Trees
using SpecialFunctions
using LinearMaps
using FunctionZeros

import BEAST: kernelvals, KernelValsHelmholtz2D
import SpecialFunctions: hankelh2
import CompScienceMeshes: normal
import LinearAlgebra: norm, dot
import AdaptiveCrossApproximation: permute

import H2Trees: permutation

function _tree_permutations(ctx::CircleEFIECtx)
    testtree  = getfield(ctx.tree, 1)
    trialtree = getfield(ctx.tree, 2)

    ptest  = collect(Int, H2Trees.permutation(testtree))
    ptrial = collect(Int, H2Trees.permutation(trialtree))

    return ptest, ptrial
end

function permute(space::BEAST.LagrangeBasis, p::AbstractVector{<:Integer})
    pI = collect(Int, p)
    mesh = getfield(space, 1)
    fns  = getfield(space, 2)
    pos  = getfield(space, 3)
    return typeof(space)(mesh, fns[pI], pos[pI])
end

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
    X = BEAST.lagrangecx(Γ, order=2)
    D = length(X)
    k_true = besselj_zero(0, 1) / radius

    @printf("D = %d\n", D)
    @printf("Analytic k_(0,1) ≈ %.10f\n", k_true)
    @printf("Dense contour center z0 = %+.12f%+.12fi, R = %.4f\n", real(z0), imag(z0), R)

    UserFun = (k, Vhat, MIV) -> Beyn_circle_D_dense(k, Vhat, MIV, X; quadstrat=quadstrat)
    ApplyT  = (k, v)         -> ApplyT_circle_D_dense(k, v, X; quadstrat=quadstrat)

    @time λ, V, ok = BeynSolve_auto(z0, R, UserFun, ApplyT, D;
        L = L, N = N, Kmax = Kmax, tol_svd = tol_svd, tol_res = tol_res, seed = seed, verbose = verbose)

    println("ok = ", ok)

    p = sortperm(λ, by = x -> abs(x - k_true))
    λ = λ[p]

    println("Eigenvalues inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e %+.12ei   |k-k_true| ≈ %.3e\n", j, real(val), imag(val), abs(val - k_true))
    end

    return λ, V, ok
end

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
    X = BEAST.lagrangecx(Γ, order=2)
    tree = build_blocktree(X; minvalues=minvalues)
    return CircleEFIECtx(X, tree)
end

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

function _tree_permutations(ctx::CircleEFIECtx)
    testtree  = getfield(ctx.tree, 1)
    trialtree = getfield(ctx.tree, 2)
    ptest  = collect(Int, permutation(testtree))
    ptrial = collect(Int, permutation(trialtree))
    return ptest, ptrial
end

function _make_permuted_hmatrix(H, perms::Tuple)
    try
        return AdaptiveCrossApproximation.PermutedHMatrix(H, perms)
    catch
        try
            return AdaptiveCrossApproximation.PermutedHMatrix(perms, H)
        catch err
            error("Could not construct PermutedHMatrix from HMatrix and permutations: $(sprint(showerror, err))")
        end
    end
end

function build_fmm_S_perm(ctx::CircleEFIECtx, k::ComplexF64;
                          tol::Float64 = 1e-4,
                          η::Float64   = 0.25,
                          scheduler     = SerialScheduler(),
                          nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                          farquadstrat  = BEAST.DoubleNumQStrat(5,4))
    Sop = Helmholtz2D.singlelayer(; wavenumber = k)
    ptest, ptrial = _tree_permutations(ctx)
    Xtest_p  = permute(ctx.X, ptest)
    Xtrial_p = permute(ctx.X, ptrial)

    H = AdaptiveCrossApproximation.HMatrix(
        Sop, Xtest_p, Xtrial_p, ctx.tree;
        compressor    = ACA(; convergence = FNormEstimator(tol)),
        isnear        = AdaptiveCrossApproximation.isnear(; η = η),
        nearquadstrat = nearquadstrat,
        farquadstrat  = farquadstrat,
        scheduler     = scheduler,
        perm          = false,
    )

    return _make_permuted_hmatrix(H, (ptest, ptrial))
end

function Beyn_circle_D_ACA(k::ComplexF64, VHat::AbstractMatrix{ComplexF64}, MInvVHat::AbstractMatrix{ComplexF64}, ctx::CircleEFIECtx;
                           tol::Float64 = 1e-4, η::Float64 = 0.25, scheduler = SerialScheduler(),
                           nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                           farquadstrat  = BEAST.DoubleNumQStrat(5,4),
                           gmres_rtol::Float64 = 1e-10, gmres_itmax::Int = 1500, verbose::Bool = false)
    S = build_fmm_S(ctx, k; tol = tol, η = η, scheduler = scheduler, nearquadstrat = nearquadstrat, farquadstrat = farquadstrat)
    nrhs = size(VHat, 2)
    for j in 1:nrhs
        b = view(VHat, :, j)
        x, st = Krylov.gmres(S, b; rtol=gmres_rtol, itmax=gmres_itmax, verbose=0)
        if verbose && !st.solved
            @printf("GMRES warning at k=%+.6e%+.6ei, rhs=%d, status=%s, niter=%d\n", real(k), imag(k), j, string(st.status), st.niter)
        end
        MInvVHat[:, j] .= x
    end
    return MInvVHat
end

function ApplyT_circle_D_ACA(k::ComplexF64, v::AbstractVecOrMat{ComplexF64}, ctx::CircleEFIECtx;
                             tol::Float64 = 1e-4, η::Float64 = 0.25, scheduler = SerialScheduler(),
                             nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                             farquadstrat  = BEAST.DoubleNumQStrat(5,4))
    S = build_fmm_S(ctx, k; tol = tol, η = η, scheduler = scheduler, nearquadstrat = nearquadstrat, farquadstrat = farquadstrat)
    return S * v
end

function Beyn_circle_D_ACA_perm(k::ComplexF64, VHat::AbstractMatrix{ComplexF64}, MInvVHat::AbstractMatrix{ComplexF64}, ctx::CircleEFIECtx;
                                tol::Float64 = 1e-4, η::Float64 = 0.25, scheduler = SerialScheduler(),
                                nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                                farquadstrat  = BEAST.DoubleNumQStrat(5,4),
                                gmres_rtol::Float64 = 1e-10, gmres_itmax::Int = 1500, verbose::Bool = false)
    S = build_fmm_S_perm(ctx, k; tol=tol, η=η, scheduler=scheduler, nearquadstrat=nearquadstrat, farquadstrat=farquadstrat)
    nrhs = size(VHat, 2)
    for j in 1:nrhs
        b = view(VHat, :, j)
        x, st = Krylov.gmres(S, b; rtol=gmres_rtol, itmax=gmres_itmax, verbose=0)
        if verbose && !st.solved
            @printf("GMRES warning at k=%+.6e%+.6ei, rhs=%d, status=%s, niter=%d\n", real(k), imag(k), j, string(st.status), st.niter)
        end
        MInvVHat[:, j] .= x
    end
    return MInvVHat
end

function ApplyT_circle_D_ACA_perm(k::ComplexF64, v::AbstractVecOrMat{ComplexF64}, ctx::CircleEFIECtx;
                                  tol::Float64 = 1e-4, η::Float64 = 0.25, scheduler = SerialScheduler(),
                                  nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                                  farquadstrat  = BEAST.DoubleNumQStrat(5,4))
    S = build_fmm_S_perm(ctx, k; tol=tol, η=η, scheduler=scheduler, nearquadstrat=nearquadstrat, farquadstrat=farquadstrat)
    return S * v
end

function contour_preflight_aca(ctx::CircleEFIECtx; z0::ComplexF64, Rx::Float64, Ry::Float64, N::Int = 12,
                               tol::Float64 = 1e-4, η::Float64 = 0.25, scheduler = SerialScheduler(),
                               nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                               farquadstrat  = BEAST.DoubleNumQStrat(5,4),
                               seed::Int = 1, use_perm::Bool = false)
    println(use_perm ? "=== ACA perm=true contour preflight ===" : "=== ACA contour preflight ===")
    @printf("z0 = %+.12f%+.12fi, Rx = %.4f, Ry = %.4f, N = %d\n", real(z0), imag(z0), Rx, Ry, N)

    Random.seed!(seed)
    D = length(ctx.X)
    Δθ = 2π / N

    for nn in 0:(N-1)
        θ = nn * Δθ
        z = z0 + Rx*cos(θ) + im*Ry*sin(θ)
        @printf("point %2d: k = %+.12f%+.12fi ... ", nn+1, real(z), imag(z))
        try
            S = use_perm ?
                build_fmm_S_perm(ctx, z; tol=tol, η=η, scheduler=scheduler, nearquadstrat=nearquadstrat, farquadstrat=farquadstrat) :
                build_fmm_S(ctx, z; tol=tol, η=η, scheduler=scheduler, nearquadstrat=nearquadstrat, farquadstrat=farquadstrat)
            v = randn(ComplexF64, D)
            y = S * v
            if all(isfinite, y)
                println("OK")
            else
                println("NON-FINITE matvec")
                return false, z
            end
        catch err
            println("FAIL")
            @show err
            @printf("failed at k = %+.12f%+.12fi\n", real(z), imag(z))
            return false, z
        end
    end

    println("all contour points passed")
    return true, z0
end

function test_circle_cavity_Dirichlet_Beyn_ACA_circle(; radius::Float64 = 1.0, Nel::Int = 512, minvalues::Int = 100,
        tol::Float64 = 1e-4, η::Float64 = 0.25, scheduler = SerialScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9), farquadstrat  = BEAST.DoubleNumQStrat(5,4),
        z0::ComplexF64 = ComplexF64(2.404825557695773, 0.0), R::Float64 = 0.05, L::Int = 8, N::Int = 32, Kmax::Int = 2,
        tol_svd::Float64 = 1e-8, tol_res::Float64 = 1e-5, gmres_rtol::Float64 = 1e-10, gmres_itmax::Int = 1500,
        seed::Int = 1, verbose::Bool = true)
    println("=== Beyn eigenvalue test: circular cavity (Dirichlet, interior, EFIE, ACA circle) ===")
    ctx = make_ctx_circle(; radius=radius, Nel=Nel, minvalues=minvalues)
    D = length(ctx.X)
    k_true = besselj_zero(0, 1) / radius
    @printf("D = %d\n", D)
    @printf("Analytic k_(0,1) ≈ %.10f\n", k_true)
    @printf("ACA params: tol = %.1e, η = %.3f, scheduler = %s\n", tol, η, string(typeof(scheduler)))
    UserFun = (k, Vhat, MIV) -> Beyn_circle_D_ACA(k, Vhat, MIV, ctx;
        tol=tol, η=η, scheduler=scheduler, nearquadstrat=nearquadstrat, farquadstrat=farquadstrat,
        gmres_rtol=gmres_rtol, gmres_itmax=gmres_itmax, verbose=false)
    ApplyT = (k, v) -> ApplyT_circle_D_ACA(k, v, ctx;
        tol=tol, η=η, scheduler=scheduler, nearquadstrat=nearquadstrat, farquadstrat=farquadstrat)
    @time λ, V, ok = BeynSolve_auto(z0, R, UserFun, ApplyT, D;
        L=L, N=N, Kmax=Kmax, tol_svd=tol_svd, tol_res=tol_res, seed=seed, verbose=verbose)
    println("ok = ", ok)
    p = sortperm(λ, by = x -> abs(x - k_true))
    λ = λ[p]
    println("Eigenvalues inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e %+.12ei   |k-k_true| ≈ %.3e\n", j, real(val), imag(val), abs(val - k_true))
    end
    return λ, V, ok
end

function test_circle_cavity_Dirichlet_Beyn_ACA_ellipse(; radius::Float64 = 1.0, Nel::Int = 512, minvalues::Int = 100,
        tol::Float64 = 1e-4, η::Float64 = 0.25, scheduler = SerialScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9), farquadstrat  = BEAST.DoubleNumQStrat(5,4),
        z0::ComplexF64 = ComplexF64(2.4067799554075626, 0.0), Rx::Float64 = 0.05, Ry::Float64 = 0.02, L::Int = 8, N::Int = 32, Kmax::Int = 2,
        tol_svd::Float64 = 1e-8, tol_res::Float64 = 1e-5, gmres_rtol::Float64 = 1e-10, gmres_itmax::Int = 1500,
        seed::Int = 1, verbose::Bool = true, do_preflight::Bool = true)
    println("=== Beyn eigenvalue test: circular cavity (Dirichlet, interior, EFIE, ACA ellipse) ===")
    ctx = make_ctx_circle(; radius=radius, Nel=Nel, minvalues=minvalues)
    D = length(ctx.X)
    k_true = besselj_zero(0, 1) / radius
    @printf("D = %d\n", D)
    @printf("Analytic k_(0,1) ≈ %.10f\n", k_true)
    @printf("ACA params: tol = %.1e, η = %.3f, scheduler = %s\n", tol, η, string(typeof(scheduler)))
    if do_preflight
        ok_pre, _ = contour_preflight_aca(ctx; z0 = z0, Rx = Rx, Ry = Ry, N = min(N, 12), tol = tol, η = η, scheduler = scheduler,
            nearquadstrat = nearquadstrat, farquadstrat = farquadstrat, seed = seed, use_perm = false)
        if !ok_pre
            println("Preflight failed.")
            return ComplexF64[], zeros(ComplexF64, D, 0), false
        end
    end
    UserFun = (k, Vhat, MIV) -> Beyn_circle_D_ACA(k, Vhat, MIV, ctx;
        tol=tol, η=η, scheduler=scheduler, nearquadstrat=nearquadstrat, farquadstrat=farquadstrat,
        gmres_rtol=gmres_rtol, gmres_itmax=gmres_itmax, verbose=false)
    ApplyT = (k, v) -> ApplyT_circle_D_ACA(k, v, ctx;
        tol=tol, η=η, scheduler=scheduler, nearquadstrat=nearquadstrat, farquadstrat=farquadstrat)
    @time λ, V, ok = BeynSolve_auto(z0, Rx, Ry, UserFun, ApplyT, D;
        L=L, N=N, Kmax=Kmax, tol_svd=tol_svd, tol_res=tol_res, seed=seed, verbose=verbose)
    println("ok = ", ok)
    p = sortperm(λ, by = x -> abs(x - k_true))
    λ = λ[p]
    println("Eigenvalues inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e %+.12ei   |k-k_true| ≈ %.3e\n", j, real(val), imag(val), abs(val - k_true))
    end
    return λ, V, ok
end

function test_circle_cavity_Dirichlet_Beyn_ACA_perm_circle(; radius::Float64 = 1.0, Nel::Int = 512, minvalues::Int = 100,
        tol::Float64 = 1e-4, η::Float64 = 0.25, scheduler = SerialScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9), farquadstrat  = BEAST.DoubleNumQStrat(5,4),
        z0::ComplexF64 = ComplexF64(2.404825557695773, 0.0), R::Float64 = 0.05, L::Int = 8, N::Int = 32, Kmax::Int = 2,
        tol_svd::Float64 = 1e-8, tol_res::Float64 = 1e-5, gmres_rtol::Float64 = 1e-10, gmres_itmax::Int = 1500,
        seed::Int = 1, verbose::Bool = true)
    println("=== Beyn eigenvalue test: circular cavity (Dirichlet, interior, EFIE, ACA perm=true circle) ===")
    ctx = make_ctx_circle(; radius=radius, Nel=Nel, minvalues=minvalues)
    D = length(ctx.X)
    k_true = besselj_zero(0, 1) / radius
    @printf("D = %d\n", D)
    @printf("Analytic k_(0,1) ≈ %.10f\n", k_true)
    @printf("ACA perm=true params: tol = %.1e, η = %.3f, scheduler = %s\n", tol, η, string(typeof(scheduler)))
    UserFun = (k, Vhat, MIV) -> Beyn_circle_D_ACA_perm(k, Vhat, MIV, ctx;
        tol=tol, η=η, scheduler=scheduler, nearquadstrat=nearquadstrat, farquadstrat=farquadstrat,
        gmres_rtol=gmres_rtol, gmres_itmax=gmres_itmax, verbose=false)
    ApplyT = (k, v) -> ApplyT_circle_D_ACA_perm(k, v, ctx;
        tol=tol, η=η, scheduler=scheduler, nearquadstrat=nearquadstrat, farquadstrat=farquadstrat)
    @time λ, V, ok = BeynSolve_auto(z0, R, UserFun, ApplyT, D;
        L = L, N = N, Kmax = Kmax, tol_svd = tol_svd, tol_res = tol_res, seed = seed, verbose = verbose)
    println("ok = ", ok)
    p = sortperm(λ, by = x -> abs(x - k_true))
    λ = λ[p]
    println("Eigenvalues inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e %+.12ei   |k-k_true| ≈ %.3e\n", j, real(val), imag(val), abs(val - k_true))
    end
    return λ, V, ok
end

function test_circle_cavity_Dirichlet_Beyn_ACA_perm_ellipse(; radius::Float64 = 1.0, Nel::Int = 512, minvalues::Int = 100,
        tol::Float64 = 1e-4, η::Float64 = 0.25, scheduler = SerialScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9), farquadstrat  = BEAST.DoubleNumQStrat(5,4),
        z0::ComplexF64 = ComplexF64(2.4067799554075626, 0.0), Rx::Float64 = 0.05, Ry::Float64 = 0.02, L::Int = 8, N::Int = 32, Kmax::Int = 2,
        tol_svd::Float64 = 1e-8, tol_res::Float64 = 1e-5, gmres_rtol::Float64 = 1e-10, gmres_itmax::Int = 1500,
        seed::Int = 1, verbose::Bool = true, do_preflight::Bool = true)
    println("=== Beyn eigenvalue test: circular cavity (Dirichlet, interior, EFIE, ACA perm=true ellipse) ===")
    ctx = make_ctx_circle(; radius=radius, Nel=Nel, minvalues=minvalues)
    D = length(ctx.X)
    k_true = besselj_zero(0, 1) / radius
    @printf("D = %d\n", D)
    @printf("Analytic k_(0,1) ≈ %.10f\n", k_true)
    @printf("ACA perm=true params: tol = %.1e, η = %.3f, scheduler = %s\n", tol, η, string(typeof(scheduler)))
    if do_preflight
        ok_pre, _ = contour_preflight_aca(ctx; z0 = z0, Rx = Rx, Ry = Ry, N = min(N, 12), tol = tol, η = η, scheduler = scheduler,
            nearquadstrat = nearquadstrat, farquadstrat = farquadstrat, seed = seed, use_perm = true)
        if !ok_pre
            println("Preflight failed.")
            return ComplexF64[], zeros(ComplexF64, D, 0), false
        end
    end
    UserFun = (k, Vhat, MIV) -> Beyn_circle_D_ACA_perm(k, Vhat, MIV, ctx;
        tol=tol, η=η, scheduler=scheduler, nearquadstrat=nearquadstrat, farquadstrat=farquadstrat,
        gmres_rtol=gmres_rtol, gmres_itmax=gmres_itmax, verbose=false)
    ApplyT = (k, v) -> ApplyT_circle_D_ACA_perm(k, v, ctx;
        tol=tol, η=η, scheduler=scheduler, nearquadstrat=nearquadstrat, farquadstrat=farquadstrat)
    @time λ, V, ok = BeynSolve_auto(z0, Rx, Ry, UserFun, ApplyT, D;
        L = L, N = N, Kmax = Kmax, tol_svd = tol_svd, tol_res = tol_res, seed = seed, verbose = verbose)
    println("ok = ", ok)
    p = sortperm(λ, by = x -> abs(x - k_true))
    λ = λ[p]
    println("Eigenvalues inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e %+.12ei   |k-k_true| ≈ %.3e\n", j, real(val), imag(val), abs(val - k_true))
    end
    return λ, V, ok
end

function test_perm_wrapper_efie(ctx::CircleEFIECtx, k::ComplexF64;
                                tol::Float64 = 1e-4, η::Float64   = 0.25, scheduler = SerialScheduler(),
                                nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9), farquadstrat  = BEAST.DoubleNumQStrat(5,4),
                                q_dense = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15), nvec::Int = 7, seed::Int = 1)
    println("=== test_perm_wrapper_efie ===")
    Sd = build_dense_S(ctx, k; quadstrat=q_dense)
    Sp = build_fmm_S_perm(ctx, k; tol=tol, η=η, scheduler=scheduler, nearquadstrat=nearquadstrat, farquadstrat=farquadstrat)
    @show typeof(Sp)
    if !hasproperty(Sp, :permutation) || !hasproperty(Sp, :M)
        error("Expected PermutedHMatrix with fields :permutation and :M")
    end
    ptest, ptrial = Sp.permutation
    Sd_perm = Sd[ptest, ptrial]
    Random.seed!(seed)
    e_int   = zeros(Float64, nvec)
    e_wrap  = zeros(Float64, nvec)
    e_dense = zeros(Float64, nvec)
    e_apply = zeros(Float64, nvec)
    for i in 1:nvec
        v = randn(ComplexF64, size(Sd, 2))
        xp = similar(v)
        xp .= v[ptrial]
        yh = Sp.M * xp
        yd = Sd_perm * xp
        e_int[i] = norm(yh - yd) / norm(yd)
        y_manual = zeros(ComplexF64, size(Sd, 1))
        y_manual[ptest] .= yh
        y_sp = Sp * v
        e_wrap[i] = norm(y_sp - y_manual) / norm(y_manual)
        y_dense_perm = zeros(ComplexF64, size(Sd, 1))
        y_dense_perm[ptest] .= Sd_perm * v[ptrial]
        y_dense_true = Sd * v
        e_dense[i] = norm(y_dense_perm - y_dense_true) / norm(y_dense_true)
        e_apply[i] = norm(Sp * v - Sd * v) / norm(Sd * v)
    end
    @printf("mean(Sp.M vs Sd[ptest,ptrial]) = %.12g   max = %.12g\n", mean(e_int), maximum(e_int))
    @printf("mean(Sp vs manual wrapper)     = %.12g   max = %.12g\n", mean(e_wrap), maximum(e_wrap))
    @printf("mean(dense perm wrapper)       = %.12g   max = %.12g\n", mean(e_dense), maximum(e_dense))
    @printf("mean(Sp vs Sd)                 = %.12g   max = %.12g\n", mean(e_apply), maximum(e_apply))
    return e_int, e_wrap, e_dense, e_apply
end

function sort_by_distance(λ, kref)
    p = sortperm(λ, by = x -> abs(x - kref))
    return λ[p], abs.(λ[p] .- kref)
end

ctx = make_ctx_circle(; radius=1.0, Nel=512, minvalues=100)
k0  = 2.404825557695773 + 0im

e_int, e_wrap, e_dense, e_apply = test_perm_wrapper_efie(ctx, k0;
    tol = 1e-4,
    η = 0.25,
    scheduler = SerialScheduler(),
    nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
    farquadstrat  = BEAST.DoubleNumQStrat(5,4),
    q_dense = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
    nvec = 7,
    seed = 1,
)

println("\n================ ACA perm=true RUN (circle) ================\n")
λap1, Vap1, okap1 = test_circle_cavity_Dirichlet_Beyn_ACA_perm_circle(
    radius = 1.0,
    Nel = 512,
    minvalues = 100,
    tol = 1e-4,
    η = 0.25,
    scheduler = SerialScheduler(),
    nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
    farquadstrat  = BEAST.DoubleNumQStrat(5,4),
    z0 = 2.404825557695773 + 0im,
    R = 0.05,
    L = 8,
    N = 32,
    Kmax = 2,
    tol_svd = 1e-8,
    tol_res = 1e-5,
    gmres_rtol = 1e-10,
    gmres_itmax = 1500,
    seed = 1,
    verbose = true,
)

println("\n================ ACA perm=true RUN (ellipse) ================\n")
λap2, Vap2, okap2 = test_circle_cavity_Dirichlet_Beyn_ACA_perm_ellipse(
    radius = 1.0,
    Nel = 512,
    minvalues = 100,
    tol = 1e-4,
    η = 0.25,
    scheduler = SerialScheduler(),
    nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
    farquadstrat  = BEAST.DoubleNumQStrat(5,4),
    z0 = 2.4067799554075626 + 0im,
    Rx = 0.05,
    Ry = 0.02,
    L = 8,
    N = 32,
    Kmax = 2,
    tol_svd = 1e-8,
    tol_res = 1e-5,
    gmres_rtol = 1e-10,
    gmres_itmax = 1500,
    seed = 1,
    verbose = true,
    do_preflight = true,
)

@show okap1
@show okap2
@show λap1
@show λap2

println("\n================ DENSE RUN 1 ================\n")
λd1, Vd1, okd1 = test_circle_cavity_Dirichlet_Beyn_dense(
    radius = 1.0,
    Nel = 512,
    z0 = 2.404825557695773 + 0im,
    R = 0.4,
    L = 8,
    N = 32,
    Kmax = 2,
    tol_svd = 1e-8,
    tol_res = 1e-8,
    seed = 1,
    verbose = true,
)

println("\n================ DENSE RUN 2 ================\n")
λd2, Vd2, okd2 = test_circle_cavity_Dirichlet_Beyn_dense(
    radius = 1.0,
    Nel = 512,
    z0 = 2.4067799554075626 + 0im,
    R = 0.05,
    L = 8,
    N = 32,
    Kmax = 2,
    tol_svd = 1e-8,
    tol_res = 1e-8,
    seed = 1,
    verbose = true,
)

println("\n================ SUMMARY ================\n")
@show okd1
@show okd2
@show λd1
@show λd2

kref = 2.404825557695773 + 0im
λd1s, errd1 = sort_by_distance(λd1, kref)
λd2s, errd2 = sort_by_distance(λd2, kref)

println("\n--- DENSE RUN 1 sorted by distance to kref ---")
@show λd1s
@show errd1

println("\n--- DENSE RUN 2 sorted by distance to kref ---")
@show λd2s
@show errd2
