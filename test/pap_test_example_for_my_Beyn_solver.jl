using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test
using SpecialFunctions
using Gmsh
using FunctionZeros


##

using BEAST
using CompScienceMeshes
using ParallelKMeans
using H2Trees
using AdaptiveCrossApproximation
using OhMyThreads
using Krylov

# Build BlockTree only once (geometry only)
function build_blocktree(X; minvalues=100)
    testtree  = KMeansTree(X.pos, 2; minvalues=minvalues)
    trialtree = KMeansTree(X.pos, 2; minvalues=minvalues)
    return BlockTree(testtree, trialtree)
end

# Fast assembly of K(k) as an HMatrix (ACA/H2)
function hassemble_with_tree(op, X, tree;
                             nearquadstrat = BEAST.DoubleNumSauterQstrat(3,4,0,4,5,6),
                             farquadstrat  = BEAST.DoubleNumQStrat(3,4),
                             scheduler     = DynamicScheduler())
    return AdaptiveCrossApproximation.HMatrix(
        op, X, X, tree;
        nearquadstrat=nearquadstrat, farquadstrat=farquadstrat, scheduler=scheduler
    )
end

# Make a cache keyed by k (Complex)
mutable struct CircleDirichletCtx
    X
    tree
    I   # assembled Identity matrix (k-independent)
    Kcache::Dict{ComplexF64, Any}  # Any because HMatrix type can be long
end

function CircleDirichletCtx(X; minvalues=100)
    tree = build_blocktree(X; minvalues=minvalues)
    I    = assemble(BEAST.Identity(), X, X)   # k-independent
    return CircleDirichletCtx(X, tree, I, Dict{ComplexF64,Any}())
end

# Get or build K(k)
function getK!(ctx::CircleDirichletCtx, k::ComplexF64)
    get!(ctx.Kcache, k) do
        𝒟 = Helmholtz2D.doublelayer(; wavenumber=k)
        @time hassemble_with_tree(𝒟, ctx.X, ctx.tree)
    end
end

# Matrix-free ApplyT
function ApplyT_circle_D(k::ComplexF64, v::AbstractVecOrMat{ComplexF64}, ctx::CircleDirichletCtx)
    K = getK!(ctx, k)
    return (-0.5) .* (ctx.I * v) .+ (K * v)
end

function Beyn_circle_D(k::ComplexF64,
                       VHat::AbstractMatrix{ComplexF64},
                       MInvVHat::AbstractMatrix{ComplexF64},
                       ctx::CircleDirichletCtx;
                       rtol=1e-8, itmax=200, verbose=false)

    n, L = size(VHat)
    @assert size(MInvVHat) == (n, L)

    # Define A(x) = T(k)*x without forming T
    A_mul = x -> ApplyT_circle_D(k, x, ctx)

    # Solve each RHS (simple, robust). Can be upgraded to block-GMRES later.
    for j in 1:L
        b = view(VHat, :, j)
        # Krylov.gmres can accept a function as operator in many setups,
        # but safest is to wrap as LinearOperator if your version needs it.
        x = Krylov.gmres(A_mul, b; rtol=rtol, itmax=itmax, verbose=verbose)
        MInvVHat[:, j] .= x
    end
    return MInvVHat
end

## test

function assemble_T_circle_Dirichlet(k::Number, X)
    𝒟 = Helmholtz2D.doublelayer(; wavenumber = k)
    D  = assemble(𝒟, X, X)
    I  = assemble(BEAST.Identity(), X, X)
    return Matrix(-0.5 .* I .+ D)
end

#test 1 funkcji
a       = 1.0                      # radius
p_order = 2                        # polynomial order of boundary elements
h       = 2π * a / 64 #5000

k = besselj_zero(0, 1) / a   # ≈ 2.4048255577 dla a=1

circle = CompScienceMeshes.meshcircle(a, h)

X = lagrangecx(circle, order=p_order)

k0 = ComplexF64(besselj_zero(0,1)/a)

T_dense = assemble_T_circle_Dirichlet(k0, X)         # your old dense
ctx     = CircleDirichletCtx(X; minvalues=100)

D = size(T_dense,1)
μ_true = randn(ComplexF64, D)
g1 = T_dense * μ_true
g2 = ApplyT_circle_D(k0, μ_true, ctx)

@show norm(g2 - g1)/norm(g1)   # should be small (ACA tolerance dependent)

K = getK!(ctx, k0)              # HMatrix
T = (-0.5) * ctx.I + K          # obiekt “sumy operatorów”
μ_gmres, stats = Krylov.gmres(T, g1; rtol=1e-8, itmax=400, verbose=1)

@show norm(μ_gmres - μ_true)/norm(μ_true)

# test 2
T = assemble_T_circle_Dirichlet(k0, X)
μ_true = randn(ComplexF64, size(T,1))
g = T * μ_true
μ = T \ g
@show norm(μ-μ_true)/norm(μ_true)

K = getK!(ctx, k0)
T̃ = (-0.5)*ctx.I + K

μ_true = randn(ComplexF64, size(ctx.I,1))
g̃ = T̃ * μ_true

μ_gmres, st = Krylov.gmres(T̃, g̃; rtol=1e-10, itmax=400, verbose=1)
@show st.solved st.status
@show norm(μ_gmres - μ_true)/norm(μ_true)

k0 = ComplexF64(2.404825557695773)
𝒟 = Helmholtz2D.doublelayer(; wavenumber=k0)
D_dense = assemble(𝒟, X, X)
K_h = getK!(ctx, k0)

v = randn(ComplexF64, size(D_dense,2))
@show norm(K_h*v - D_dense*v)/norm(D_dense*v)


using Printf
using Statistics

function relerr_K_vs_dense(ctx, X, k; trials=5)
    kk = ComplexF64(k)
    𝒟 = Helmholtz2D.doublelayer(; wavenumber=kk)
    D_dense = assemble(𝒟, X, X)
    K_h = getK!(ctx, kk)
    errs = Float64[]
    for _ in 1:trials
        v = randn(ComplexF64, size(D_dense,2))
        push!(errs, norm(K_h*v - D_dense*v)/norm(D_dense*v))
    end
    return minimum(errs), median(errs), maximum(errs)
end

for k in (2.2, 2.35, 2.4048255577, 2.45, 2.6)
    mn, md, mx = relerr_K_vs_dense(ctx, X, k; trials=5)
    @printf("k=%.4f  relerr(min/med/max) = %.2e / %.2e / %.2e\n", k, mn, md, mx)
end

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

cd, cf, diff, idx = inspect_one_column_efie(ctx, X; k=khat, col=901)