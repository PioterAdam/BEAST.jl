# minimal_singleshot.jl
using LinearAlgebra, Random
using BEAST, CompScienceMeshes

#FMM
using BEAST
using CompScienceMeshes
using ParallelKMeans
using H2Trees
using AdaptiveCrossApproximation
using OhMyThreads
using Krylov

using Printf
using Statistics

"""
Minimal single-shot:
- buduje Tfmm (z ctx + getK!)
- buduje Tdense_lo i Tdense_hi
- liczy rel. błędy na losowych v:
    ||(Tfmm - Td_hi)v|| / ||Td_hi v||
    ||(Td_lo - Td_hi)v|| / ||Td_hi v||
"""

# Build BlockTree only once (geometry only)
function build_blocktree(X; minvalues=100)
    testtree  = KMeansTree(X.pos, 2; minvalues=minvalues)
    trialtree = KMeansTree(X.pos, 2; minvalues=minvalues)
    return BlockTree(testtree, trialtree)
end

# Fast assembly of K(k) as an HMatrix (ACA/H2)
#"Only quadratures of degree 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25 and 30 available."
function hassemble_with_tree(op, X, tree;
                             #nearquadstrat = BEAST.DoubleNumSauterQstrat(3,4,0,4,5,6),
                             #nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                             nearquadstrat = BEAST.DoubleNumSauterQstrat(7,7,0,4,9,9),
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

## test dla kwadratury

function singleshot(ctx, X, khat::ComplexF64;
        q_lo = BEAST.DoubleNumSauterQstrat(3,4,0,4,5,6),
        q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
        nvec::Int = 5, seed::Int = 1)

    # FMM operator
    Tfmm = (-0.5)*ctx.I + getK!(ctx, khat)

    # Dense operators
    Dop = Helmholtz2D.doublelayer(; wavenumber = khat)
    Klo = assemble(Dop, X, X; quadstrat=q_lo)
    Khi = assemble(Dop, X, X; quadstrat=q_hi)
    Td_lo = (-0.5)*ctx.I + Klo
    Td_hi = (-0.5)*ctx.I + Khi

    D = size(Td_hi, 2)
    Random.seed!(seed)

    ef = zeros(Float64, nvec)
    ed = zeros(Float64, nvec)

    for j in 1:nvec
        v = randn(ComplexF64, D)
        hiv = Td_hi * v
        ef[j] = norm((Tfmm  - Td_hi) * v) / (norm(hiv) + eps())
        ed[j] = norm((Td_lo - Td_hi) * v) / (norm(hiv) + eps())
    end

    println("mean(fmm vs dense_hi)  = ", mean(ef), "   max = ", maximum(ef))
    println("mean(dense_lo vs hi)   = ", mean(ed), "   max = ", maximum(ed))

    return ef, ed
end




a       = 1.0                      # radius
p_order = 2                        # polynomial order of boundary elements
h       = 2π * a / (2*246)            # mesh size along boundary

#circle  = circle_curvilinear(a, p_order; h = h)
circle = CompScienceMeshes.meshcircle(a, h)

# Use your C0 scalar basis on the boundary
# (If you prefer the "d1 linear" variant, swap the next line)
#X = lagrangec0_curvilinear(circle)
#X = lagrangec0d1_curvilinear_exact_linear(circle)
X = lagrangecx(circle, order=p_order)
#@assert hasproperty(X, :pos)
#@assert size(X.pos, 1) == 2   # dla 2D

# Build context (tree + I + cache)
ctx = CircleDirichletCtx(X; minvalues=100)

khat = ComplexF64(2.404825557695773 + 0im)  # albo Twoje khat
ef, ed = singleshot(ctx, X, khat; nvec=7, seed=1)


#test dla ACA
using LinearAlgebra, Random, Statistics, Printf
using BEAST, CompScienceMeshes
using AdaptiveCrossApproximation
using OhMyThreads   # StaticScheduler / DynamicScheduler

using AdaptiveCrossApproximation
using OhMyThreads

function hassemble_with_tree(op, X, tree;
    compressor = ACA(; tol=1e-4),                 # <-- DODANE
    nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
    farquadstrat  = BEAST.DoubleNumQStrat(4,4),
    scheduler     = StaticScheduler(),
    perm::Bool    = false,
)
    return AdaptiveCrossApproximation.HMatrix(op, X, X, tree;
        compressor    = compressor,               # <-- PRZEKAZANIE
        nearquadstrat = nearquadstrat,
        farquadstrat  = farquadstrat,
        scheduler     = scheduler,
        perm          = false,     # <-- dodaj na sztywno,
    )
end

# ------------------------------------------------------------
# Helper: build dense operator (control)
# ------------------------------------------------------------
function build_dense_T(X, k::ComplexF64; quadstrat)
    Dop = Helmholtz2D.doublelayer(; wavenumber = k)
    Kdense = assemble(Dop, X, X; quadstrat = quadstrat)
    Tdense = (-0.5) * I(size(Kdense,1)) + Kdense
    return Tdense
end

# ------------------------------------------------------------
# Helper: apply FMM/ACA operator (you MUST adapt this to your code)
#   Option A: you already have hassemble_with_tree(op, X, tree; ...)
#   Option B: you create HMatrix(...) directly
# ------------------------------------------------------------
"""
Return Tdense-like object TFMM = -0.5*I + Kfmm (Kfmm from ACA/HMatrix).
You MUST edit the body to call your own hassemble_with_tree / HMatrix constructor.
"""
function build_fmm_T(ctx, X, k::ComplexF64;
                     tol::Float64,
                     scheduler,
                     nearquadstrat,
                     farquadstrat)

    Dop = Helmholtz2D.doublelayer(; wavenumber = k)

    # --- EDIT THIS PART to match your codebase ---
    # If you have: Kfmm = hassemble_with_tree(Dop, X, ctx.tree; compressor=..., scheduler=...)
    Kfmm = hassemble_with_tree(Dop, X, ctx.tree;
        compressor = ACA(; tol = tol),
        nearquadstrat = nearquadstrat,
        farquadstrat  = farquadstrat,
        scheduler     = scheduler
    )

    # If Kfmm is an HMatrix-like object supporting *(vector), this works:
    TFMM = (-0.5) * ctx.I + Kfmm
    return TFMM
end

# ------------------------------------------------------------
# Main single-shot metric:
#   relerr(TA, TB) on random vectors:
#     ||(TA - TB)v|| / ||TB v||
# ------------------------------------------------------------
function relerrs(TA, TB; nvec=7, seed=1)
    Random.seed!(seed)
    n = size(TB, 1)
    errs = Vector{Float64}(undef, nvec)
    for i in 1:nvec
        v = randn(ComplexF64, n)
        num = norm((TA - TB) * v)
        den = norm(TB * v)
        errs[i] = num / den
    end
    return errs
end

# ------------------------------------------------------------
# Sweep ACA tolerance for fixed k, X, tree, quadstrats
# ------------------------------------------------------------
function sweep_aca(ctx, X, k::ComplexF64;
                   tols = [1e-3, 1e-4, 1e-5, 1e-6],
                   schedulers = [StaticScheduler(), DynamicScheduler()],
                   q_lo = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                   q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
                   nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                   farquadstrat  = BEAST.DoubleNumQStrat(4,4),
                   nvec = 7,
                   seed = 1)

    @printf("k = %+.12f%+.12fi\n", real(k), imag(k))

    # Dense controls (once)
    Td_lo = build_dense_T(X, k; quadstrat=q_lo)
    Td_hi = build_dense_T(X, k; quadstrat=q_hi)

    e_dense = relerrs(Td_lo, Td_hi; nvec=nvec, seed=seed)
    @printf("Dense self-check: mean=%.3e  max=%.3e\n", mean(e_dense), maximum(e_dense))

    # Sweep
    for sch in schedulers
        @printf("\n=== scheduler = %s ===\n", typeof(sch))
        for tol in tols
            TF = build_fmm_T(ctx, X, k;
                tol=tol, scheduler=sch,
                nearquadstrat=nearquadstrat,
                farquadstrat=farquadstrat
            )

            e_fmm = relerrs(TF, Td_hi; nvec=nvec, seed=seed)

            @printf("tol=%-8.1e  mean(FMM vs dense_hi)=%.3e  max=%.3e\n",
                    tol, mean(e_fmm), maximum(e_fmm))
        end
    end

    return nothing
end

# ------------------------------------------------------------
# Example call (YOU set ctx, X, k)
# ------------------------------------------------------------
k = ComplexF64(2.404825557695773, 0.0)
 sweep_aca(ctx, X, k;
     tols=[1e-3,1e-4,1e-5,1e-6],
     schedulers=[StaticScheduler(), DynamicScheduler()],
     nvec=7, seed=123
)



using LinearAlgebra, Statistics, Random, Printf

function op_mismatch_random(ctx, X; k::ComplexF64, ntests::Int=5, seed::Int=1,
                            nearquad_dense = BEAST.DoubleNumSauterQstrat(8,9,0,4,9,10))
    rng = MersenneTwister(seed)
    D = size(ctx.I, 1)

    # --- FMM/ACA operator (action) ---
    Kfmm = getK!(ctx, k)
    Tfmm_mul = v -> begin
        y = similar(v)
        mul!(y, Kfmm, v)                 # y = Kfmm*v
        mul!(y, ctx.I, v, -0.5+0im, 1+0im) # y += (-0.5)*I*v
        y
    end

    # --- Dense operator (assemble once) ---
    Dop = Helmholtz2D.doublelayer(; wavenumber=k)
    Kdense = assemble(Dop, X, X; quadstrat=nearquad_dense)
    Tdense = (-0.5) * ctx.I + Kdense

    errs = zeros(Float64, ntests)
    for t in 1:ntests
        v = randn(rng, ComplexF64, D)
        Tv_dense = Tdense * v
        Tv_fmm   = Tfmm_mul(v)
        errs[t]  = norm(Tv_fmm - Tv_dense) / max(norm(Tv_dense), 1e-300)
        @printf("test %d: err = %.3e\n", t, errs[t])
    end

    @printf("err mean = %.3e, max = %.3e  (k = %.10f%+.10fi)\n",
            mean(errs), maximum(errs), real(k), imag(k))
    return errs
end

# ---- call ----
k = ComplexF64(khat + 0im)   # albo khat jeśli chcesz
errs = op_mismatch_random(ctx, X; k=k, ntests=5, seed=1)