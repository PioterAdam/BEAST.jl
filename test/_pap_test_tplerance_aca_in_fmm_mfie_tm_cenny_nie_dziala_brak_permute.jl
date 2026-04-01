# pap_test_aca_sweep.jl
#
#nie dziala bo problemy z perm, mimo ustalenie false nizaleznie od scheulera ACA chce permutowac
#


using LinearAlgebra, SparseArrays, Random, Statistics, Printf
using BEAST, CompScienceMeshes
using AdaptiveCrossApproximation
using OhMyThreads   # StaticScheduler(), DynamicScheduler()
using Krylov 

using BEAST
using CompScienceMeshes
using ParallelKMeans
using H2Trees
using AdaptiveCrossApproximation
using OhMyThreads
using Krylov

using Printf
using Statistics
# ============================================================
# Context: geometry + tree + (sparse) identity + cache
# ============================================================

#import AdaptiveCrossApproximation: permute

# 1) Najpierw próbujemy użyć getindex (często space[idx] w BEAST działa jako reordering)
#function permute(space::BEAST.LagrangeBasis, p::AbstractVector{<:Integer})
#    return space[p]   # jeśli BEAST ma getindex dla basis -> działa od razu
#end

#import AdaptiveCrossApproximation: permute

#function permute(space::BEAST.LagrangeBasis, p::AbstractVector{<:Integer})
#    return BEAST.subspace(space, collect(p))
#end

#function permute(space::BEAST.LagrangeBasis, p::AbstractVector{<:Integer})
#    return BEAST.restrict(space, collect(p))
#end


#using BEAST
#import AdaptiveCrossApproximation: permute

#=
"""
Fallback permute for BEAST.LagrangeBasis used by ACA/HMatrix.
Tries: getindex (space[p]) -> BEAST.subspace -> BEAST.restrict
"""
function permute(space::BEAST.LagrangeBasis, p::AbstractVector{<:Integer})
    # 1) Try getindex: space[p]
    try
        return space[p]
    catch
        # ignore and try other options
    end

    # 2) Try BEAST.subspace(space, p)
    if isdefined(BEAST, :subspace) && hasmethod(BEAST.subspace, Tuple{typeof(space), Vector{Int}})
        return BEAST.subspace(space, collect(Int, p))
    end

    # 3) Try BEAST.restrict(space, p)
    if isdefined(BEAST, :restrict) && hasmethod(BEAST.restrict, Tuple{typeof(space), Vector{Int}})
        return BEAST.restrict(space, collect(Int, p))
    end

    error("permute(::BEAST.LagrangeBasis, p) not supported: no space[p], BEAST.subspace, or BEAST.restrict available.")
end
=#
#=
using BEAST
import AdaptiveCrossApproximation: permute

"""
Permute a BEAST.LagrangeBasis by a DOF permutation vector p.
We reorder only those fields that look "dof-indexed" (length == ndofs).
This is a pragmatic adapter for ACA/HMatrix.
"""
function permute(space::BEAST.LagrangeBasis, p::AbstractVector{<:Integer})
    pI = collect(Int, p)
    n  = length(pI)

    # Fast path: identity permutation -> no-op
    if all(pI .== 1:n)
        return space
    end

    T   = typeof(space)
    fns = fieldnames(T)
    vals = Any[getfield(space, i) for i in 1:length(fns)]

    # Heuristic: any Vector with length == ndofs is considered dof-indexed and permuted.
    # This covers typical BEAST internals (dofs list, dof points, etc.).
    for i in eachindex(vals)
        v = vals[i]
        if v isa AbstractVector && length(v) == n
            # Only permute if it is an "indexable" container; skip scalars / matrices etc.
            try
                vals[i] = v[pI]
            catch
                # leave as-is if indexing is not supported
            end
        end
    end

    # Rebuild basis (works if LagrangeBasis is a plain struct with default outer ctor)
    return T(vals...)
end

import AdaptiveCrossApproximation: PermutedHMatrix, HMatrix

# jeśli w paczce jest PermutedHMatrix(p, q, H), a wołane jest PermutedHMatrix((p,q), H)
function PermutedHMatrix(pq::Tuple, H::HMatrix)
    return PermutedHMatrix(pq..., H)
end
=#

#=
using LinearAlgebra

import AdaptiveCrossApproximation: PermutedHMatrix

"""
Concrete wrapper that represents P' * A * Q, where
(Q*x) == x[q] and (P*y) == y[p].
So y = P'*(A*(Q*x)) is implemented via invperm(p).
"""
struct PermutedHMatrixImpl{T,MT<:AbstractMatrix{T}} <: AbstractMatrix{T}
    A::MT
    p::Vector{Int}
    q::Vector{Int}
    pinv::Vector{Int}
end

Base.size(PH::PermutedHMatrixImpl) = size(PH.A)
Base.eltype(PH::PermutedHMatrixImpl{T}) where {T} = T

function LinearAlgebra.mul!(y::AbstractVector, PH::PermutedHMatrixImpl, x::AbstractVector)
    @assert length(x) == size(PH,2)
    @assert length(y) == size(PH,1)

    # x2 = Q*x  (permute trial)
    x2 = x[PH.q]

    # y2 = A*x2  (in permuted test ordering)
    y2 = PH.A * x2

    # y = P'*y2  (back to original test ordering)
    y .= y2[PH.pinv]
    return y
end

Base.:*(PH::PermutedHMatrixImpl, x::AbstractVector) = mul!(similar(x, eltype(PH), size(PH,1)), PH, x)

# ---- The missing constructor ACA tries to call: PermutedHMatrix((p,q), A) ----
function (::Type{PermutedHMatrix})(pq::Tuple{<:AbstractVector{<:Integer},<:AbstractVector{<:Integer}}, A::AbstractMatrix)
    p = collect(Int, pq[1])
    q = collect(Int, pq[2])
    return PermutedHMatrixImpl(A, p, q, invperm(p))
end

# Optional: if ACA sometimes calls PermutedHMatrix(p, q, A) instead:
function (::Type{PermutedHMatrix})(p::AbstractVector{<:Integer}, q::AbstractVector{<:Integer}, A::AbstractMatrix)
    return PermutedHMatrix((p,q), A)
end
=#

using LinearAlgebra
import AdaptiveCrossApproximation: PermutedHMatrix

"""
Wrapper reprezentujący P' * A * Q.
p, q to permutacje (wektory indeksów), pinv = invperm(p).
Działa nawet gdy A nie jest AbstractMatrix (np. ACA.HMatrix),
o ile A*x działa (fallback) albo mul!(y,A,x) jest zdefiniowane.
"""
struct PermutedHMatrixImpl{AType,T} <: AbstractMatrix{T}
    A::AType
    p::Vector{Int}
    q::Vector{Int}
    pinv::Vector{Int}
end

Base.size(PH::PermutedHMatrixImpl) = size(PH.A)
Base.eltype(::PermutedHMatrixImpl{AType,T}) where {AType,T} = T

# --- mul! w stylu BLAS: y := α*(PH*x) + β*y ---
function LinearAlgebra.mul!(y::AbstractVector{T}, PH::PermutedHMatrixImpl{AType,T},
                            x::AbstractVector{T}, α::Bool=true, β::Bool=false) where {AType,T}
    @assert length(x) == size(PH,2)
    @assert length(y) == size(PH,1)

    αv = α ? one(T) : zero(T)
    βv = β ? one(T) : zero(T)

    # x2 = Q*x
    x2 = x[PH.q]

    # y2 = A*x2  (fallback: A*x2)
    y2 = Vector{T}(undef, size(PH,1))
    try
        mul!(y2, PH.A, x2)           # jeśli istnieje mul!(y,A,x)
    catch
        y2 .= PH.A * x2              # fallback jeśli jest tylko *(A,x)
    end

    # ytmp = P' * y2
    ytmp = y2[PH.pinv]

    # y = α*ytmp + β*y
    if βv == zero(T)
        y .= αv .* ytmp
    else
        @inbounds for i in eachindex(y)
            y[i] = αv*ytmp[i] + βv*y[i]
        end
    end
    return y
end

# standardowy operator *
function Base.:*(PH::PermutedHMatrixImpl{AType,T}, x::AbstractVector{T}) where {AType,T}
    y = Vector{T}(undef, size(PH,1))
    LinearAlgebra.mul!(y, PH, x, true, false)
    return y
end

# ============================================================
# KLUCZ: konstruktor PermutedHMatrix łapie też A::HMatrix (nie-AbstractMatrix)
# ============================================================
function (::Type{PermutedHMatrix})(pq::Tuple{<:AbstractVector{<:Integer},<:AbstractVector{<:Integer}}, A)
    p = collect(Int, pq[1])
    q = collect(Int, pq[2])
    T = eltype(A)   # ważne: typ elementów
    return PermutedHMatrixImpl{typeof(A),T}(A, p, q, invperm(p))
end

# na wszelki wypadek: PermutedHMatrix(p,q,A)
function (::Type{PermutedHMatrix})(p::AbstractVector{<:Integer}, q::AbstractVector{<:Integer}, A)
    return PermutedHMatrix((p,q), A)
end



struct CircleDirichletCtx
    X
    tree
    I::SparseMatrixCSC{Float64,Int}
    cache::Dict{ComplexF64,Any}
end

"""
Create circle mesh with Nel segments and Lagrange basis X.
Build a blocktree once.
Use sparse identity I (same "I" used in dense and FMM builds).
"""
function make_ctx_circle(; radius::Float64=1.0, Nel::Int=512)
    # Mesh of a circle boundary (CompScienceMeshes)
    # NOTE: depending on your CMS version, circle signature can be circle(radius, Nel)
    h = 2π*radius / Nel
    Γ = CompScienceMeshes.meshcircle(radius, h)

    # Basis on boundary
    X = BEAST.lagrangec0d1(Γ)

    # Block tree for HMatrix (your environment already has build_blocktree)
    tree = build_blocktree(X) ## tree = build_blocktree(X; perm=false) tree = build_blocktree(X; reorder=false) 
    #tree = build_blocktree(X; perm=false)
    #tree = build_blocktree(X; reorder=false)


    n = numfunctions(X)
    I = spdiagm(0 => ones(Float64, n))  # sparse identity

    return CircleDirichletCtx(X, tree, I, Dict{ComplexF64,Any}())
end

# ============================================================
# Assemblers
# ============================================================

"""
Dense operator T(k) = -0.5*I + K, where K is double-layer (assembled dense).
quadstrat: e.g. DoubleNumSauterQstrat(...)
"""
function build_dense_T(ctx::CircleDirichletCtx, k::ComplexF64; quadstrat)
    Dop = Helmholtz2D.doublelayer(; wavenumber = k)
    K = assemble(Dop, ctx.X, ctx.X; quadstrat=quadstrat)
    return (-0.5) * ctx.I + K
end

"""
HMatrix operator T_FMM/ACA(k) = -0.5*I + Kfmm, where Kfmm is ACA/HMatrix.
nearquadstrat: Sauter strategy for near-singular
farquadstrat : regular quadrature (non-singular)
scheduler     : StaticScheduler() or DynamicScheduler()
tol           : ACA tolerance (Frobenius-like target)
perm          : MUST be false unless ACA has permute implemented for this space
"""
#function build_fmm_T(ctx::CircleDirichletCtx, k::ComplexF64;
#                     tol::Float64=1e-4,
#                     scheduler=StaticScheduler(),
#                     nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
#                     farquadstrat  = BEAST.DoubleNumQStrat(4,4),
#                     perm::Bool=false)
#
#    Dop = Helmholtz2D.doublelayer(; wavenumber = k)
#
#    # HMatrix assembly
#    Kfmm = AdaptiveCrossApproximation.HMatrix(Dop, ctx.X, ctx.X, ctx.tree;
#        compressor    = ACA(; tol=tol),
#        nearquadstrat = nearquadstrat,
#        farquadstrat  = farquadstrat,
#        scheduler     = scheduler,
#        perm          = false,      # <<< KLUCZ,
#    )
#
#    return (-0.5) * ctx.I + Kfmm
#end

function build_fmm_T(ctx::CircleDirichletCtx, k::ComplexF64;
                     tol::Float64=1e-4,
                     #scheduler=StaticScheduler(),
                     scheduler=DynamicScheduler(),
                     nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                     farquadstrat  = BEAST.DoubleNumQStrat(4,4))

    Dop = Helmholtz2D.doublelayer(; wavenumber = k)

    Kfmm = AdaptiveCrossApproximation.HMatrix(Dop, ctx.X, ctx.X, ctx.tree;
        compressor    = ACA(; tol=tol),
        nearquadstrat = nearquadstrat,
        farquadstrat  = farquadstrat,
        scheduler     = scheduler,
        perm          = false,   # twardo
    )

    return (-0.5) * ctx.I + Kfmm
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

# ============================================================
# Single-shot "source of error" + sweep over ACA tol/schedulers
# ============================================================

function singleshot_compare(ctx::CircleDirichletCtx, k::ComplexF64;
                            q_lo = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                            q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
                            tol = 1e-4,
                            scheduler = StaticScheduler(),
                            nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                            farquadstrat  = BEAST.DoubleNumQStrat(4,4),
                            nvec::Int=7, seed::Int=1)

    @printf("\n=== singleshot_compare ===\n")
    @printf("k = %+.12f%+.12fi\n", real(k), imag(k))
    @printf("scheduler = %s, ACA tol = %.1e\n", string(typeof(scheduler)), tol)

    Td_lo = build_dense_T(ctx, k; quadstrat=q_lo)
    Td_hi = build_dense_T(ctx, k; quadstrat=q_hi)

    TF = build_fmm_T(ctx, k;
        tol=tol,
        scheduler=scheduler,
        nearquadstrat=nearquadstrat,
        farquadstrat=farquadstrat,
        #perm=false,
    )

    e_fmm   = relerrs_apply(TF,   Td_hi; nvec=nvec, seed=seed)
    e_dense = relerrs_apply(Td_lo, Td_hi; nvec=nvec, seed=seed)

    @printf("mean(fmm vs dense_hi)  = %.12g   max = %.12g\n", mean(e_fmm), maximum(e_fmm))
    @printf("mean(dense_lo vs hi)   = %.12g   max = %.12g\n", mean(e_dense), maximum(e_dense))

    return e_fmm, e_dense
end

function sweep_aca(ctx::CircleDirichletCtx, k::ComplexF64;
                   tols = [1e-3, 1e-4, 1e-5, 1e-6],
                   schedulers = [StaticScheduler(), DynamicScheduler()], #  
                   q_lo = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                   q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
                   nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
                   farquadstrat  = BEAST.DoubleNumQStrat(4,4),
                   nvec::Int=7, seed::Int=1)

    @printf("\n============================================\n")
    @printf("ACA sweep (single k, fixed X/tree)\n")
    @printf("k = %+.12f%+.12fi\n", real(k), imag(k))
    @printf("============================================\n")

    # Dense reference once
    Td_lo = build_dense_T(ctx, k; quadstrat=q_lo)
    Td_hi = build_dense_T(ctx, k; quadstrat=q_hi)
    e_dense = relerrs_apply(Td_lo, Td_hi; nvec=nvec, seed=seed)
    @printf("Dense self-check: mean=%.3e  max=%.3e\n", mean(e_dense), maximum(e_dense))

    for sch in schedulers
        @printf("\n=== scheduler = %s ===\n", string(typeof(sch)))
        for tol in tols
            try
                TF = build_fmm_T(ctx, k;
                    tol=tol,
                    scheduler=sch,
                    nearquadstrat=nearquadstrat,
                    farquadstrat=farquadstrat,
                    #perm=false,
                )
                e_fmm = relerrs_apply(TF, Td_hi; nvec=nvec, seed=seed)
                @printf("tol=%-8.1e  mean(FMM vs dense_hi)=%.3e  max=%.3e\n",
                        tol, mean(e_fmm), maximum(e_fmm))
            catch err
                @printf("tol=%-8.1e  ERROR: %s\n", tol, sprint(showerror, err))
            end
        end
    end

    return nothing
end

# Build BlockTree only once (geometry only)
function build_blocktree(X; minvalues=100)
    testtree  = KMeansTree(X.pos, 2; minvalues=minvalues)
    trialtree = KMeansTree(X.pos, 2; minvalues=minvalues)
    return BlockTree(testtree, trialtree)
end

# ============================================================
# MAIN
# ============================================================

function main()
    # choose discretization here:
    Nel = 512  # try 246 / 512 / 1024
    ctx = make_ctx_circle(; radius=1.0, Nel=Nel)

    k = ComplexF64(2.404825557695773, 0.0)

    # quick single-shot (your “source-of-error decider”)
    singleshot_compare(ctx, k;
        q_lo = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        q_hi = BEAST.DoubleNumSauterQstrat(10,10,0,4,15,15),
        tol = 1e-4,
        scheduler = StaticScheduler(),
        nearquadstrat = BEAST.DoubleNumSauterQstrat(6,7,0,4,8,9),
        farquadstrat  = BEAST.DoubleNumQStrat(4,4),
        nvec=7, seed=123
    )

    # full sweep ACA tol + schedulers
    sweep_aca(ctx, k;
        tols=[1e-3, 1e-4, 1e-5, 1e-6],
        schedulers=[StaticScheduler(), DynamicScheduler()],
        nvec=7, seed=123
    )
end

main()