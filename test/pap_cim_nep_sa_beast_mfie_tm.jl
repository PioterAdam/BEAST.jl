
using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test
using SpecialFunctions
using Gmsh


#import CompScienceMeshes: measure, nodes, Simplex
import LinearAlgebra: norm
#import CompScienceMeshes: indices, circle, chart

#import BEAST: HelmholtzOperator2D, kernelvals, KernelValsHelmholtz2D
#import SpecialFunctions: hankelh2
#import CompScienceMeshes: normal
#using LinearAlgebra: norm, dot

const Shape = BEAST.Shape
const LagrangeBasis = BEAST.LagrangeBasis

#using LinearAlgebra
using Printf
using Base.Threads

using FunctionZeros
j01 = besselj_zero(0, 1)

#=
# in your version of BEAST probably not required
# measure for straight 1D edges embedded in R^U
function measure(s::Simplex{U,1,<:Any,<:Any,T}) where {U,T}
    p0, p1 = nodes(s)
    return norm(p1 - p0)
end

# in your version of BEAST probably not required
# helper functions 
function circle_curvilinear(radius::Real, porder::Integer; h::Real = 2π*radius/64)
    @assert porder ≥ 1 "porder must be ≥ 1"
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Verbosity", 0)

        gmsh.model.add("circle_p$(porder)")
        s = gmsh.model.occ.addDisk(0.0, 0.0, 0.0, radius, radius)
        gmsh.model.occ.synchronize()

        # sizing + high-order
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", h)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", h)
        gmsh.option.setNumber("Mesh.ElementOrder", porder)
        gmsh.option.setNumber("Mesh.HighOrderOptimize", porder >= 2 ? 2 : 0)
        gmsh.option.setNumber("Mesh.SecondOrderLinear", porder >= 2 ? 0 : 1)

        # physicals
        gmsh.model.addPhysicalGroup(2, [s], 1)
        gmsh.model.setPhysicalName(2, 1, "Domain")
        curves = [t[2] for t in gmsh.model.getBoundary([(2, s)], false, false, false) if t[1] == 1]
        gmsh.model.addPhysicalGroup(1, curves, 2)
        gmsh.model.setPhysicalName(1, 2, "Boundary")

        gmsh.model.mesh.generate(2)

        # nodes on boundary physical
        nb = gmsh.model.mesh.getNodesForPhysicalGroup(1, 2)
        nodeTags   = nb[1]
        nodeCoords = nb[2]
        @assert length(nodeCoords) == 3*length(nodeTags) "unexpected coords size"

        tag2idx = Dict{Int,Int}(Int(nodeTags[i]) => i for i in eachindex(nodeTags))

        verts = Vector{SVector{2,Float64}}(undef, length(nodeTags))
        @inbounds for i in eachindex(nodeTags)
            x = nodeCoords[3(i-1)+1]; y = nodeCoords[3(i-1)+2]
            verts[i] = SVector{2,Float64}(x, y)
        end

        NF = porder + 1
        faces = SVector{NF,Int}[]
        for c in curves
            types, elemTags, elemNodeTags = gmsh.model.mesh.getElements(1, c)
            @inbounds for k in eachindex(types)
                tags  = elemTags[k]
                nodes = elemNodeTags[k]
                ne = length(tags)
                ne == 0 && continue

                nNode_block = div(length(nodes), ne)
                nNode_block == NF || continue

                @assert length(nodes) == ne * nNode_block
                for e in 1:ne
                    off = (e-1) * nNode_block
                    ids = ntuple(j -> tag2idx[Int(nodes[off + j])], NF)
                    push!(faces, SVector{NF,Int}(ids))
                end
            end
        end
        @assert !isempty(faces) "No boundary line elements with p=$(porder) found."

        return circle(verts, faces, porder)
    finally
        gmsh.finalize()
    end
end

# in your version of BEAST probably not required
function lagrangec0d1_curvilinear_exact_linear(mesh::circle)
    T = coordtype(mesh)
    S = BEAST.Shape{T}

    # Corner-only vertex DOFs (ignore mid nodes on curved edges)
    dofs = Dict{Int, Vector{S}}()
    for (c, cell) in pairs(mesh.faces)
        @assert length(cell) ≥ 2
        g1, g2 = cell[1], cell[2]   # corners (s,e)
        push!(get!(dofs, g1, Vector{S}()), S(c, 1, T(1)))
        push!(get!(dofs, g2, Vector{S}()), S(c, 2, T(1)))
    end

    gids = sort!(collect(keys(dofs)))
    fns  = [dofs[g] for g in gids]
    pos  = [mesh.vertices[g] for g in gids]

    return BEAST.LagrangeBasis{1, 0, 2}(mesh, fns, pos)
end
=#

function BeynSolve_auto(z0, Rx, Ry, UserBeynFunction, ApplyT, D;
                        L::Int       = D,
                        N::Int       = 64,
                        Kmax::Int    = 3,
                        tol_svd::Float64 = 1e-8,
                        tol_res::Float64 = 1e-8,
                        verbose::Bool = true)

    # ----------------------------------------------------------
    # 0) Fix probing matrix VHat ONCE for all K
    # ----------------------------------------------------------
    VHat = randn(D, L) .+ im * randn(D, L)

    #=
    ## one thread
    # Precompute contour moments up to 2*Kmax-1:
    #   Ap[p] ≈ (1/(2πi)) ∮ (z-z0)^(p-1) T(z)^{-1} VHat dz
    # p = 1..2*Kmax
    # ----------------------------------------------------------
    nMom_max = 2*Kmax
    Ap       = [zeros(ComplexF64, D, L) for _ in 1:nMom_max]
    MInvVHat = zeros(ComplexF64, D, L)

    Δθ = 2π / N
    for n = 0:N-1
        θ  = n * Δθ
        cθ = cos(θ); sθ = sin(θ)

        w  = Rx * cθ + im * Ry * sθ      # z - z0
        z  = z0 + w
        dz = (-Rx * sθ + im * Ry * cθ) * Δθ

        # T(z)^{-1} VHat:
        MInvVHat = UserBeynFunction(z, VHat, MInvVHat)

        # accumulate powers of w up to nMom_max
        wp = one(w)
        for p = 1:nMom_max
            @inbounds Ap[p] .+= wp * dz .* MInvVHat
            wp *= w
        end
    end
    =#

    ## multi-threads version: using Base.Threads
    # Precompute contour moments up to 2*Kmax-1:
    #   Ap[p] ≈ (1/(2πi)) ∮ (z-z0)^(p-1) T(z)^{-1} VHat dz
    # p = 1..2*Kmax
    # ----------------------------------------------------------
    # Parallelization note:
    # We now parallelize the most expensive loop: the N solves T(z_j)^{-1} * VHat
    # along the contour. This works both for dense matrices and for FMM-based
    # matrix-free operators, provided that:
    #
    #   * `UserBeynFunction` is thread-safe (here each thread gets its own
    #     `MInvVHat` workspace, so there is no race on that array).
    #
    # If, at some point, the FMM / linear solver becomes internally multithreaded,
    # then you should either:
    #   - run Julia with JULIA_NUM_THREADS=1 and let the FMM use its own threads,
    #     or
    #   - disable threading inside the FMM and use threading at the Beyn level
    #     (as done here), to avoid oversubscription and performance loss.
    # ----------------------------------------------------------
    nMom_max = 2*Kmax

    # Global moments (filled after reduction)
    Ap = [zeros(ComplexF64, D, L) for _ in 1:nMom_max]

    # Thread-local accumulators to avoid write conflicts
    nT           = nthreads()
    Ap_thread    = [ [zeros(ComplexF64, D, L) for _ in 1:nMom_max] for _ in 1:nT ]
    MInv_thread  = [ zeros(ComplexF64, D, L) for _ in 1:nT ]

    Δθ = 2π / N

    @threads for nn in 0:(N-1)
        tid = threadid()
        Ap_loc   = Ap_thread[tid]
        MInvVHat = MInv_thread[tid]

        θ  = nn * Δθ
        cθ = cos(θ); sθ = sin(θ)

        w  = Rx * cθ + im * Ry * sθ      # z - z0
        z  = z0 + w
        dz = (-Rx * sθ + im * Ry * cθ) * Δθ

        # T(z)^{-1} VHat:
        MInvVHat = UserBeynFunction(z, VHat, MInvVHat)

        # accumulate powers of w up to nMom_max into thread-local Ap_loc
        wp = one(w)
        @inbounds for p = 1:nMom_max
            Ap_loc[p] .+= wp * dz .* MInvVHat
            wp *= w
        end
    end

    # Reduce thread-local moments into the global Ap
    for t in 1:nT
        Ap_loc = Ap_thread[t]
        @inbounds for p in 1:nMom_max
            Ap[p] .+= Ap_loc[p]
        end
    end


    # ----------------------------------------------------------
    # 1) Internal solver for a fixed K using the precomputed Ap
    # ----------------------------------------------------------
    function Beyn_from_moments(K::Int)
        nMom = 2*K
        m    = D
        ℓ    = L
        Km   = K*m
        Kl   = K*ℓ

        # build block Hankel B0, B1 of size (K*m)×(K*ℓ)
        B0 = zeros(ComplexF64, Km, Kl)
        B1 = zeros(ComplexF64, Km, Kl)

        for i in 1:K
            rows = (i-1)*m+1 : i*m
            for j in 1:K
                cols = (j-1)*ℓ+1 : j*ℓ
                p0 = i + j - 2          # 0..2K-2
                p1 = i + j - 1          # 1..2K-1
                @inbounds begin
                    B0[rows, cols] .= Ap[p0+1]   # Ap is 1-based
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

        # first block of size m×k contains the physical space
        Vblock1 = U0[1:m, :]
        Vfull   = Vblock1 * S
        λ       = μ .+ z0

        # residuals and good indices
        r  = zeros(Float64, length(λ))
        good_idxs = Int[]

        if ApplyT === nothing
            # can't check
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

    # ----------------------------------------------------------
    # 2) Adaptive loop over K (based on #good eigenpairs)
    # ----------------------------------------------------------
    λ_best    = ComplexF64[]
    V_best    = zeros(ComplexF64, D, 0)
    good_best = 0
    good_prev = 0
    K_used    = 0

    for K in 1:Kmax
        λ, V, good_idxs, r = Beyn_from_moments(K)
        ngood = length(good_idxs)

        if verbose
            @printf("  K=%d: %d good eigenpairs (res ≤ %.1e)\n",
                    K, ngood, tol_res)
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
                @printf("  #good eigenpairs stabilized at %d for K=%d.\n",
                        ngood, K)
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

    # sort the good eigenpairs
    p = sortperm(λ_best, by = x -> (real(x), imag(x)))
    λ_best = λ_best[p]
    V_best = V_best[:, p]

    if verbose
        @printf("Using K=%d with %d good eigenpairs.\n", K_used, good_best)
    end

    return λ_best, V_best, true
end


# circular wrapper
function BeynSolve_auto(z0, R, UserBeynFunction, ApplyT, D; kwargs...)
    return BeynSolve_auto(z0, R, R, UserBeynFunction, ApplyT, D; kwargs...)
end

"""
    AsakuraSolve_circle_auto(z0, R, UserFun, ApplyT, D;
                             L        = 8,
                             Kblocks  = 4,
                             N        = 64,
                             δ_svd    = 1e-12,
                             tol_res  = 1e-8,
                             verbose  = true)

Contour-integral eigenvalue solver in the Asakura / block Sakurai–Sugiura style
for NEPs F(λ)x = 0 on a circular contour

    Γ = { z : z = z0 + R * exp(iθ), 0 ≤ θ < 2π }.

Arguments
---------
- `z0::Complex`: center of the contour (γ in the paper).
- `R::Float64`: radius of the contour (ρ in the paper).
- `UserFun(z, VHat, MInvVHat)`: must return `F(z) \\ VHat` (size D×L).
- `ApplyT(z, v)`: must return `F(z) * v` for residual checks.
- `D::Int`: dimension of the NEP (size of F(z)).

Important keyword parameters
----------------------------
- `L`: block size (number of right-hand sides V).
- `Kblocks`: number of block moments (called K or m̃ in the paper);
            total reduced dimension is `Kblocks * L`.
- `N`: number of quadrature points on the contour (trapezoidal rule).
- `δ_svd`: SVD truncation threshold: keep σ_j ≥ δ_svd * max(σ).
- `tol_res`: residual tolerance ‖F(λ̂_j)x̂_j‖ / ‖x̂_j‖ ≤ tol_res.
- `verbose`: print basic diagnostic info.

Returns
-------
`λ_good::Vector{ComplexF64}, X_good::Matrix{ComplexF64}, ok::Bool`

- `λ_good`: eigenvalues inside the contour that passed the residual test.
- `X_good`: corresponding eigenvectors (columns).
- `ok`: `true` if at least one eigenpair passed the residual test.
"""
function AsakuraSolve_circle_auto(z0::Complex, R::Float64,
                                  UserFun, ApplyT, D::Int;
                                  L::Int         = 8,
                                  Kblocks::Int   = 4,
                                  N::Int         = 64,
                                  δ_svd::Float64 = 1e-12,
                                  tol_res::Float64 = 1e-8,
                                  verbose::Bool  = true)

    # --- basic sanity checks -------------------------------------------------
    Kblocks ≥ 1 || error("Kblocks must be ≥ 1")
    L       ≥ 1 || error("L must be ≥ 1")
    N       ≥ 4 || error("N must be ≥ 4")

    KL = Kblocks * L
    nMom = 2*Kblocks              # number of moments S_0,...,S_{2K-1}

    verbose && println("=== Asakura / block-SS contour solver (circle) ===")
    verbose && @printf("Center z0 = %+.6e%+.6ei, R = %.3f, D = %d, L = %d, Kblocks = %d, N = %d\n",
                       real(z0), imag(z0), R, D, L, Kblocks, N)

    # --- Step 1: random probing block V -------------------------------------
    V   = randn(ComplexF64, D, L)
    S_blocks = [zeros(ComplexF64, D, L) for _ in 1:nMom]
    MInv = zeros(ComplexF64, D, L)

    # --- Step 2: contour integration (same parameterisation as Beyn) --------
    #
    # Circle:
    #   z(θ) = z0 + R*exp(iθ), θ ∈ [0,2π)
    # dz = i R exp(iθ) dθ
    #
    # S_k = (1/(2π i)) ∮ (z-z0)^k F(z)^{-1} V dz
    #     ≈ Σ_j (w_j^k) F(z_j)^{-1} V * (dz_j/(2π i))
    #
    Δθ = 2π / N

    for n = 0:N-1
        θ  = (n + 0.5) * Δθ
        ζ  = exp(im * θ)
        z  = z0 + R * ζ
        dz = im * R * ζ * Δθ        # dz = z'(θ) dθ

        # weight for (1/(2π i)) dz
        w_dθ = dz / (2π * im)

        # w = z - z0
        w = z - z0

        # F(z)^{-1} V
        G = UserFun(z, V, MInv)     # D×L

        # accumulate moments
        wpow = one(w)
        for k = 0:(nMom-1)
            S_blocks[k+1] .+= (wpow * w_dθ) .* G
            wpow *= w
        end
    end

    # --- Step 3: build block Hankel H, H< from S_k --------------------------
    #
    # H, H< are (Kblocks*D) × (Kblocks*L):
    #   H_ij   = S_{i+j-2}
    #   H<_ij  = S_{i+j-1}
    #
    rows_H = Kblocks * D
    cols_H = Kblocks * L

    H   = zeros(ComplexF64, rows_H, cols_H)
    Hlt = zeros(ComplexF64, rows_H, cols_H)

    for i in 1:Kblocks
        for j in 1:Kblocks
            Sk   = S_blocks[i + j - 1]   # S_{i+j-2}
            Sk1  = S_blocks[i + j]       # S_{i+j-1}
            rows = (i-1)*D+1 : i*D
            cols = (j-1)*L+1 : j*L
            @inbounds begin
                H[rows, cols]   .= Sk
                Hlt[rows, cols] .= Sk1
            end
        end
    end

    # --- Step 4: SVD(H) and truncation --------------------------------------
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

    UH1 = UH[:, keep]             # (Kblocks*D) × m'
    VH1 = VH[:, keep]             # (Kblocks*L) × m'
    Σ1  = Diagonal(σ[keep])       # m'×m'

    # --- Step 5: reduced small pencil ---------------------------------------
    # G = Σ1^{-1} * (UH1' * Hlt * VH1)
    Hlt_tilde = UH1' * Hlt * VH1  # m'×m'
    G        = Σ1 \ Hlt_tilde     # m'×m'

    eigG = eigen(G)
    μ    = eigG.values            # eigenvalues in shifted coordinate (w-plane)
    Y    = eigG.vectors           # m'×m'

    # --- Step 6: map μ to λ and reconstruct eigenvectors --------------------
    # λ = z0 + μ
    λ_hat = z0 .+ μ

    # Build S = [S_0 S_1 ... S_{Kblocks-1}] ∈ ℂ^{D × (Kblocks*L)}
    ##S_matrix = hcat(S_blocks[1:Kblocks]...)   # D×(Kblocks*L)
    S_matrix = hcat(S_blocks[1:Kblocks]...)   # D×(Kblocks*L)

    # Project S into truncated subspace:
    # Z = S_matrix * VH1   (D×m'), then X = Z * Σ1^{-1} * y_j
    ##Z = S_matrix * VH1    # D×m'
    W_orig = VH1 * Y    # (Kblocks*L) × m'

    #X_hat = Matrix{ComplexF64}(undef, D, mprime)
    #for j in 1:mprime
    #    xj = Z * (Σ1 \ Y[:, j])
    #    n = norm(xj)
    #    X_hat[:, j] .= n > 0 ? xj ./ n : xj
    #end
    X_hat = Matrix{ComplexF64}(undef, D, mprime)
    for j in 1:mprime
        xj = S_matrix * W_orig[:, j]   # x_j ≈ S * w_j
        n  = norm(xj)
        X_hat[:, j] .= n > 0 ? xj ./ n : xj
    end   

    # --- Step 7: residual check ---------------------------------------------
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


#=
# ============================================================

# Override / redefine the complex-γ Helmholtz 2D kernelvals

function kernelvals(biop::HelmholtzOperator2D{T, K}, tgeo, bgeo) where {T, K <: Complex}

    # Local numeric gamma (avoid BEAST.gamma function!)
    gamma = biop.gamma

    # Decide wavenumber k from gamma = i k  (so k = -i*gamma)
    # Optimised real case if gamma is purely imaginary
    if iszero(real(gamma))
        k = imag(gamma)           # gamma = i k ⇒ k real
    else
        k = -im * gamma           # general complex gamma
    end

    r = tgeo.cart - bgeo.cart
    R = norm(r)

    kr = k * R
    hankels = hankelh2.([0, 1], kr)
    green = -im / 4 * hankels[1]

    # Gradient wrt observation point
    gradgreen = k * im / 4 * hankels[2] * r / R

    txty = dot(normal(tgeo), normal(bgeo))

    return KernelValsHelmholtz2D(gamma, k, r, R, green, gradgreen, txty)
end
=#

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



# analytical value: k01 = 2.4048255577

"""
    assemble_T_circle_Dirichlet(k, circle, X)

Build the boundary operator T(k) = -0.5 I + K_k for the
2D interior Dirichlet Helmholtz problem on the circle.
"""
#function assemble_T_circle_Dirichlet(k::Number,
#                                     circle::circle,
#                                     X)
#    𝒟 = Helmholtz2D.doublelayer(; wavenumber = k)
#    D  = assemble(𝒟, X, X)
#    I  = assemble(BEAST.Identity(), X, X)
#
#    # interior Dirichlet: T(k) = -0.5 I + K_k
#    M  = -0.5 .* I .+ D
#    return Matrix(M)
#end

function assemble_T_circle_Dirichlet(k::Number, X)
    𝒟 = Helmholtz2D.doublelayer(; wavenumber = k)
    D  = assemble(𝒟, X, X)
    I  = assemble(BEAST.Identity(), X, X)
    return Matrix(-0.5 .* I .+ D)
end


"""
    Beyn_circle_D(k, VHat, MInvVHat, circle, X)

User function for Beyn: solves T(k) * X = VHat.
"""
#function Beyn_circle_D(k::ComplexF64,
#                       VHat::Matrix{ComplexF64},
#                       MInvVHat::Matrix{ComplexF64},
#                       circle,
#                       X)
#    T = assemble_T_circle_Dirichlet(k, circle, X)
#    MInvVHat .= T \ VHat   # solve
#    return MInvVHat
#end

function Beyn_circle_D(k, VHat, MInvVHat, X)
    T = assemble_T_circle_Dirichlet(k, X)
    MInvVHat .= T \ VHat
    return MInvVHat
end

"""
    ApplyT_circle_D(k, v, circle, X)

Apply T(k) to vector v.
"""
#function ApplyT_circle_D(k::ComplexF64,
#                         v::AbstractVecOrMat{ComplexF64},
#                         circle,
#                         X)
#    #T = assemble_T_circle_Dirichlet(k, circle, X)
#    T = assemble_T_circle_Dirichlet(k, X)
#    return T * v
#end
function ApplyT_circle_D(k::ComplexF64,
    v::AbstractVecOrMat{ComplexF64},
    X)
    T = assemble_T_circle_Dirichlet(k, X)
    return T * v
end


### dal FMM
#=

function Beyn_circle_D(k, VHat, MInvVHat, X)
    T = assemble_T_circle_Dirichlet(k, X)  # <-- składanie macierzy (gęste!)
    MInvVHat .= T \ VHat                   # <-- solve na macierzy
    return MInvVHat
end

function UserBeynFunction_FMM(k, VHat, MInvVHat, cache; tol=1e-10, maxiter=200)
    D = size(VHat,1)

    A = LinearMap{ComplexF64}((y,x)->(y .= ApplyT_FMM(k, x, cache)), D, D)

    # solve kolumna po kolumnie (na start) — potem można block-GMRES
    for j in 1:size(VHat,2)
        rhs = view(VHat,:,j)
        sol = view(MInvVHat,:,j)
        fill!(sol, 0)
        gmres!(sol, A, rhs; tol=tol, maxiter=maxiter, log=false)
    end
    return MInvVHat
end


# to ma liczyć y = T(k) * x  (FMM w środku)
function ApplyT_FMM(k, x, cache)
    # 1) near-field (lokalne) + 2) far-field FMM (kernel zależny od k)
    # y = Near(k)*x + FMM_apply(k, x, cache)
    return y
end

UserFun = (k, Vhat, MIV) -> Beyn_circle_D(k, Vhat, MIV, X)
ApplyT  = (k, v)         -> ApplyT_circle_D(k, v, X)


cache = ... # precomputy FMM: drzewo, grupy, geometrię, near-field itp.

UserFun = (k, Vhat, MIV) -> UserBeynFunction_FMM(k, Vhat, MIV, cache)
ApplyT  = (k, v)         -> ApplyT_FMM(k, v, cache)


 function Beyn_circle_D(k, VHat, MInvVHat, X;
                       cache=nothing,
                       tol::Float64=1e-10,
                       maxiter::Int=300,
                       verbose::Bool=false)

    D, L = size(VHat)
    @assert size(MInvVHat) == (D, L)

    # --- operator A: x -> T(k)*x (tu ma siedzieć FMM) -----------------------
    function apply_T!(y, x)
        # y .= ApplyT_FMM(k, x, X, cache)
        # >>> WSTAW SWÓJ KOD <<<

        # Minimalny placeholder (zabroniony docelowo): użycie macierzy:
        # T = assemble_T_circle_Dirichlet(k, X)
        # mul!(y, T, x)

        error("apply_T! not implemented: tutaj musi być ApplyT_FMM / matrix-free matvec")
    end

    A = LinearMap{ComplexF64}(apply_T!, D, D; ismutating=true)

    # --- solve: każda kolumna RHS osobno (start) ----------------------------
    for j in 1:L
        rhs = view(VHat, :, j)
        sol = view(MInvVHat, :, j)
        fill!(sol, 0)

        it = gmres!(sol, A, rhs; tol=tol, maxiter=maxiter, log=false)
        if verbose
            # it to może być liczba iteracji / obiekt zależnie od wersji IterativeSolvers
            @info "GMRES done for rhs=$j"
        end
    end

    return MInvVHat
end   
=#








function test_circle_cavity_Dirichlet_Beyn()
    println("=== Beyn eigenvalue test: circular cavity (Dirichlet, interior) ===")

    # Geometry & discretization
    a       = 1.0                      # radius
    p_order = 2                        # polynomial order of boundary elements
    h       = 2π * a / 64              # mesh size along boundary

    #circle  = circle_curvilinear(a, p_order; h = h)
    circle = CompScienceMeshes.meshcircle(a, h)

    # Use your C0 scalar basis on the boundary
    # (If you prefer the "d1 linear" variant, swap the next line)
    #X = lagrangec0_curvilinear(circle)
    #X = lagrangec0d1_curvilinear_exact_linear(circle)
    X = lagrangecx(circle, order=p_order)

    # Dimension D from actual matrix:
    #T0 = assemble_T_circle_Dirichlet(2.0, circle, X)  # dummy k
    T0 = assemble_T_circle_Dirichlet(2.0, X)
    D  = size(T0, 1)
    @printf("D = %d\n", D)

    # Analytic first Dirichlet eigenvalue (n=0, m=1): k = j_{0,1}/a
    #j01 = besseljzero(0, 1)
    j01 = besselj_zero(0, 1)
    k_true = j01 / a
    @printf("Analytic k_{0,1} ≈ %.10f\n", k_true)

    # Contour in k-plane: circle around k_true
    z0 = k_true + 0.0im
    R  = 0.4        # small enough to (hopefully) contain only k_{0,1}

    # Beyn parameters
    L       = 8     # probing subspace dimension (≥ expected # of eigenvalues)
    N       = 32    # quadrature points on contour
    Kmax    = 2
    tol_svd = 1e-8
    tol_res = 1e-8
    verbose = true

    # Closures capturing circle and X
    #UserFun = (k, Vhat, MIV) -> Beyn_circle_D(k, Vhat, MIV, circle, X)
    #ApplyT  = (k, v)         -> ApplyT_circle_D(k, v, circle, X)
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

    # Sort by real part for nicer printing
    p = sortperm(λ, by = x -> (real(x), imag(x)))
    λ = λ[p]

    println("Eigenvalues k inside contour:")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e, %+.12e   |k - k_true| ≈ %.3e\n",
                j, real(val), imag(val), abs(val - k_true))
    end

    return λ, V, ok
end

λ, V, ok = test_circle_cavity_Dirichlet_Beyn()


function test_circle_cavity_Dirichlet_Asakura(;
        a::Float64      = 1.0,            # radius
        p_order::Int    = 2,              # polynomial order on boundary
        h::Float64      = 2π * a / 64,    # mesh size along boundary
        L::Int          = 8,              # block size
        Kblocks::Int    = 2,              # number of block moments
        N::Int          = 64,             # quadrature points on contour
        δ_svd::Float64  = 1e-12,
        tol_res::Float64 = 1e-8,
        verbose::Bool   = true)

    println("=== Asakura eigenvalue test: circular cavity (Dirichlet, interior) ===")

    # --- Geometry and boundary space (same as in Beyn test) -----------------
    #circle = circle_curvilinear(a, p_order; h = h)
    circle = CompScienceMeshes.meshcircle(a, h)

    # Use the same basis as in test_circle_cavity_Dirichlet_Beyn
    #X = lagrangec0d1_curvilinear_exact_linear(circle)
    #order = 1 #[2; 3]
    X = lagrangecx(circle, order=p_order)

    # Dimension D from an assembled operator
    #T0 = assemble_T_circle_Dirichlet(2.0, circle, X)  # dummy k
    T0 = assemble_T_circle_Dirichlet(2.0, X)
    D  = size(T0, 1)
    @printf("D = %d\n", D)

    # Analytic first Dirichlet eigenvalue: k_{0,1} = j_{0,1} / a
    j01    = besselj_zero(0, 1)
    k_true = j01 / a
    @printf("Analytic k_{0,1} ≈ %.10f\n", k_true)

    # --- Contour in k-plane (same philosophy as Beyn test) ------------------
    z0 = k_true + 0.0im    # center near the first eigenvalue
    R  = 0.4               # small radius, should contain only k_{0,1}

    # --- Closures capturing circle and X ------------------------------------
    # UserFun: T(k)^{-1} * Vhat
    #UserFun = (k, Vhat, MInvVHat) -> Beyn_circle_D(k, Vhat, MInvVHat, circle, X)
    UserFun = (k, Vhat, MInvVHat) -> Beyn_circle_D(k, Vhat, MInvVHat, X)

    # ApplyT: T(k) * v
    #ApplyT  = (k, v) -> ApplyT_circle_D(k, v, circle, X)
    ApplyT  = (k, v) -> ApplyT_circle_D(k, v, X)

    # --- Call Asakura solver -------------------------------------------------
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

    # Sort by real part for nicer printing
    p = sortperm(λ, by = x -> (real(x), imag(x)))
    λ = λ[p]

    println("Eigenvalues k inside contour (Asakura):")
    for (j, val) in enumerate(λ)
        @printf("%2d: k ≈ %+.12e, %+.12e   |k - k_true| ≈ %.3e\n",
                j, real(val), imag(val), abs(val - k_true))
    end

    return λ, V, ok
end


λA, XA, okA = test_circle_cavity_Dirichlet_Asakura()


