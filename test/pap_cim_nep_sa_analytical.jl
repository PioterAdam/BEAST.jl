
using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test
using SpecialFunctions
using Gmsh


import CompScienceMeshes: measure, nodes, Simplex
import LinearAlgebra: norm
import CompScienceMeshes: indices, CurvilinearMesh, chart

const Shape = BEAST.Shape
const LagrangeBasis = BEAST.LagrangeBasis

#using LinearAlgebra
using Printf
using Base.Threads

# in your version of BEAST probably not required
# measure for straight 1D edges embedded in R^U
function measure(s::Simplex{U,1,<:Any,<:Any,T}) where {U,T}
    p0, p1 = nodes(s)
    return norm(p1 - p0)
end

#=
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

        return CurvilinearMesh(verts, faces, porder)
    finally
        gmsh.finalize()
    end
end

# in your version of BEAST probably not required
function lagrangec0d1_curvilinear_exact_linear(mesh::CurvilinearMesh)
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

# ============================================================

# Examples: analytical benchmarks 

# 1) Single place where we define the matrix T(z)
function F_example2(z)
    F = Matrix{ComplexF64}(undef, 3, 3)

    F[1,1] = z + 1
    F[1,2] = 6*z^2 - 6*z
    F[1,3] = 0

    F[2,1] = 2*z
    F[2,2] = 6*z^2 - 7*z + 1
    F[2,3] = 0

    F[3,1] = 0
    F[3,2] = 0
    F[3,3] = z^2 + 1

    return F
end

# Beyns interfaces
# 2) User function for Beyn: T(z)^{-1} * VHat
function Beyn_Example2(z, VHat, MInvVHat)
    F = F_example2(z)
    MInvVHat .= F \ VHat
    return MInvVHat
end

# 3) ApplyT for residual check: T(z) * v
function ApplyT_Example2(z, v)
    F = F_example2(z)
    return F * v
end

D = 3
λ, V, ok = BeynSolve_auto(0.0 + 0im, 1.5,
                          Beyn_Example2, ApplyT_Example2, D;
                          L = 3, N = 64, Kmax = 3,
                          tol_svd = 1e-10, tol_res = 1e-10,
                          verbose = true)

println("ok = ", ok)
println("λ = ", λ)


# Build M(z) for Example 4.11 in-place
function build_M_Beyn411!(M, z)
    D = size(M, 1)
    DE  =  2.0*D - 4.0*z/(6.0*D)   # diagonal entry
    ODE = -1.0*D -     z/(6.0*D)   # off-diagonal entry

    fill!(M, 0)                    # zero everything

    # top row
    M[1,1] = DE
    M[1,2] = ODE

    # interior rows 2..D-1
    for d = 2:D-1
        M[d,d]   = DE
        M[d,d-1] = ODE
        M[d,d+1] = ODE
    end

    # bottom row
    M[D,D-1] = ODE
    M[D,D]   = 0.5*DE + z/(z-1.0)

    return M
end

# T(z)^{-1} * VHat  (used by Beyn)
function Beyn411_user(z, VHat, MInvVHat, M)
    build_M_Beyn411!(M, z)
    MInvVHat .= M \ VHat
    return MInvVHat
end

# T(z) * v (used for residual checks)
function ApplyT_Beyn411(z, v, M)
    build_M_Beyn411!(M, z)
    return M * v
end

function SolveBeyn411_auto(; z0 = 150.0 + 0.0im,
                             Rx = 148.0,
                             Ry = 148.0,
                             L  = 10,
                             N  = 50,
                             Kmax = 2,
                             tol_svd = 1e-8,
                             tol_res = 1e-8)

    D = 400
    M = zeros(ComplexF64, D, D)  # workspace reused every call

    # closures capturing M
    UserFun = (z, VHat, MInvVHat) -> Beyn411_user(z, VHat, MInvVHat, M)
    ApplyT  = (z, v)             -> ApplyT_Beyn411(z, v, M)

    λ, V, ok = BeynSolve_auto(z0, Rx, Ry,
                          (z,Vhat,MIV)->Beyn411_user(z,Vhat,MIV,M),
                          (z,v)->ApplyT_Beyn411(z,v,M),
                          D;
                          L=10, N=50, Kmax=2,
                          tol_svd=1e-8, tol_res=1e-8,
                          verbose=true)


    @printf("ok = %s\n", ok ? "true" : "false")
    for n in 1:length(λ)
        @printf("%2d: {%+.12e, %+.12e}\n",
                n, real(λ[n]), imag(λ[n]))
    end

    return λ, V, ok
end

λ, V, ok = SolveBeyn411_auto()

println("ok = ", ok)
println("λ = ", λ)



#Asakura interfaces

UserFun = (z, VHat, MInvVHat) -> Beyn411_user(z, VHat, MInvVHat, M)
ApplyT  = (z, v)              -> ApplyT_Beyn411(z, v, M)


function SolveBeyn411_Asakura_circle(; z0     = 150.0 + 0.0im,
                                      R      = 148.0,
                                      L      = 8,
                                      Kblocks = 4,
                                      N      = 64,
                                      δ_svd  = 1e-12,
                                      tol_res = 1e-8,
                                      verbose = true)

    D = 400
    M = zeros(ComplexF64, D, D)

    UserFun = (z, VHat, MInvVHat) -> Beyn411_user(z, VHat, MInvVHat, M)
    ApplyT  = (z, v)              -> ApplyT_Beyn411(z, v, M)

    λ, V, ok = AsakuraSolve_circle_auto(z0, R,
                                        UserFun, ApplyT, D;
                                        L        = L,
                                        Kblocks  = Kblocks,
                                        N        = N,
                                        δ_svd    = δ_svd,
                                        tol_res  = tol_res,
                                        verbose  = verbose)

    @printf("ok = %s\n", ok ? "true" : "false")
    for (n, λn) in enumerate(λ)
        @printf("%2d: {%+.12e, %+.12e}\n", n, real(λn), imag(λn))
    end

    return λ, V, ok
end


λ_a, V_a, ok_a = SolveBeyn411_Asakura_circle()


println("ok = ", ok_a)
println("λ = ", λ_a)

# 2) UserFun: F(z) \ VHat  (block solve)
function Asakura_Example2_UserFun(z, VHat, MInvVHat)
    F = F_example2(z)
    MInvVHat .= F \ VHat
    return MInvVHat
end

# 3) ApplyT: F(z) * v  (for residual check)
function Asakura_Example2_ApplyT(z, v)
    F = F_example2(z)
    return F * v
end

# 4) Wrapper calling AsakuraSolve_circle_auto for Example 2
function Solve_Asakura_Example2(; 
        z0      = 0.0 + 0.0im,   # center of contour
        R       = 1.2,           # radius (covers -i, i, 1/3, 1/2, 1)
        L       = 3,             # block size (≥ expected multiplicity)
        Kblocks = 2,             # number of block moments
        N       = 64,            # quadrature points on circle
        δ_svd   = 1e-12,
        tol_res = 1e-10,
        verbose = true)

    D = 3  # dimension of F(z)

    # closures in the Asakura interface format
    UserFun = (z, VHat, MInvVHat) -> Asakura_Example2_UserFun(z, VHat, MInvVHat)
    ApplyT  = (z, v)              -> Asakura_Example2_ApplyT(z, v)

    λ, X, ok = AsakuraSolve_circle_auto(z0, R,
                                        UserFun, ApplyT, D;
                                        L        = L,
                                        Kblocks  = Kblocks,
                                        N        = N,
                                        δ_svd    = δ_svd,
                                        tol_res  = tol_res,
                                        verbose  = verbose)

    @printf("ok = %s\n", ok ? "true" : "false")
    for (j, λj) in enumerate(λ)
        @printf("%2d: λ[%d] = % .15e %+.15ei\n",
                j, j, real(λj), imag(λj))
    end

    return λ, X, ok
end

λ_asakura, X_asakura, ok_asakura = Solve_Asakura_Example2()


println("ok = ", ok_asakura)
println("λ = ", λ_asakura)
