using LinearAlgebra
using StaticArrays
using CompScienceMeshes
using BEAST
import BEAST: Shape, LagrangeBasis

# ---------- robust endpoint clustering (handles duplicated Gmsh IDs) ----------
@inline function _bbox_diag(mesh)
    U = universedimension(mesh)
    T = coordtype(mesh)
    lo = fill(typemax(T), U)
    hi = fill(typemin(T), U)
    @inbounds for p in mesh.vertices
        @inbounds for i in 1:U
            x = p[i]
            x < lo[i] && (lo[i] = x)
            x > hi[i] && (hi[i] = x)
        end
    end
    return norm(SVector{U,T}(hi) - SVector{U,T}(lo))
end

function _endpoint_clusters(mesh::CompScienceMeshes.CurvilinearMesh{U,N,T,O};
                            rtol::T = T(1e-6), atol::T = T(0)) where {U,N,T,O}
    reps    = SVector{U,T}[]
    members = Vector{Vector{Tuple{Int,Int}}}()
    L = _bbox_diag(mesh)
    tol = max(atol, rtol * (L > 0 ? L : one(T)))

    @inbounds for (c, face) in enumerate(mesh.faces)
        Nloc = length(face)
        @inbounds for r in (1, Nloc)             # only endpoints
            q = mesh.vertices[face[r]]
            idx = 0
            @inbounds for i in eachindex(reps)
                if norm(q - reps[i]) ≤ tol
                    idx = i; break
                end
            end
            if idx == 0
                push!(reps, q)
                push!(members, Tuple{Int,Int}[(c, r)])
            else
                push!(members[idx], (c, r))
            end
        end
    end
    return reps, members
end

function _add_endpoint_dofs!(fns::Vector{<:Vector}, pos::Vector,
                             mesh::CompScienceMeshes.CurvilinearMesh{U,N,T,O};
                             dirichlet::Bool = true,
                             rtol::T = T(1e-6), atol::T = T(0)) where {U,N,T,O}
    S = Shape{T}
    reps, members = _endpoint_clusters(mesh; rtol=rtol, atol=atol)
    keepdeg = dirichlet ? 2 : 1
    @inbounds for i in eachindex(reps)
        length(members[i]) < keepdeg && continue   # drop degree-1 ends on open chains
        shapes = S[]
        @inbounds for (c, r) in members[i]
            push!(shapes, S(c, r, one(T)))
        end
        push!(fns, shapes)
        push!(pos, reps[i])
    end
    return nothing
end

# ------------------------ main basis constructor ------------------------
function lagrangec0_curvilinear_(mesh::CompScienceMeshes.CurvilinearMesh{U,N,T,O};
                                order::Integer = O,
                                dirichlet::Bool = true,
                                rtol::T = T(1e-6), atol::T = T(0)) where {U,N,T,O}
    Nloc = length(first(mesh.faces))
    @assert all(length(f) == Nloc for f in mesh.faces)
    @assert order == Nloc - 1  # geometry order must match basis order

    S = Shape{T}
    P = SVector{U,T}

    fns = Vector{Vector{S}}()
    pos = Vector{P}()

    # 1) share endpoints globally (geometric clusters)
    _add_endpoint_dofs!(fns, pos, mesh; dirichlet=dirichlet, rtol=rtol, atol=atol)

    # 2) keep interior nodes element-local (r = 2..Nloc-1)
    if Nloc > 2
        @inbounds for (c, face) in enumerate(mesh.faces)
            @inbounds for r in 2:(Nloc-1)
                push!(fns, S[S(c, r, one(T))])
                push!(pos, mesh.vertices[face[r]])
            end
        end
    end

    return LagrangeBasis{order,0,Nloc}(mesh, fns, pos)
end

function lagrangec0_curvilinear(mesh::CompScienceMeshes.CurvilinearMesh{U,N,T,O};
                                order::Integer = O,
                                dirichlet::Bool = true,
                                rtol::T = T(1e-4), atol::T = T(0)) where {U,N,T,O}
    Nloc = order + 1
    S = BEAST.Shape{T}
    P = SVector{U,T}
    fns = Vector{Vector{S}}()
    pos = Vector{P}()

    # 1) endpoints (shared) by geometric clustering
    _add_endpoint_dofs!(fns, pos, mesh; dirichlet=dirichlet, rtol=rtol, atol=atol)

    # 2) interior DoFs (element-local), independent of faces length
    if Nloc > 2
        S = BEAST.Shape{coordtype(mesh)}
        @inbounds for (c, face) in enumerate(mesh.faces)
            if length(face) == Nloc
                for r in 2:(Nloc-1)
                    push!(fns, S[S(c, r, one(coordtype(mesh)))])
                    push!(pos, mesh.vertices[face[r]])   # exact control point
                end
            else
                # fallback if faces lack interiors:
                ch = chart(mesh, face)
                for r in 2:(Nloc-1)
                    push!(fns, S[S(c, r, one(coordtype(mesh)))])
                    push!(pos, cartesian(center(ch)))
                end
            end
        end
    end

    return BEAST.LagrangeBasis{order,0,Nloc}(mesh, fns, pos)
end


# expected DoFs for testing (uses the same clustering)
function expected_dofs_curvilinear(mesh::CompScienceMeshes.CurvilinearMesh{U,N,T,O};
                                   dirichlet::Bool = true,
                                   rtol::T = T(1e-6), atol::T = T(0)) where {U,N,T,O}
    _, members = _endpoint_clusters(mesh; rtol=rtol, atol=atol)
    keepdeg = dirichlet ? 2 : 1
    n_end = count(m -> length(m) ≥ keepdeg, members)
    Nloc = length(first(mesh.faces))
    E    = length(mesh.faces)
    return n_end + E * max(Nloc - 2, 0)
end



"""
    circle_curvilinear(radius, porder; h = 2π*radius/64)

Return a 1D-in-2D CurvilinearMesh of a circle with polynomial order `porder`
(using Gmsh high-order line elements). `h` controls target edge length.
"""
function circle_curvilinear(radius::Real, porder::Integer; h::Real = 2π*radius/64)
    @assert porder ≥ 1 "porder must be ≥ 1"
    gmsh.initialize()
    try
        gmsh.model.add("circle_p$(porder)")
        s = gmsh.model.occ.addDisk(0.0, 0.0, 0.0, radius, radius)
        gmsh.model.occ.synchronize()

        # sizing + high-order
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", h)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", h)
        gmsh.option.setNumber("Mesh.ElementOrder", porder)
        gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)
        gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)

        # physicals
        gmsh.model.addPhysicalGroup(2, [s], 1)
        gmsh.model.setPhysicalName(2, 1, "Domain")
        curves = [t[2] for t in gmsh.model.getBoundary([(2, s)], false, false, false) if t[1] == 1]
        gmsh.model.addPhysicalGroup(1, curves, 2)
        gmsh.model.setPhysicalName(1, 2, "Boundary")

        gmsh.model.mesh.generate(2)

        # ---- nodes on boundary physical (robust to 2- or 3-value return) ----
        nb = gmsh.model.mesh.getNodesForPhysicalGroup(1, 2)
        nodeTags   = nb[1]
        nodeCoords = nb[2]
        @assert length(nodeCoords) == 3*length(nodeTags) "unexpected coords size"

        # map tag -> local index
        tag2idx = Dict{Int,Int}(Int(nodeTags[i]) => i for i in eachindex(nodeTags))

        # vertices in any order (faces will reference via tag2idx)
        verts = Vector{SVector{2,Float64}}(undef, length(nodeTags))
        @inbounds for i in eachindex(nodeTags)
            x = nodeCoords[3(i-1)+1]; y = nodeCoords[3(i-1)+2]
            verts[i] = SVector{2,Float64}(x, y)
        end

        # ---- faces from boundary curves, only keep the type matching porder ----
        NF = porder + 1
        faces = SVector{NF,Int}[]
        for c in curves
            types, elemTags, elemNodeTags = gmsh.model.mesh.getElements(1, c)
            @inbounds for k in eachindex(types)
                nodes = elemNodeTags[k]
                # keep blocks whose stride matches NF (this selects the right line type)
                if length(nodes) % NF != 0; continue; end
                ne = div(length(nodes), NF)
                for e in 1:ne
                    ids = ntuple(j -> tag2idx[Int(nodes[NF*(e-1)+j])], NF)
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


using Test
using StaticArrays
using CompScienceMeshes: indices
using BEAST: Shape, LagrangeBasis
using Gmsh

import CompScienceMeshes: indices

using CompScienceMeshes
using StaticArrays


# -------- 0D vertex mesh for CurvilinearMesh skeleton --------

struct CurviVertexMesh{U,T} <: CompScienceMeshes.AbstractMesh{U,0,T}
    vertices::Vector{SVector{U,T}}
end

# minimal mesh interface
CompScienceMeshes.dimension(::CurviVertexMesh) = 0
CompScienceMeshes.universedimension(::CurviVertexMesh{U}) where {U} = U
CompScienceMeshes.coordtype(::CurviVertexMesh{U,T}) where {U,T} = T
CompScienceMeshes.vertextype(::CurviVertexMesh{U,T}) where {U,T} = SVector{U,T}
CompScienceMeshes.numvertices(m::CurviVertexMesh) = length(m.vertices)
CompScienceMeshes.cells(m::CurviVertexMesh) = Base.OneTo(length(m.vertices))
CompScienceMeshes.cell(m::CurviVertexMesh, i::Int) = i  # vertex "cell" = its index

# charts for vertices
struct CurviVertexChart{U,T}
    p::SVector{U,T}
end
CompScienceMeshes.chart(m::CurviVertexMesh{U,T}, v::Int) where {U,T} = CurviVertexChart{U,T}(m.vertices[v])
CompScienceMeshes.center(ch::CurviVertexChart) = ch.p
CompScienceMeshes.cartesian(p::SVector) = p  # passthrough (safe even if CSM defines it)

# -------- skeleton overrides (avoid skeleton_fast/celltype) --------
import CompScienceMeshes: skeleton
skeleton(m::CompScienceMeshes.CurvilinearMesh{U,N,T,O}, dim::Int) where {U,N,T,O} =
    dim == 1 ? m :
    dim == 0 ? CurviVertexMesh{U,T}(m.vertices) :
    throw(ArgumentError("skeleton dimension $dim out of range for 1D mesh"))

# optional Val{…} overloads if your code uses them
skeleton(m::CompScienceMeshes.CurvilinearMesh{U,N,T,O}, ::Type{Val{1}}) where {U,N,T,O} = m
skeleton(m::CompScienceMeshes.CurvilinearMesh{U,N,T,O}, ::Type{Val{0}}) where {U,N,T,O} = CurviVertexMesh{U,T}(m.vertices)

## -------- indices helpers (unblock vertextocellmap/connectivity if needed) --------
##import CompScienceMeshes: indices
##indices(m::CompScienceMeshes.CurvilinearMesh{U,N,T,O}, i::Int) where {U,N,T,O} = m.faces[i]
##indices(::CompScienceMeshes.CurvilinearMesh{U,N,T,O}, cell::SVector{N,Int}) where {U,N,T,O} = cell

import CompScienceMeshes: cells, cell, indices

# Cell iteration/helpers (safe shims)
cells(m::CurvilinearMesh) = Base.OneTo(length(m.faces))
cell(m::CurvilinearMesh, i::Int) = m.faces[i]

# Connectivity queries used by vertextocellmap:
#  1) by cell id (Int) -> return the SVector of vertex ids
indices(m::CurvilinearMesh{U,N,T,O}, i::Int) where {U,N,T,O} = m.faces[i]
#  2) by cell object (SVector{N,Int}) -> identity for skeleton
##indices(::CurvilinearMesh{U,N,T,O}, c::SVector{N,Int}) where {U,N,T,O} = c


@testset "Curvi C0 Lagrange (counts on circle)" begin
    R = 1.0
    E = 8
    for p in (1,2,3,4)
        Γ = circle_curvilinear(R, p; h = 2π*R/E)   # geometric order = p (Nloc = p+1)
        X = lagrangec0_curvilinear(Γ; order=p, dirichlet=true)  # closed loop

        pos = getfield(X, :pos)
        expected = expected_dofs_curvilinear(Γ; dirichlet=true)

        @info "p=$p" E=E Nloc=length(first(Γ.faces)) length_pos=length(pos) expected
        @test length(pos) == expected
    end
end

using Test
using LinearAlgebra
using StaticArrays
using CompScienceMeshes
using BEAST

# --- helpers (same as before, but returns tol so we reuse it consistently) ---
function _endpoint_clusters_test(mesh; rtol=1e-4, atol=0.0)
    U = universedimension(mesh); T = coordtype(mesh)
    lo = fill(typemax(T), U); hi = fill(typemin(T), U)
    @inbounds for p in mesh.vertices, i in 1:U
        x = p[i]; x < lo[i] && (lo[i] = x); x > hi[i] && (hi[i] = x)
    end
    L = norm(SVector{U,T}(hi) - SVector{U,T}(lo))
    tol = max(T(atol), T(rtol) * (L > 0 ? L : one(T)))

    reps    = SVector{U,T}[]
    members = Vector{Vector{Tuple{Int,Int}}}()
    @inbounds for (c, face) in enumerate(mesh.faces)
        Nloc = length(face)
        for r in (1, Nloc)
            q = mesh.vertices[face[r]]
            idx = 0
            for i in eachindex(reps)
                if norm(q - reps[i]) ≤ tol
                    idx = i; break
                end
            end
            if idx == 0
                push!(reps, q); push!(members, Tuple{Int,Int}[(c, r)])
            else
                push!(members[idx], (c, r))
            end
        end
    end
    return reps, members, tol
end

_contains_coord(pos::Vector, q; atol=1e-10) = any(p -> norm(p .- q) ≤ atol, pos)

# -------- test --------
@testset "lagrangec0_curvilinear — counts & positions on a circle" begin
    R = 1.0
    Nseg = 8
    dirichlet = true
    RTOL = 1e-4   # use the SAME tolerance in basis and test

    for p in (1,2,3,4)
        Γ = circle_curvilinear(R, p; h = 2π*R/Nseg)

        # build basis with the same RTOL
        X = lagrangec0_curvilinear(Γ; order=p, dirichlet=dirichlet, rtol=RTOL)

        pos   = getfield(X, :pos)
        E     = length(Γ.faces)
        Nloc  = length(first(Γ.faces))
        reps, members, tol = _endpoint_clusters_test(Γ; rtol=RTOL)
        keepdeg = dirichlet ? 2 : 1
        n_end = count(m -> length(m) ≥ keepdeg, members)
        expected_ndofs = n_end + E * max(Nloc - 2, 0)

        @info "p=$p" E Nloc length_pos=length(pos) expected_ndofs
        @test length(pos) == expected_ndofs

        # endpoint DoFs present (only check kept clusters)
        for (i, q) in pairs(reps)
            length(members[i]) < keepdeg && continue
            @test _contains_coord(pos, q; atol=tol)
        end

        # interior DoFs present
        if Nloc > 2
            # collect interior control points from geometry
            interiors = SVector{universedimension(Γ), coordtype(Γ)}[]
            @inbounds for face in Γ.faces, r in 2:(Nloc-1)
                push!(interiors, Γ.vertices[face[r]])
            end
            @test length(interiors) == E * (Nloc - 2)
            for q in interiors
                @test _contains_coord(pos, q; atol=1e-10)
            end
        end
    end
end

using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test
using Plots

# ---- 2D manufactured test on a curvilinear circle ----
hs     = [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625]  # non-dim arc step
#orders = [1, 2, 3, 4]   # isoparametric: geometry order = basis order
orders = [1, 2]   # isoparametric: geometry order = basis order

accs   = Dict(o => Float64[] for o in orders)

for order in orders
    for h in hs
        # geometry & k
        R = 10.0
        λ = 20R
        k = 2π / λ

        # choose number of segments so arc step ≈ h*R
        Nseg = max(8, ceil(Int, 2π / h))

        # curvilinear circle (order = p)
        Γ = circle_curvilinear(R, order; h = 2π*R/Nseg)

        # IMPORTANT: use your curvilinear basis builder (avoid BEAST.lagrangec0 default)
        X = lagrangec0_curvilinear(Γ; order=order, dirichlet=true, rtol=1e-4)
        @show order numfunctions(X)

        E = length(Γ.faces)                 # number of edges
        @show order E numfunctions(X)
        @assert numfunctions(X) == E*order  "Expected $(E*order) DoFs for p=$order; got $(numfunctions(X))"

        # operators & incident field (two distant monopoles)
        S = Helmholtz2D.singlelayer(; gamma = im*k)
        q, ϵ = 100.0, 1.0
        pos1, pos2 = SVector(R*30.0, 0.0), SVector(-R*30.0, 0.0)
        charge1 = Helmholtz2D.monopole(position=pos1, amplitude=q/(4π*ϵ), wavenumber=k)
        charge2 = Helmholtz2D.monopole(position=pos2, amplitude=-q/(4π*ϵ), wavenumber=k)
        Φ_inc(x) = charge1(x) + charge2(x)

        # Dirichlet data & solve Sρ = -gD
        gD = assemble(DirichletTrace(charge1), X) + assemble(DirichletTrace(charge2), X)
        @time M = assemble(S, X, X)
        ρ = M \ (-gD)

        # evaluate inside (total field ≈ 0)
        pts = meshcircle(0.6R, 0.24R).vertices
        pot_sc = potential(BEAST.HH2DSingleLayerNear(im*k), pts, ρ, X; type=ComplexF64)
        err = norm(pot_sc .+ Φ_inc.(pts)) / norm(Φ_inc.(pts))

        push!(accs[order], err)
    end
end

# ---- plot ----
plotlyjs()
plt = plot(
    xscale=:log10, yscale=:log10, linewidth=3,
    xlabel="h (non-dimensional arc step)", ylabel="Relative error",
    legend=:bottomright, title="Manufactured solution — circle (Helmholtz2D, curvilinear)"
)
for o in orders
    plot!(plt, hs, accs[o], label="Order $o", markershape=:auto)
end
savefig("curvi_circle_relative_error.pdf")
display(plt)
