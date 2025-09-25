using Test
using LinearAlgebra

using BEAST, CompScienceMeshes, StaticArrays

using Gmsh

"""
    lagrangec0_gmsh_exact(mesh; order)

Build a C^0 Lagrange basis on a 1D (line) mesh using the **exact Gmsh nodes** as DOFs.
- One DOF per global mesh node (endpoints and all mid-edge nodes).
- For each global node g, we collect all incident (cell, local_node) pairs and
  attach the corresponding local shape `Shape(c, r, 1)` to that DOF.
- DOF positions are **exactly** `mesh.vertices[g]` (no reparametrization).

Assumes each line cell has `order+1` local nodes in Gmsh ordering.
"""
function lagrangec0_gmsh_exact(mesh; order::Int)
    @assert order > 0 "order must be positive for high-order lines"

    # Types
    T = coordtype(mesh)
    P = vertextype(mesh)
    S = BEAST.Shape{T}

    # Map: global node id -> vector of attached local shapes
    dofs = Dict{Int, Vector{S}}()

    # Traverse all line cells; attach every local node r = 1..NF
    for (c, cell) in pairs(mesh.faces)
        NF = length(cell)
        @assert NF == order + 1 "Cell $c has $NF nodes, expected $(order+1)"
        @inbounds for r in 1:NF
            g = cell[r]                          # global node id from Gmsh
            push!(get!(dofs, g, Vector{S}()), S(c, r, T(1)))
        end
    end

    # Build basis arrays, ordered by increasing global node id (deterministic)
    gids = sort!(collect(keys(dofs)))
    fns  = Vector{Vector{S}}(undef, length(gids))
    pos  = Vector{P}(undef, length(gids))
    @inbounds for (i, g) in enumerate(gids)
        fns[i] = dofs[g]                        # all shapes attached to node g
        pos[i] = mesh.vertices[g]               # exact Gmsh node position
    end

    NFcell = order + 1
    return BEAST.LagrangeBasis{order, 0, NFcell}(mesh, fns, pos)
end

gmsh.initialize()
try
    gmsh.model.add("disk_p2")
    s_tag = gmsh.model.occ.addDisk(0.0, 0.0, 0.0, 1.0, 1.0)
    gmsh.model.occ.synchronize()

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.10)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.10)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)

    gmsh.model.addPhysicalGroup(2, [s_tag], 1)
    gmsh.model.setPhysicalName(2, 1, "Domain")

    dimtags = [(Int32(2), Int32(s_tag))]
    bnd = gmsh.model.getBoundary(dimtags, false, false, false)
    curves = [t[2] for t in bnd if t[1] == 1]
    gmsh.model.addPhysicalGroup(1, curves, 2)
    gmsh.model.setPhysicalName(1, 2, "Boundary")

    gmsh.model.mesh.generate(2)
    mshpath = joinpath(@__DIR__, "assets", "circle2d_quadratic.msh")
    isdir(dirname(mshpath)) || mkpath(dirname(mshpath))
    gmsh.write(mshpath)
finally
    gmsh.finalize()
end

#m = load_gmsh_mesh(
#    joinpath(@__DIR__, "assets", "circle2d_quadratic.msh");
#    udim   = 2,
#    element = :line,
#    order   = 2,
#    physical = "Boundary",
#)

# Float32 not working since hankelh2 returns always F64
#for T in [Float64]
    #T = Float64
    #Γto = meshsegment(T(1.0), T(0.5))
    
    Γto = load_gmsh_mesh(
        joinpath(@__DIR__, "assets", "circle2d_quadratic.msh");
        udim   = 2,
        element = :line,
        order   = 2,
        physical = "Boundary",
    )
    
    Γt = Γto
    Γs = Γt

    Xt = lagrangecxd0(Γt)
    Xs = lagrangecxd0(Γs)

    λ = 10
    k = 2π/λ

    ops = [
        Helmholtz2D.singlelayer(; wavenumber=k)
    ]

    refstrat = BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30)

    for op in ops
        Sref = assemble(op, Xt, Xs; quadstrat=refstrat)
        ref = Sref[1, 1]

        n = 25

        #Sgl = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumQStrat(n, n+1))
        Sss = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumSauterQstrat(3,3,0,4,n,n))

        #gl = Sgl[1, 1]
        ss = Sss[1, 1]

        #glrel = norm(gl - ref) /  norm(ref)
        ssrel = norm(ss - ref) /  norm(ref)
        @show typeof(op), ssrel

        @test ssrel < 1e-14
    end

    #T = Float64
    Xt = lagrangec0_gmsh_exact(Γt; order=2) #lagrangec0d1(Γt)
    Xs = lagrangec0_gmsh_exact(Γs; order=2) #lagrangec0d1(Γs)

    ops = [
        Helmholtz2D.singlelayer(; wavenumber=k)
        Helmholtz2D.hypersingular(; wavenumber=k)
    ]
    𝒮 = 

    refstrat = BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30)

    for op in ops
        Sref = assemble(op, Xt, Xs; quadstrat=refstrat)
        ref = Sref[1, 1]

        n = 25

        #Sgl = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumQStrat(n, n+1))
        Sss = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumSauterQstrat(3,3,0,4,n,n))

        #gl = Sgl[1, 1]
        ss = Sss[1, 1]

        #glrel = norm(gl - ref) /  norm(ref)
        ssrel = norm(ss - ref) /  norm(ref)
        @show typeof(op), ssrel

        tol = op isa BEAST.HH2DHyperSingularFDBIO ? 2e-6 : 1e-11

        @test ssrel < tol
    end

    ## Test accuracy of vertex integrations
  
    T = Float64
    
    verts_x = SVector{2,T}[
        SVector(0.0, 0.0),  # end 1
        SVector(1.0, 0.0),  # end 2
        SVector(0.5, 0.0),  # midpoint
    ]
    faces_x = SVector{3,Int}[ SVector(1, 2, 3) ]  # Gmsh order: [end1, end2, mid]
    Γt = CurvilinearMesh(verts_x, faces_x, 2)
        
    verts_y = SVector{2,T}[
        SVector(0.0, 0.0),  # end 1
        SVector(0.0, 1.0),  # end 2
        SVector(0.0, 0.5),  # midpoint
    ]
    faces_y = SVector{3,Int}[ SVector(1, 2, 3) ]
    Γs = CurvilinearMesh(verts_y, faces_y, 2)
    
    Xt = lagrangecxd0(Γt)
    Xs = lagrangecxd0(Γs)
    
    #=
    g = BEAST.refspace(lagrangecxd0(Γt))
    qd = BEAST.quaddata(Helmholtz2D.singlelayer(wavenumber=1.0), g, g, [chart(Γt,1)], [chart(Γs,1)],
                    BEAST.DoubleNumSauterQstrat(10,10, 0,4, 10,10))
    r  = BEAST.quadrule(Helmholtz2D.singlelayer(wavenumber=1.0), g, g, 1, chart(Γt,1), 1, chart(Γs,1), qd,
                    BEAST.DoubleNumSauterQstrat(10,10, 0,4, 10,10))
    @test r isa BEAST.SauterSchwabQuadrature1D.CommonVertex


    τ = chart(Γt,1); σ = chart(Γs,1)
    tτ = CompScienceMeshes.tangents(CompScienceMeshes.neighborhood(τ, SVector(0.0)))[:,1]
    tσ = CompScienceMeshes.tangents(CompScienceMeshes.neighborhood(σ, SVector(0.0)))[:,1]
    nτ = SVector(-tτ[2], tτ[1]); nσ = SVector(-tσ[2], tσ[1])  # rotate CCW


    i  = first(cells(Γs))      # or any cell index
    ch = chart(Γs, i)          # <- this is the CurvilinearSimplex you printed
    ζs = (0.0, 0.37, 1.0)
    for ζ in ζs
        nb = neighborhood(ch, SVector(ζ))
        t  = CompScienceMeshes.tangents(nb)[:,1]   # should be (0,1)
        J  = CompScienceMeshes.jacobian(ch, ζ)     # should be 1
        @show ζ t J
    end

    @test isapprox(CompScienceMeshes.measure(ch), 1.0; atol=1e-14, rtol=0)
    =#

    λ = 10
    k = 2π/λ

    ops = [
        Helmholtz2D.singlelayer(; wavenumber=k)
        Helmholtz2D.doublelayer(; wavenumber=k)
        Helmholtz2D.doublelayer_transposed(; wavenumber=k)
    ]
    𝒮 = 

    refstrat = BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30)

    for op in ops
        Sref = assemble(op, Xt, Xs; quadstrat=refstrat)
        ref = Sref[1, 1]

        n = 25

        #Sgl = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumQStrat(n, n+1))
        Sss = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumSauterQstrat(3,3,0,4,n,n))

        #gl = Sgl[1, 1]
        ss = Sss[1, 1]

        #glrel = norm(gl - ref) /  norm(ref)
        ssrel = norm(ss - ref) /  norm(ref)
        @show typeof(op), ssrel

        #@test ssrel < 1e-14
    end

    for op in ops
        ref = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30))
        for q in (10,15,20)
            Z = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumSauterQstrat(3,3,0,4,q,q))
            δ = norm(Z - ref)/norm(ref)
            @show δ
            #@test δ < 1e-3              # or compare δ(q) strictly decreases
        end
    end

    Xt = lagrangec0_gmsh_exact(Γt; order=2)#lagrangecxd0(Γt)
    Xs = lagrangec0_gmsh_exact(Γs; order=2)#lagrangecxd0(Γs)
    
    λ = 10
    k = 2π/λ

    ops = [
        Helmholtz2D.singlelayer(; wavenumber=k)
        Helmholtz2D.doublelayer(; wavenumber=k)
        Helmholtz2D.doublelayer_transposed(; wavenumber=k)
    ]
    𝒮 = 

    refstrat = BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30)

    for op in ops
        Sref = assemble(op, Xt, Xs; quadstrat=refstrat)
        ref = Sref[1, 1]

        n = 25

        #Sgl = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumQStrat(n, n+1))
        Sss = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumSauterQstrat(3,3,0,4,n,n))

        #gl = Sgl[1, 1]
        ss = Sss[1, 1]

        #glrel = norm(gl - ref) /  norm(ref)
        ssrel = norm(ss - ref) /  norm(ref)
        @show typeof(op), ssrel

        #@test ssrel < 1e-14
    end

    for op in ops
        ref = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30))
        for q in (10,15,20)
            Z = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumSauterQstrat(3,3,0,4,q,q))
            δ = norm(Z - ref)/norm(ref)
            @show δ
            #@test δ < 1e-3              # or compare δ(q) strictly decreases
        end
    end
#end


Γto = load_gmsh_mesh(
        joinpath(@__DIR__, "assets", "circle2d_quadratic.msh");
        udim   = 2,
        element = :line,
        order   = 2,
        physical = "Boundary",
    )
X = lagrangec0_gmsh_exact(Γto; order=2)


using Gmsh
using StaticArrays
using CompScienceMeshes  # your CurvilinearMesh ctor (verts, faces, order)

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


function lagrangec0_gmsh_exact(mesh)
    # Infer cell node count from faces
    Ns = unique(length(f) for f in mesh.faces)
    @assert length(Ns) == 1 "Non-uniform local node counts across cells: $(collect(Ns))"
    NF = only(Ns)
    order = NF - 1

    # Types
    T = coordtype(mesh)
    P = vertextype(mesh)
    S = BEAST.Shape{T}

    # Map global id -> attached local shapes
    dofs = Dict{Int, Vector{S}}()
    for (c, cell) in pairs(mesh.faces)
        @inbounds for r in 1:NF
            g = cell[r]
            push!(get!(dofs, g, Vector{S}()), S(c, r, T(1)))
        end
    end

    gids = sort!(collect(keys(dofs)))
    fns  = Vector{Vector{S}}(undef, length(gids))
    pos  = Vector{P}(undef, length(gids))
    @inbounds for (i, g) in enumerate(gids)
        fns[i] = dofs[g]
        pos[i] = mesh.vertices[g]
    end

    return BEAST.LagrangeBasis{order, 0, NF}(mesh, fns, pos)
end

N = 64
Γp2 = circle_curvilinear(1.0, 2; h =2π*1.0/N)
Γp3 = circle_curvilinear(1.0, 3; h =2π*1.0/N) #(1.0, 3; h=0.1)

using BEAST
using Test
using LinearAlgebra

using BEAST, CompScienceMeshes, StaticArrays

using Gmsh

X2 = lagrangec0_gmsh_exact(Γp2)

X3 = lagrangec0_gmsh_exact(Γp3)

    Xt = X2
    Xs = Xt
    λ = 10
    k = 2π/λ

    ops = [
        Helmholtz2D.singlelayer(; wavenumber=k)
    ]

    refstrat = BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30)

    for op in ops
        Sref = assemble(op, Xt, Xs; quadstrat=refstrat)
        ref = Sref[1, 1]

        n = 25

        #Sgl = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumQStrat(n, n+1))
        Sss = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumSauterQstrat(3,3,0,4,n,n))

        #gl = Sgl[1, 1]
        ss = Sss[1, 1]

        #glrel = norm(gl - ref) /  norm(ref)
        ssrel = norm(ss - ref) /  norm(ref)
        @show typeof(op), ssrel

        @test ssrel < 1e-14
    end

    #T = Float64
    Xt = lagrangec0_gmsh_exact(Γt; order=2) #lagrangec0d1(Γt)
    Xs = lagrangec0_gmsh_exact(Γs; order=2) #lagrangec0d1(Γs)

    ops = [
        Helmholtz2D.singlelayer(; wavenumber=k)
        Helmholtz2D.hypersingular(; wavenumber=k)
    ]
    𝒮 = 

    refstrat = BEAST.DoubleNumSauterQstrat(3,3,0,4,30,30)

    for op in ops
        Sref = assemble(op, Xt, Xs; quadstrat=refstrat)
        ref = Sref[1, 1]

        n = 25

        #Sgl = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumQStrat(n, n+1))
        Sss = assemble(op, Xt, Xs; quadstrat=BEAST.DoubleNumSauterQstrat(3,3,0,4,n,n))

        #gl = Sgl[1, 1]
        ss = Sss[1, 1]

        #glrel = norm(gl - ref) /  norm(ref)
        ssrel = norm(ss - ref) /  norm(ref)
        @show typeof(op), ssrel

        tol = op isa BEAST.HH2DHyperSingularFDBIO ? 2e-6 : 1e-11

        @test ssrel < tol
    end