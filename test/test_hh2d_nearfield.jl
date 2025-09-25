using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test

function meshsquare(side::T, delta::T, udim=2) where T<:Real
    PT = SVector{udim,T}
    CT = SVector{2,Int}

    # Total perimeter
    perim = 4 * side
    num_segments = ceil(Int, perim / delta)
    num_segments = max(num_segments, 4)  # At least one segment per side

    # Segments per side (evenly distribute)
    seg_per_side = ceil(Int, num_segments / 4)
    actual_delta = side / seg_per_side

    # Build square points (start at bottom left, CCW)
    points = SVector{udim,T}[]
    
    for i in 0:seg_per_side-1
        push!(points, PT(side * i / seg_per_side, 0))               # Bottom
    end
    for i in 0:seg_per_side-1
        push!(points, PT(side, side * i / seg_per_side))            # Right
    end
    for i in 0:seg_per_side-1
        push!(points, PT(side * (seg_per_side - i) / seg_per_side, side))  # Top
    end
    for i in 0:seg_per_side-1
        push!(points, PT(0, side * (seg_per_side - i) / seg_per_side))     # Left
    end

    # Build edges
    N = length(points)
    faces = [CT(i, i % N + 1) for i in 1:N]

    return Mesh(points, faces)
end



##
#hs = [0.8, 0.4, 0.2, 0.1, 0.05, 0.025]

hs = [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.01562]

orders = [1, 2, 3, 4, 5]

accs = Dict()

for order in orders
    accs[order] = Float64[]

    for h in hs
        r = 10.0
        λ = 20 * r
        k = 2 * π / λ
        # sphere = meshsphere(r, 0.2 * r)

        square = translate(meshsquare(2*r, h*r, 2), r*SVector(-0.5, -0.5))

        #k = 0.031415926535897934
        println(order)
        #if order == 0
            #X0 = BEAST.lagrangecxd0(square) # For comparison to zeroth order (discont.)
        #    X0 = BEAST.lagrangec0d1(square, dirichlet=true) #For comparison to linear order (cont.)
        #else
            #X0 = BEAST.lagrangecx(square, order=order)
            X0 = BEAST.lagrangec0(square, order=order, dirichlet=true)
        #end
        # X0 = BEAST.lagrangec0(square, order=order, dirichlet=true)
        #X0 = BEAST.lagrangec0d1(square, dirichlet=true) For comparison to linear order (cont.)
        @show numfunctions(X0)
        #X1 = lagrangec0d1(square)

        S = Helmholtz2D.singlelayer(; gamma=im * k)

        q = 100.0
        ϵ = 1.0

        # Interior problem
        # Formulations from Sauter and Schwab, Boundary Element Methods(2011), Chapter 3.4.1.1

        pos1 = SVector(r * 30.0, 0.0)  # positioning of point charges
        pos2 = SVector(-r * 30.0, 0.0)

        charge1 = Helmholtz2D.monopole(position=pos1, amplitude=q/(4*π*ϵ), wavenumber=k)
        charge2 = Helmholtz2D.monopole(position=pos2, amplitude=-q/(4*π*ϵ), wavenumber=k)

        # Potential of point charges

        Φ_inc(x) = charge1(x) + charge2(x)

        gD0 = assemble(DirichletTrace(charge1), X0) + assemble(DirichletTrace(charge2), X0)
    
        # Interior Dirichlet problem - compare Sauter & Schwab eqs. 3.81
        @time M_IDPSL = assemble(S, X0, X0) # Single layer (SL)
    
        ρ_IDPSL = M_IDPSL \ (-gD0)

        pts = meshcircle(0.4 * r, 0.4 * 0.6 * r).vertices # sphere inside on which the potential and field are evaluated

        pot_IDPSL = potential(BEAST.HH2DSingleLayerNear(im * k), pts, ρ_IDPSL, X0; type=ComplexF64)

        # Total field inside should be zero
        err_IDPSL_pot = norm(pot_IDPSL + Φ_inc.(pts)) / norm(Φ_inc.(pts))

        push!(accs[order], err_IDPSL_pot)
    end
end
##

using Plots
plotlyjs()
plt = Plots.plot(
    xscale=:log10,
    yscale=:log10,
    xlabel="h in m",
    ylabel="Relative error",
    legend=:bottomright,
    title="Manufactured solution for square (Helmholtz2D)",  linewidth=3)

    #Plots.plot!(plt, hs, accs[0], label="Order 0", markershape=:auto)
    Plots.plot!(plt, hs, accs[1], label="Order 1", markershape=:auto)
    Plots.plot!(plt, hs, accs[2], label="Order 2", markershape=:auto)
    Plots.plot!(plt, hs, accs[3], label="Order 3", markershape=:auto)
    Plots.plot!(plt, hs, accs[4], label="Order 4", markershape=:auto)
    Plots.plot!(plt, hs, accs[5], label="Order 5", markershape=:auto)

    savefig("clagrangecx_square_relative_error.pdf")


for i in orders
    Plots.plot!(plt, hs, accs[i], label="Order $i", markershape=:auto)
end
display(plt)
##

hs = [0.8, 0.4, 0.2, 0.1]

orders = [0, 1]

acc = [Float64[], Float64[]]
for order in orders

    for h in hs
        r = 10.0
        λ = 20 * r
        k = 2 * π / λ
        
        sphere = CompScienceMeshes.meshsphere(r, h * r)

        #k = 0.031415926535897934

        X0 = BEAST.lagrangecx(sphere, order=order)
        @show numfunctions(X0)
        #X1 = lagrangec0d1(square)

        S = Helmholtz3D.singlelayer(;)

        q = 100.0
        ϵ = 1.0

        # Interior problem
        # Formulations from Sauter and Schwab, Boundary Element Methods(2011), Chapter 3.4.1.1

        pos1 = SVector(r * 1.5, 0.0, 0.0)  # positioning of point charges
        pos2 = SVector(-r * 1.5, 0.0, 0.0)

        charge1 = Helmholtz3D.monopole(position=pos1, amplitude=q/(4*π*ϵ))
        charge2 = Helmholtz3D.monopole(position=pos2, amplitude=-q/(4*π*ϵ))

        # Potential of point charges

        Φ_inc(x) = charge1(x) + charge2(x)

        gD0 = assemble(DirichletTrace(charge1), X0) + assemble(DirichletTrace(charge2), X0)
    
        # Interior Dirichlet problem - compare Sauter & Schwab eqs. 3.81
        M_IDPSL = assemble(S, X0, X0) # Single layer (SL)
    
        ρ_IDPSL = M_IDPSL \ (-gD0)

        pts = meshsphere(0.8 * r, 0.8 * 0.6 * r).vertices # sphere inside on which the potential and field are evaluated

        pot_IDPSL = potential(BEAST.HH3DSingleLayerNear(0.0), pts, ρ_IDPSL, X0; type=ComplexF64)

        # Total field inside should be zero
        err_IDPSL_pot = norm(pot_IDPSL + Φ_inc.(pts)) / norm(Φ_inc.(pts))

        push!(acc[order+1], err_IDPSL_pot)
    end
end

##

hs = [2.0, 1.0, 0.5, 0.25]

orders = [0, 1, 2, 3]

accs = Dict()

for order in orders
    accs[order] = Float64[]
    for h in hs
        r = 1.0
        λ = 20 * r
        k = 2 * π / λ
        
        cube =  translate(CompScienceMeshes.meshcuboid(2.0, 2.0, 2.0, h), -SVector(1.0,1.0,1.0))

        #k = 0.031415926535897934

        X0 = BEAST.lagrangecx(cube, order=order)
        @show numfunctions(X0)
        #X1 = lagrangec0d1(square)

        S = Helmholtz3D.singlelayer(;)

        q = 100.0
        ϵ = 1.0

        # Interior problem
        # Formulations from Sauter and Schwab, Boundary Element Methods(2011), Chapter 3.4.1.1

        pos1 = SVector(r * 3.0, 0.0, 0.0)  # positioning of point charges
        pos2 = SVector(-r * 3.0, 0.0, 0.0)

        charge1 = Helmholtz3D.monopole(position=pos1, amplitude=q/(4*π*ϵ))
        charge2 = Helmholtz3D.monopole(position=pos2, amplitude=-q/(4*π*ϵ))

        # Potential of point charges

        Φ_inc(x) = charge1(x) + charge2(x)

        gD0 = assemble(DirichletTrace(charge1), X0) + assemble(DirichletTrace(charge2), X0)
    
        # Interior Dirichlet problem - compare Sauter & Schwab eqs. 3.81
        M_IDPSL = assemble(S, X0, X0) # Single layer (SL)
    
        ρ_IDPSL = M_IDPSL \ (-gD0)

        pts = meshsphere(0.8 * r, 0.8 * 0.6 * r).vertices # sphere inside on which the potential and field are evaluated

        pot_IDPSL = potential(BEAST.HH3DSingleLayerNear(0.0), pts, ρ_IDPSL, X0; type=ComplexF64)

        # Total field inside should be zero
        err_IDPSL_pot = norm(pot_IDPSL + Φ_inc.(pts)) / norm(Φ_inc.(pts))

        push!(accs[order], err_IDPSL_pot)
    end
end

##

using Plots
plotlyjs()
plt = Plots.plot(
    xscale=:log10,
    yscale=:log10,
    xlabel="h in m",
    ylabel="Relative error",
    legend=:bottomright,
    title="Manufactured solution for cube (Helmholtz3D)",  linewidth=3)
Plots.plot!(plt, hs, accs[0], label="Order 0", markershape=:circle)
Plots.plot!(plt, hs, accs[1], label="Order 1", markershape=:square)
Plots.plot!(plt, hs, accs[2], label="Order 2", markershape=:diamond)
Plots.plot!(plt, hs, accs[3][1] .* hs.^5, label="h", markershape=:dtriangle)
Plots.plot!(plt, hs, accs[3], label="Order 3", markershape=:utriangle)