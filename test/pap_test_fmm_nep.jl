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

function hassemble(op, X, Y)
    testtree = KMeansTree(X.pos, 2; minvalues=100)
    trialtree = KMeansTree(Y.pos, 2; minvalues=100)
    tree = BlockTree(testtree, trialtree)
    nearquadstrat = BEAST.DoubleNumSauterQstrat(3, 4, 0, 4, 5, 6)
    farquadstrat = BEAST.DoubleNumQStrat(3, 4)
    #scheduler=SerialScheduler() # single threaded
    scheduler=DynamicScheduler() # multi threaded
    @time hmat = AdaptiveCrossApproximation.HMatrix(
        op, X, Y, tree;
        nearquadstrat=nearquadstrat, farquadstrat=farquadstrat, scheduler=scheduler
    );
    return hmat
end

##

using Krylov

a       = 1.0                      # radius
p_order = 2                        # polynomial order of boundary elements
h       = 2π * a / 5000

k = besselj_zero(0, 1) / a   # ≈ 2.4048255577 dla a=1

circle = CompScienceMeshes.meshcircle(a, h)

X = lagrangecx(circle, order=p_order)

@show typeof(X)
@show size(getproperty(X, :pos))

@assert hasproperty(X, :pos)
#@assert size(X.pos, 1) == 2   # dla 2D

length(X)
𝒟 = Helmholtz2D.doublelayer(; wavenumber = k)

##

K = hassemble(𝒟, X, X)

#
I = assemble(BEAST.Identity(), X, X)
T = -0.5 * I + K
μ_true = randn(ComplexF64, size(T, 2))
T*μ_true
g = Krylov.gmres(T, μ_true;rtol=1e-4, itmax=200,verbose=1, history=true)


##
D  = assemble(𝒟, X, X)




function assemble_T_circle_Dirichlet(k::Number, X)
    𝒟 = Helmholtz2D.doublelayer(; wavenumber = k)
    D  = assemble(𝒟, X, X)
    I  = assemble(BEAST.Identity(), X, X)
    return Matrix(-0.5 .* I .+ D)
end


a       = 1.0                      # radius
p_order = 2                        # polynomial order of boundary elements
h       = 2π * a / 64  

k = besselj_zero(0, 1) / a   # ≈ 2.4048255577 dla a=1

circle = CompScienceMeshes.meshcircle(a, h)
##
X = lagrangecx(circle, order=p_order)
𝒟 = Helmholtz2D.doublelayer(; wavenumber = k)
D  = assemble(𝒟, X, X)
##




I  = assemble(BEAST.Identity(), X, X)

##
T = assemble_T_circle_Dirichlet(k, X)

D = size(T,1)
μ_true = randn(ComplexF64, D)
g = T * μ_true
μ = T \ g

@show norm(μ-μ_true)/norm(μ_true)


using Krylov, LinearMaps, IterativeSolvers, LinearSolve, IncompleteLU, AlgebraicMultigrid, Preconditioners
#using FMM2D