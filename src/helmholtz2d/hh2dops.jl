import SpecialFunctions: hankelh2

abstract type HelmholtzOperator2D <: IntegralOperator end
scalartype(::HelmholtzOperator2D) = ComplexF64

struct SingleLayer{T} <: HelmholtzOperator2D
    wavenumber::T
end

struct HyperSingular{T} <: HelmholtzOperator2D
    wavenumber::T
end

struct DoubleLayer{T} <: HelmholtzOperator2D
    wavenumber::T
end

struct DoubleLayerTransposed{T} <: HelmholtzOperator2D
    wavenumber::T
end



mutable struct KernelValsHelmholtz2D
    wavenumber
    vect
    dist
    green
    gradgreen
    txty
end


function kernelvals(biop::HelmholtzOperator2D, tgeo, bgeo)

    k = biop.wavenumber
    r = tgeo.cart - bgeo.cart
    R = norm(r)

    kr = k * R
    hankels = hankelh2.([0 1], kr)
    green = -im / 4 * hankels[1]
    gradgreen = k * im / 4 * hankels[2] * r / R

    txty = dot(normal(tgeo), normal(bgeo))

    KernelValsHelmholtz2D(k, r, R, green, gradgreen, txty)
end


shapevals(op::HelmholtzOperator2D, ϕ, ts) = shapevals(ValDiff(), ϕ, ts)


function integrand(biop::SingleLayer, kerneldata, tvals,
    tgeo, bvals, bgeo)

    gx = tvals[1]
    fy = bvals[1]

    gx * kerneldata.green * fy
end

function integrand(biop::HyperSingular, kernel, tvals, tgeo,
    bvals, bgeo)

    gx = tvals[1]
    fy = bvals[1]

    dgx = tvals[2]
    dfy = bvals[2]

    k    = kernel.wavenumber
    G    = kernel.green
    txty = kernel.txty

    (dgx * dfy - k*k * txty * gx * fy) * G
end

function integrand(biop::DoubleLayer, kernel, fp, mp, fq, mq)
    nq = normal(mq)
    fp[1] * dot(nq, -kernel.gradgreen) * fq[1]
end

function integrand(biop::DoubleLayerTransposed, kernel, fp, mp, fq, mq)
    np = normal(mp)
    fp[1] * dot(np, kernel.gradgreen) * fq[1]
end

function cellcellinteractions!(biop::HelmholtzOperator2D, tshs, bshs, tcell, bcell, z)

    regularcellcellinteractions!(biop, tshs, bshs, tcell, bcell, z)

end

function testfunc1()
    print("test function!")
end

# defaultquadstrat(op::HelmholtzOperator2D, tfs, bfs) = DoubleNumWiltonSauterQStrat(4,3,4,3,4,4,4,4)
defaultquadstrat(op::HelmholtzOperator2D, tfs, bfs) = DoubleNumQStrat(4,3)


function quaddata(op::HelmholtzOperator2D, g::LagrangeRefSpace, f::LagrangeRefSpace, tels, bels,
        qs::DoubleNumWiltonSauterQStrat)

  tqd = quadpoints(g, tels, (qs.outer_rule_far,))
  bqd = quadpoints(f, bels, (qs.inner_rule_far,))

  return (tpoints=tqd, bpoints=bqd)
end

function quadrule(op::HelmholtzOperator2D, g::LagrangeRefSpace, f::LagrangeRefSpace, i, τ, j, σ, qd,
        qs::DoubleNumWiltonSauterQStrat)

    DoubleQuadRule(
        qd.tpoints[1,i],
        qd.bpoints[1,j]
    )

end


