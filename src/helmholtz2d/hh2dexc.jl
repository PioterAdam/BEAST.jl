
struct HH2DPlaneWave{P,K,T} <: Functional
    direction::P
    gamma::K
    amplitude::T
end

function (f::HH2DPlaneWave)(r)
    d = f.direction
    γ = f.gamma
    a = f.amplitude
    return a * exp(-γ*dot(d,r))
end

function (f::HH2DPlaneWave)(r::CompScienceMeshes.MeshPointNM)
    return f(cartesian(r))
end

scalartype(f::HH2DPlaneWave{P,K,T}) where {P,K,T} = promote_type(eltype(P), K, T)

struct gradHH2DPlaneWave{P,K,T} <: Functional
    direction::P
    gamma::K
    amplitude::T
end

function (f::gradHH2DPlaneWave)(r)
    d = f.direction
    γ = f.gamma
    a = f.amplitude

    return -γ * d * exp(-γ * dot(d, r))
end

function (f::gradHH2DPlaneWave)(mp::CompScienceMeshes.MeshPointNM)
    r = cartesian(mp)
    return dot(normal(mp), f(r))
end

scalartype(f::gradHH2DPlaneWave{P,K,T}) where {P,K,T} = promote_type(eltype(P), K, T)

mutable struct ScalarTrace{T,F} <: Functional
    field::F
end

ScalarTrace(f::F) where {F} = ScalarTrace{scalartype(f), F}(f)
ScalarTrace{T}(f::F) where {T,F} = ScalarTrace{T,F}(f)

strace(f, mesh::Mesh) = ScalarTrace(f)

(s::ScalarTrace)(x) = s.field(cartesian(x))
integrand(s::ScalarTrace, tx, fx) = dot(tx.value, fx)
scalartype(s::ScalarTrace{T}) where {T} = T

shapevals(f::Functional, ϕ, ts) = shapevals(ValOnly(), ϕ, ts)

function (f::NormalDerivative{T,F})(manipoint) where {T,F<:HH2DPlaneWave}
    d = f.field.direction
    γ = f.field.gamma
    a = f.field.amplitude
    n = normal(manipoint)
    r = cartesian(manipoint)
    return -γ*a * dot(d,n) * exp(-γ*dot(d,r))
end

*(a::Number, m::HH2DPlaneWave) = HH2DPlaneWave(m.direction, m.gamma, a * m.amplitude)
*(a::Number, m::gradHH2DPlaneWave) = gradHH2DPlaneWave(m.direction, m.gamma, a * m.amplitude)

dot(::NormalVector, m::gradHH2DPlaneWave) = NormalDerivative(HH2DPlaneWave(m.direction, m.gamma, m.amplitude))