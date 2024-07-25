

mutable struct PlaneWaveDirichlet{T,P} <: Functional
    wavenumber::T
    direction::P
end

scalartype(x::PlaneWaveDirichlet{T}) where {T} = complex(T)

mutable struct PlaneWaveNeumann{T,P} <: Functional
    wavenumber::T
    direction::P
end

scalartype(x::PlaneWaveNeumann{T}) where {T} = complex(T)

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

function (field::PlaneWaveDirichlet)(mp)

    wavenumber = field.wavenumber
    direction  = field.direction

    cart = cartesian(mp)
    exp(-im * wavenumber * dot(direction, cart))
end


function(field::PlaneWaveNeumann)(mp)

    wavenumber = field.wavenumber
    direction  = field.direction

    cart = cartesian(mp)
    norm = normal(mp)

    wave = exp(-im * wavenumber * dot(direction, cart))
    grad = -im * wavenumber * direction * wave

    d = norm[1] * grad[1]
    for i in 2:length(norm)  d += norm[i] * grad[i]  end
    return d
end

function integrand(pw::PlaneWaveDirichlet, sv, fx)
    tx = sv[1]
    return dot(tx, fx)
end

function integrand(pw::PlaneWaveNeumann, sv, fx)
    tx = sv[1]
    d = tx[1] * fx[1]
    for i in 2:length(tx) d += tx[i]*fx[i] end
    return d
end