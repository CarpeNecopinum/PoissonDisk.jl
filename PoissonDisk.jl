cd(@__DIR__); pwd()
import Pkg
Pkg.activate(".")

using GeometryTypes
using Distributions
using Base.Cartesian: @nloops, @nref, @ncall
using LinearAlgebra: norm


function sample_sphere_annulus(::Type{Vec{N,T}}) where {N,T}
    e = log(1-T(0.5)^N) / log(0.5)

    dir = normalize(randn(Vec{N,T}))
    radius = T((1 - rand(T) * T(0.5)^e) ^ (1/N))
    T.(dir .* radius)
end

function disksampling(dom_min::Vec2{T}, dom_max::Vec2{T}, r::Real, k::Integer = 30) where {T <: Real}
    N = 2
    cell_size = r / sqrt(N)
    dom_size = dom_max - dom_min
    bg_grid = zeros(Int, Tuple(ceil.(Int, (dom_max - dom_min) ./ cell_size)))
    points = [rand(Vec{N,T}) .* dom_size + dom_min]
    active = [1]

    cell_for(x) = ceil.(Int, (x - dom_min) ./ cell_size)

    function is_okay(x)
        all(dom_min .< x .< dom_max) || (return false)
        center = cell_for(x)
        @nloops 2 i i->-1:1 begin
            coord = (@ncall 2 tuple i) .+ center
            (checkbounds(Bool, bg_grid, coord...)) || continue
            j = bg_grid[coord...]
            (j == 0) && continue
            (norm(x - points[j]) <= r) && return false
        end
        return true
    end

    while !isempty(active)
        idx = rand(eachindex(active))
        i = active[idx]
        deleteat!(active, idx)
        x = points[i]

        for _ in 1:k
            y = sample_sphere_annulus(Vec{N,T}) .* 2r + x
            if is_okay(y)
                push!(points, y)
                push!(active, length(points))
                bg_grid[cell_for(y)...] = length(points)
                println("Added point $(y)")
                break
            end
        end
    end
    points
end


using Plots


second(x) = x[2]

xs = disksampling(Vec2f0(0), Vec2f0(100), 5f0, 20000)
Plots.scatter(first.(xs), second.(xs))
