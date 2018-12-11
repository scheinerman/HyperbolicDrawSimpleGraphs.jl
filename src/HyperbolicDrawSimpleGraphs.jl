
module HyperbolicDrawSimpleGraphs

using SimpleGraphs, HyperbolicPlane, SimpleDrawing

import Base: show
import HyperbolicPlane: polar

export HyperbolicGraphEmbedding, hdraw, hembed

mutable struct HyperbolicGraphEmbedding
    G::SimpleGraph
    locs::Dict{Any,HPoint}
end

function show(io::IO, X::HyperbolicGraphEmbedding)
    print(io, "Hyperbolic embedding of $(X.G)")
end

polar(x::Real,y::Real) = (sqrt(x*x+y*y), atan(x,y) )

# convert SimpleGraph to Hyperbolic embedding
function h_convert(GG::SimpleGraph)::HyperbolicGraphEmbedding
    loc = getxy(GG)
    d = Dict{Any,HPoint}()

    for v in GG.V
        xy = loc[v]
        (r,theta) = polar(xy[1],xy[2])
        P = HPoint(r,theta)
        d[v] = P
    end
    return HyperbolicGraphEmbedding(GG,d)
end

function h_circular(GG::SimpleGraph)::HyperbolicGraphEmbedding
    n = NV(GG)
    VV = vlist(GG)
    r = sqrt(n)
    angs = [2*pi*t/n for t=0:n-1]
    d = Dict{Any,HPoint}()
    for j=1:n
        v = VV[j]
        P = HPoint(r,angs[j])
        d[v] = P
    end
    return HyperbolicGraphEmbedding(GG,d)
end

function h_random(G::SimpleGraph)::HyperbolicGraphEmbedding
    d = Dict{Any,HPoint}()
    for v in G.V
        d[v] = RandomHPoint()
    end
    return HyperbolicGraphEmbedding(G,d)
end

function private_adj(G::SimpleGraph)
    n = NV(G)
    A = zeros(Int,n,n)
    vv = vlist(G)
    for i=1:n-1
        u = vv[i]
        for j=(i+1):n
            w = vv[j]
            if has(G,u,w)
                A[i,j] = 1
                A[j,i] = 1
            end
        end
    end
    return A,vv
end

include("myspring.jl")

"""
`spring!(X)` gives the graph held in `X` with a spring embedding
(based on code in the `GraphLayout` module). If runs a default number of
iterations (100) of that algorithm. To change the number of
iterations, use `spring!(X,nits)`.
"""
function h_spring(G::SimpleGraph, nits::Int=100)::HyperbolicGraphEmbedding
    GG = deepcopy(G)
    embed(GG,:spring)
    X = h_convert(GG)
    locs_r = Array{Float64,1}
    angs = Array{Float64,1}
    i = 1
    for v in G.V
        xy = X.locs[v]
        (r,theta) = polar(xy)
        locs_r[i] = r
        angs[i] = theta
        i = i + 1
    end
    n = NV(G)
    A,vv = private_adj(G)

    locs_r, angs = h_layout_spring_adj(locs_r,angs,A,MAXITER=nits)

    d = Dict{Any,HPoint}()
    for i = 1:n
        v = vv[i]
        P = HPoint(locs_r[i],angs[i])
        d[v] = P
    end

    return HyperbolicGraphEmbedding(G,d)
end

"""
`hembed(G,method)` embeds `G` in the hyperbolic plane by placing vertices
around a circle.

Options for `method`:
+ `:circular` places vertices around a circle
+ `:random` places vertices at random
+ `:convert` converts the Euclidean embedding of `G` into a hyperbolic embedding
"""
function hembed(G::SimpleGraph, method::Symbol = :circular)

    if method == :circular
        X = h_circular(G)
    elseif method == :random
        X = h_random(G)
    elseif method == :convert
        X = h_convert(G)
    elseif method == :spring
        X = h_spring(G)
    else
        error("Unknown method $method")
    end
    cache_save(G,:HyperbolicGraphEmbedding,X)
end

"""
`hdraw(G)` draws the graph `G` in its current hyperoblic embedding
(or creates a default, circular embedding if need be).
"""
function hdraw(G::SimpleGraph)
    if !cache_check(G,:HyperbolicGraphEmbedding)
        hembed(G)
    end
    X = cache_recall(G,:HyperbolicGraphEmbedding)
    newdraw()
    for e in G.E
        u,v = e
        A = X.locs[u]
        B = X.locs[v]
        draw(A+B)
    end
    for v in G.V
        A = X.locs[v]
        A = HPoint(A)
        set_radius(A,2)
        draw(A)
    end
    # draw(HPlane()) # let the user decide to show the line at infinity
    finish()
end

end  # end of module
