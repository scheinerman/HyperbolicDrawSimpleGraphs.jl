
module HyperbolicDrawSimpleGraphs

using SimpleGraphs, HyperbolicPlane, SimpleDrawing

import Base: show

export HyperbolicGraphEmbedding, hdraw, hembed, convertEmbed

mutable struct HyperbolicGraphEmbedding
    G::SimpleGraph
    locs::Dict{Any,HPoint}
end

function show(io::IO, X::HyperbolicGraphEmbedding)
    print(io, "Hyperbolic embedding of $(X.G)")
end

# convert embedded SimpleGraph to Hyperbolic embedding
function convertEmbed(GG::SimpleGraph)
    n = NV(GG)
    VV = vlist(GG)
    loc = getxy(GG)
    d = Dict{Any,HPoint}()
    for i = 1:n
        l = loc[i]
        v = VV[i]
        r = sqrt(l[1]*l[1] + l[2]*l[2])
        angs = angle(l)
        P = HPoint(r,angs)
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

"""
`hembed(G,method)` embeds `G` in the hyperbolic plane by placing vertices
around a circle.

Options for `method`:
+ `:circular` places vertices around a circle
+ `:random` places vertices at random
"""
function hembed(G::SimpleGraph, method::Symbol = :circular)

    if method == :circular
        X = h_circular(G)
    elseif method == :random
        X = h_random(G)
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
    draw(HPlane())
    finish()
end

end  # end of module
