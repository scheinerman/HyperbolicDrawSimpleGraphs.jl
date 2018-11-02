
module HyperbolicDrawSimpleGraphs

using SimpleGraphs, HyperbolicPlane, SimpleDrawing


export HyperbolicGraphEmbedding, hdraw, hembed

mutable struct HyperbolicGraphEmbedding
    G::SimpleGraph
    locs::Dict{Any,HPoint}
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

"""
`hembed(G)` embeds `G` in the hyperbolic plane by placing vertices
around a circle.
"""
function hembed(G::SimpleGraph)
    X = h_circular(G)
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
    for e in G.E
        u,v = e
        A = X.locs[u]
        B = X.locs[v]
        draw(A+B)
    end
    for v in G.V
        A = X.locs[v]
        draw(A)
    end
    finish()
end

end
