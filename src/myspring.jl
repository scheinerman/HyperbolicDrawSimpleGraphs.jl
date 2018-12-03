function h_layout_spring_adj(adj_matrix::Array{T,2}; C=2.0, MAXITER=100, INITTEMP=2.0) where T

    size(adj_matrix, 1) != size(adj_matrix, 2) && error("Adj. matrix must be square.")
    N = size(adj_matrix, 1)

    # Initial layout is random in the unit circle
    locs_r = rand(N)
    angs = [2*pi*t/N for t=0:N-1]

    # The optimal distance bewteen vertices
    K = C * sqrt(4.0 / N)

    # Store forces and apply at end of iteration all at once
    force = zeros(N)

    # Iterate MAXITER times
    @inbounds for iter = 1:MAXITER
        # Calculate forces
        for i = 1:N
            force_vec = 0.0
            for j = 1:N
                i == j && continue
                P = HPoint(locs_r[i],angs[i])
                Q = HPoint(locs_r[j],angs[j])
                d = dist(P,Q)
                if adj_matrix[i,j] != zero(eltype(adj_matrix)) || adj_matrix[j,i] != zero(eltype(adj_matrix))
                    # F = d^2 / K - K^2 / d
                    F_d = d / K - K^2 / d^2
                else
                    # Just repulsive
                    # F = -K^2 / d^
                    F_d = -K^2 / d^2
                end
                # d  /          sin θ = d_y/d = fy/F
                # F /| dy fy    -> fy = F*d_y/d
                #  / |          cos θ = d_x/d = fx/F
                # /---          -> fx = F*d_x/d
                # dx fx
                force_vec += F_d*d

            end
            force[i] = force_vec
        end
        # Cool down
        TEMP = INITTEMP / iter
        # Now apply them, but limit to temperature
        for i = 1:N
            force_mag  = force[i]
            scale      = min(force_mag, TEMP)/force_mag
            locs_r[i] += force_r[i] * scale
        end
    end

    # Scale to unit circle
    min_r, max_r = minimum(locs_r), maximum(locs_r)
    function scaler(z, a, b)
        2.0*((z - a)/(b - a)) - 1.0
    end
    locs_r = map(z -> scaler(z, min_r, max_r), locs_r)

    return locs_r
end
