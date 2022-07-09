"""
    ABCDHParams

A structure holding parameters for ABCD graph generator. Fields:
* is_simple::Bool:             if hypergraph is simple
* w::Vector{Int}:             list of vertex degrees
* y::Vector{Int}:             community degree
* z::Vector{Int}:             network degree
* s::Vector{Int}:             list of cluster sizes
* ξ::Union{Float64, Nothing}: background graph fraction
* q::Vector{Float64}:         distribution of hyperedge sizes
* wcd::Matrix{Float64}: desired composition of hyperedges
* maxiter: maximum number of iterations trying to resolve collisions per collision
"""
struct ABCDHParams
    is_simple::Bool
    w::Vector{Int}
    y::Vector{Int}
    z::Vector{Int}
    s::Vector{Int}
    ξ::Float64
    q::Vector{Float64}
    wcd::Matrix{Float64}
    maxiter::Int

    function ABCDHParams(w, s, ξ, q, wcd, is_simple, maxiter)
        length(w) == sum(s) || throw(ArgumentError("inconsistent data"))
        all(>=(0), w) || throw(ArgumentError("negative degree passed"))
        all(>=(0), s) || throw(ArgumentError("negative community size passed"))
        all(>=(0), q) || throw(ArgumentError("negative hyperedge proportion passed"))
        all(>=(0), wcd) || throw(ArgumentError("negative hyperedge composition passed"))
        0 ≤ ξ ≤ 1 || throw(ArgumentError("ξ=$ξ not in [0,1] interval"))
        q = Vector{Float64}(q)
        sq = sum(q)
        if sq != 1
            @info "distribution of hyperedge proportions does not add up to 1. Fixing."
            sq < eps() && throw(ArgumentError("sum of hyperedge proportions $sq is abnormally small"))
            q ./= sq
        end
        wcd = Matrix{Float64}(wcd)
        size(wcd) == (length(q), length(q)) || throw(ArgumentError("incorrect dimension of hyperedge composition matrix"))
        for d in 1:length(q)
            for c in 1:div(d, 2)
                if wcd[c, d] != 0
                    @info "weight for c <= d/2 is positive where c=$c and d=$d. Fixing."
                    wcd[c, d] = 0.0
                end
            end
            for c in d+1:length(q)
                if wcd[c, d] != 0
                    @info "weight for c > d is positive where c=$c and d=$d. Fixing."
                    wcd[c, d] = 0.0
                end
            end
            swcd = sum(wcd[:, d])
            if swcd != 1
                @info "distribution of hyperedge composition for d=$d does not add up to 1. Fixing."
                swcd < eps() && throw(ArgumentError("sum of hyperedge composition $swcd is abnormally small for d=$d"))
                wcd[:, d] ./= swcd
            end
        end
        return new(is_simple, w, fill(-1, length(w)), fill(-1, length(w)), s, ξ, q, wcd, maxiter)
    end
end

function randround(x)
    d = floor(Int, x)
    return d + (rand() < x - d)
end

function generate_1hyperedges(params::ABCDHParams)
    n = length(params.w)
    ww = Weights(params.w)
    if params.is_simple
        m1 = min(randround(params.q[1] * sum(params.w)), n)
        he1 = sample(1:n, ww, m1, replace=false)
        @assert allunique(he1)
    else
        m1 = randround(params.q[1] * sum(params.w))
        he1 = sample(1:n, ww, m1, replace=true)
    end
    for i in he1
        @assert params.w[i] > 0
        params.w[i] -= 1
    end
    return [[i] for i in he1]
end

function populate_clusters(params::ABCDHParams)
    n = length(params.w)
    clusters = fill(-1, n)
    slots = copy(params.s)
    community_allowed = fill(true, length(slots))
    for i in sortperm(params.w, rev=true)
        loc = -1
        wts_all = Weights(slots)
        for _ in 1:10
            pos = sample(wts_all)
            choice_allowed = true
                cj = big(params.s[pos])
                for d in 2:size(params.wcd, 2), c in div(d, 2)+1:d
                    xi = sum(div(d, 2)+1:c) do f
                        params.y[i] *
                        params.q[d] *
                        params.wcd[f, d] *
                        binomial(big(d-f), big(c-f)) *
                        (cj/n)^(c-f) *
                        (1-cj/n)^(d-c) +
                        params.z[i] *
                        params.q[d] * binomial(big(d-1), big(c-1)) *
                        (cj/n)^(c-1) *
                        (1-cj/n)^(d-c)
                    end
                    if xi > binomial(cj-1, c-1) * binomial(n-cj, d-c)
                        choice_allowed = false
                        break
                    end
                end
            if choice_allowed
                loc = pos
                break
            end
        end
        if loc == -1
            community_allowed .= true
            for j in 1:length(slots)
                cj = big(params.s[j])
                for d in 2:size(params.wcd, 2), c in div(d, 2)+1:d
                    xi = sum(div(d, 2)+1:c) do f
                        params.y[i] *
                        params.q[d] *
                        params.wcd[f, d] *
                        binomial(big(d-f), big(c-f)) *
                        (cj/n)^(c-f) *
                        (1-cj/n)^(d-c) +
                        params.z[i] *
                        params.q[d] * binomial(big(d-1), big(c-1)) *
                        (cj/n)^(c-1) *
                        (1-cj/n)^(d-c)
                    end
                    if xi > binomial(cj-1, c-1) * binomial(n-cj, d-c)
                        community_allowed[j] = false
                        break
                    end
                end
                comm_idxs = findall(community_allowed)
                wts = slots[comm_idxs]
                sum(wts) == 0 && throw(ArgumentError("hypergraph is too tight. Failed to find community for node $i"))
                loc = sample(comm_idxs, Weights(wts))
            end
        end
        clusters[i] = loc
        slots[loc] -= 1
    end
    @assert sum(slots) == 0
    @assert extrema(clusters) == (1, length(params.s))
    return clusters
end

function config_model(clusters, params, he1)
    L = length(params.q)

    edges = Vector{Int}[]

    # community graphs
    community_stumps = Int[] # stumps left for filling d-c slots
    edges_with_missing_size = Tuple{Vector{Int}, Int}[] # partial edges with missing d-c slots

    for j in 1:length(params.s)
        cluster_idxs = findall(==(j), clusters)

        md = zeros(Int, L)
        pj = sum(params.y[cluster_idxs])
        for d in L:-1:2
            sumq2 = sum(params.q[2:d])
            if sumq2 > 0
                md[d] = floor(Int, params.q[d] / sumq2 * (pj - sum((d+1:L) .* md[d+1:L])) / d)
            end
        end
        sumpj = sum(d * md[d] for d in 1:L)
        if pj > sumpj
            @info "Moving $(pj - sumpj) stumps from community $j to background graph"
        end

        while pj > sumpj
            dec_idx = sample(cluster_idxs, Weights(params.y[cluster_idxs]))
            params.y[dec_idx] -= 1
            params.z[dec_idx] += 1
            pj -= 1
        end
        @assert pj == sumpj

        if pj == 0
            @assert sum(params.y[cluster_idxs]) == 0
        else
            mcd = zeros(Int, L, L)
            for d in 2:L
                for c in d:-1:div(d, 2)+1
                    sumwfd = sum(params.wcd[div(d, 2)+1:c, d])
                    if sumwfd > 0
                        mcd[c, d] = randround(params.wcd[c, d] / sumwfd * (md[d] - sum(mcd[c+1:d, d])))
                    end
                end
                @assert md[d] == sum(mcd[:, d])
            end

            @assert sum(d * mcd[c, d] for d in 2:L for c in d:-1:div(d, 2)+1) == pj

            pjc = sum(c * mcd[c, d] for d in 2:L for c in d:-1:div(d, 2)+1)
            # TODO - fixe the case when pj is 0
            yc = Float64.(params.y[cluster_idxs] * (big(pjc) / big(pj))) # this is to improve rounding behavior
            yc_base = floor.(Int, yc)
            yc_rem = yc .- yc_base
            tail_size = pjc - sum(yc_base)

            @assert 0 <= tail_size <= length(yc_rem)
            if tail_size > 0
                additional_points = sample(1:length(yc_rem), Weights(yc_rem), tail_size, replace=false)
                for point in additional_points
                    yc_base[point] += 1
                end
            end
            @assert sum(yc_base) == pjc
            @assert floor.(yc) <= yc_base <= ceil.(yc)

            @assert length(cluster_idxs) == length(yc_base)
            internal_stumps = Int[]
            for (yci, index) in zip(yc_base, cluster_idxs)
                append!(internal_stumps, fill(index, yci))
            end
            shuffle!(internal_stumps)
            @assert sum(mcd) <= length(internal_stumps)
            @assert pjc == length(internal_stumps)
            stump_idx = 1
            for d in 2:L
                for c in div(d, 2)+1:d
                    for _ in 1:mcd[c, d]
                        stump_end = stump_idx + c - 1
                        push!(edges_with_missing_size, (internal_stumps[stump_idx:stump_end], d-c))
                        stump_idx += c
                    end
                end
            end
            @assert stump_idx == length(internal_stumps) + 1

            start_len = length(community_stumps)
            for (yri, index) in zip(params.y[cluster_idxs] .- yc_base, cluster_idxs)
                append!(community_stumps, fill(index, yri))
            end
            end_len = length(community_stumps)
            @assert end_len - start_len == pj - pjc
        end
    end

    shuffle!(community_stumps)

    stump_idx = 1
    for (he, dc) in edges_with_missing_size
        if dc > 0
            stump_end = stump_idx + dc - 1
            append!(he, community_stumps[stump_idx:stump_end])
            stump_idx += dc
        end
        push!(edges, sort!(he))
    end

    # background graph

    # stumps transferred from community graphs to background graph
    background_stumps = community_stumps[stump_idx:end]

    if !isempty(background_stumps)
        @info "moved $(length(bacground_stumps)) non-matched stumps from community graphs to background graph"
    end

    md = zeros(Int, L)
    p = sum(params.z) + length(background_stumps)
    for d in L:-1:2
        sumq2 = sum(params.q[2:d])
        if sumq2 > 0
            md[d] = floor(Int, params.q[d] / sumq2 * (p - sum((d+1:L) .* md[d+1:L])) / d)
        end
    end

    for index in 1:length(params.z)
        append!(background_stumps, fill(index, params.z[index]))
    end
    shuffle!(background_stumps)
    @assert sum(md) <= length(background_stumps)
    stump_idx = 1
    for d in 2:L
        for _ in 1:md[d]
            stump_end = stump_idx + d - 1
            push!(edges, sort!(background_stumps[stump_idx:stump_end]))
            stump_idx += d
        end
    end

    if stump_idx <= length(background_stumps)
        left_stumps = background_stumps[stump_idx:end]
        # note that these size-1 hyperedges might be duplicate with he1 generated earlier
        if params.q[1] > 0 && (!params.is_simple || allunique(left_stumps) && length(intersect(only.(he1), left_stumps)) == 0)
            append!(edges, [[x] for x in left_stumps])
        else
            targetq = 1 + findfirst(>(0), params.q[2:end])
            to_add = targetq - length(left_stumps)
            @assert to_add > 0
            extra_stumps = sample(1:length(params.w), Weights((params.w)), to_add)
            append!(left_stumps, extra_stumps)
            @assert length(left_stumps) == targetq
            push!(edges, sort!(left_stumps))
            @info "added degree to the following nodes due to parity issues: $(extra_stumps) "
        end
    end

    append!(edges, he1)
    @assert all(issorted, edges)

    if params.is_simple
        bad_multi = 0
        bad_dup = 0
        good_edges = Set{Vector{Int}}()
        bad_edges = Vector{Int}[]
        for he in edges
            not_multi = allunique(he)
            not_dup = !(he in good_edges)
            bad_multi += !not_multi
            bad_dup += !not_dup
            if not_multi && not_dup
                push!(good_edges, he)
            else
                push!(bad_edges, he)
            end
        end

        if bad_multi > 0
            @info "fixing $bad_multi hyperedges that were multisets"
        end
        if bad_dup > 0
            @info "fixing $bad_dup hyperedges that were duplicated"
        end

        shuffle!(bad_edges)

        iters = params.maxiter * length(bad_edges)
        for _ in 1:iters
            isempty(bad_edges) && break

            bad_edge = pop!(bad_edges)
            initial_badness = (bad_edge in good_edges) + (length(bad_edge) - length(unique(bad_edge)))
            if initial_badness == 0
                push!(good_edges, bad_edge)
            else
                other_edge = rand(good_edges)
                pop!(good_edges, other_edge)
                new_split = [bad_edge; other_edge]
                shuffle!(new_split)
                new1 = sort!(new_split[1:length(bad_edge)])
                new2 = sort!(new_split[length(bad_edge)+1:end])
                final_bandess = (new1 in good_edges) + (length(new1) - length(unique(new1))) +
                                (new2 in good_edges) + (length(new2) - length(unique(new2)))
                if final_bandess < initial_badness
                    if allunique(new1) && !(new1 in good_edges)
                        push!(good_edges, new1)
                    else
                        push!(bad_edges, new1)
                    end
                    if allunique(new2) && !(new2 in good_edges)
                        push!(good_edges, new2)
                    else
                        push!(bad_edges, new2)
                    end
                else
                    pushfirst!(bad_edges, bad_edge)
                    push!(good_edges, other_edge)
                end
            end
        end

        if !isempty(bad_edges)
            @info "Failed to fix all bad edges in $(params.maxiter) rounds. " *
                  "Dropping $(length(bad_edges)) bad edges that violated " *
                  "simple hypergraph condition."
            @show bad_edges
        end

        edges = collect(good_edges)
    end

    @assert all(issorted, edges)

    return edges
end

"""
    gen_hypergraph(params::ABCDHParams)

Generate ABCD hypergraph following parameters specified in `params`.

Return a named tuple containing a vector of hyperedges of the graph and a list of cluster
assignments of the vertices.
The ordering of vertices and clusters is in descending order (as in `params`).
"""
function gen_hypergraph(params::ABCDHParams)
    he1 = generate_1hyperedges(params)

    for i in eachindex(params.w, params.z)
        params.z[i] = randround(params.ξ * params.w[i])
    end
    params.y .= params.w .- params.z

    clusters = populate_clusters(params)
    hyperedges = config_model(clusters, params, he1)

    return (hyperedges=hyperedges, clusters=clusters)
end
