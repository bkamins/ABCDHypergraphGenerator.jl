using Pkg
Pkg.instantiate()

using Random
using ArgParse
using DataFrames
using StatsBase

# example call (assuming that the code is run in /utils folder):
# julia --project abcdh.jl -n 100000 -d 2.5,5,100 -c 1.5,1000,10000 -x 0.5 -q 0.0,0.4,0.3,0.2,0.1 -w :linear -s 1234 -o result --stats
#
# for verbose help run:
# julia --project abcdh.jl --help

include("../src/pl_sampler.jl")
include("../src/hypergraph_sampler.jl")

function parse_commandline()
    s = ArgParse.ArgParseSettings()

    add_arg_group!(s, "required arguments", "REQ");

    ArgParse.@add_arg_table! s begin
        "-n"
            help = "number of vertices"
            required = true
            group = "REQ"
        "-d"
            help = "either a tuple γ,δ,D or file with degree sequence " *
                   "(filename cannot contain two commas)"
            required = true
            group = "REQ"
        "-c"
            help = "either a tuple β,s,S or file with community size sequence " *
                   "(filename cannot contain two commas)"
            required = true
            group = "REQ"
        "-x"
            help = "mixing parameter ξ"
            required = true
            group = "REQ"
        "-q"
            help = "either a sequence q₁,q₂,...,qₖ of weights of hyperedges of sizes " *
                   "from 1 to k or a file name with such a sequence in the same format " *
                   "(if a single number is passed it is considered to be a filename)"
            required = true
            group = "REQ"
        "-w"
            help = "either one of values: ':strict', ':linear', ':majority' or a file " *
                   "name with weights for w_cd where in i-th line are comma separated " *
                   "weights for hyperedge of size i+2 ranging from floor(i+2,2)+1 to i; " *
                   "the :strict value assumes that only c=d weight is non zero, " *
                   "the :linear value assumes weight equal to c, " *
                   "the :majority value assumes all weights are 1"
            required = true
            group = "REQ"
        "-s"
            help = "seed value for generator"
            required = true
            group = "REQ"
        "-o"
            help = "prefix for output file names; the generated file names are " *
                   "[prefix]_deg.txt for degree sequence, [prefix]_comm.txt for " *
                   "community size sequence, [prefix]_assign.txt for assignment " *
                   "of vertices to communities, [prefix]_he.txt for hyperedges"
            group = "REQ"
            group = "optional"
        "-m"
            help = "if this flag passed a multi-hypergraph is generated; by default " *
                   "a simple hypergraph is generated"
            action = :store_true
            group = "optional"
        "--stats"
            help = "if this flag passed print generated hypergraph statistics"
            action = :store_true
            group = "optional"
    end
    s.usage = "julia --project abcdh.jl -n N -d D -c C -x X -q Q -w W -s S [-o O] [-m] [-h] [--stats]"
    return ArgParse.parse_args(s)
end

function main(n, dss, css, x, q, ws, seed; m=false, stats=false, prefix=nothing)

    seed = Int(seed)
    isnothing(seed) && throw(ArgumentError("seed must be an integer"))
    Random.seed!(seed)

    n = Int(n)
    isnothing(n) && throw(ArgumentError("Number of vertices must be an integer"))
    n > 0 || throw(ArgumentError("Number of vertices must be positive"))
    
    if length(dss) == 3 
        γ = Float64(dss[1])
        isnothing(γ) && throw(ArgumentError("γ must be a number"))
        2 < γ < 3 || throw(ArgumentError("γ must be in (2, 3)"))
        δ = Int(dss[2])
        isnothing(δ) && throw(ArgumentError("Number of vertices must be an integer"))
        D = Int(dss[3])
        isnothing(D) && throw(ArgumentError("Number of vertices must be an integer"))
        0 < δ <= D || throw(ArgumentError("Condition 0 < d <= D not met"))
        degs = sample_degrees(γ, δ, D, n)
    end

    if length(css) == 3 
        β = Float64(css[1])
        isnothing(β) && throw(ArgumentError("β must be a number"))
        1 < β < 2 || throw(ArgumentError("β must be in (1, 2)"))
        s = Int(css[2])
        isnothing(s) && throw(ArgumentError("Number of vertices must be an integer"))
        S = Int(css[3])
        isnothing(S) && throw(ArgumentError("Number of vertices must be an integer"))
        δ <= s <= S || throw(ArgumentError("Condition δ <= s <= S not met"))
        coms = sample_communities(β, s, S, n, 1000)
    end

    n != sum(coms) && throw(ArgumentError("number of vertices does not match the sum of community sizes"))

    ξ = Float64(x)
    isnothing(ξ) && throw(ArgumentError("ξ must be a number"))
    0 <= ξ <= 1 || throw(ArgumentError("ξ must be in [0, 1]"))

    q = Float64.(q)

    w = zeros(Float64, length(q), length(q))

    if ws == ":strict"
        for d in 1:length(q)
            w[d, d] = 1.0
        end
    elseif ws == ":linear"
        for d in 1:length(q)
            for c in div(d, 2)+1:d
                w[c, d] = c
            end
            w[:, d] ./= sum(w[:, d])
        end
    elseif ws == ":majority"
        for d in 1:length(q)
            for c in div(d, 2)+1:d
                w[c, d] = 1.0 / (d-div(d, 2))
            end
        end
    end

    hparams = ABCDHParams(degs, coms, ξ, q, w, !m, 100)
    hyperedges, clusters = gen_hypergraph(hparams)

    if prefix !== nothing
        degree_out = "$(prefix)_deg.txt"
        community_out = "$(prefix)_comm.txt"
        assignment_out = "$(prefix)_assign.txt"
        hyperedge_out = "$(prefix)_he.txt"

        open(degree_out, "w") do io
            for d in degs
                println(io, d)
            end
        end

        open(community_out, "w") do io
            for c in coms
                println(io, c)
            end
        end

        open(assignment_out, "w") do io
            for c in clusters
                println(io, c)
            end
        end

        open(hyperedge_out, "w") do io
            for h in hyperedges
                println(io, join(h, ","))
            end
        end
    else
        @info "skipping saving generated graph"
    end

    if stats
        println()
        @info "Degrees"
        dg = zeros(Int, length(degs))
        for he in hyperedges
            for v in he
                dg[v] += 1
            end
        end
        df1 = DataFrame(degs_wanted=degs, degs_generated=dg)
        res = combine(groupby(df1, 1:2), nrow)
        filter!(row -> row[1] != row[2], res)
        println("generated degree distribution")
        describe(dg)
        if nrow(res) > 0
            println("deviations from wanted degrees")
            println(res)
        end

        println()
        println()
        @info "Communities"
        df2 = combine(groupby(DataFrame(x=clusters, sort=true), :x), nrow => :coms_generated, keepkeys=false) |> x -> insertcols!(x, 1, :coms_wanted => coms)
        res = combine(groupby(df2, 1:2), nrow)
        filter!(row -> row[1] != row[2], res)
        println("generated community size distribution")
        describe(coms)
        if nrow(res) > 0
            println("deviations from wanted community sizes")
            println(res)
        end

        println()
        println()
        @info "Hyperedges"
        df3 = leftjoin!(DataFrame(x=1:length(q), q=q), combine(groupby(DataFrame(x=length.(hyperedges)), :x), :x => (x -> first(x)*length(x)) => :vol), on=:x)
        df3.vol = df3.vol ./ sum(skipmissing(df3.vol))
        println(df3)

        @show allunique(hyperedges)
        @show all(allunique, hyperedges)

        hyper_coms = [classify_cluster(clusters[h]) for h in hyperedges]
        ξ_emp = count(x -> x[1] == 0, hyper_coms) / length(hyper_coms)
        @show ξ, ξ_emp

        hyper_sizes = unique(getindex.(hyper_coms, 3))

        xi_df = DataFrame()
        for s in sort(hyper_sizes)
            hc = filter(x -> x[3] == s, hyper_coms)
            ξ_emp = count(x -> x[1] == 0, hc) / length(hc)
            push!(xi_df, (he_size=s, ξ_emp=ξ_emp))
        end

        println()
        display(xi_df)
        println()

        commu = [c[2:3] for c in hyper_coms if c[1] > 0]

        w_emp = zero(w)
        for (c, d) in commu
            w_emp[c, d] += 1
        end
        for i in axes(w_emp, 2)
            w_emp[:, i] ./= sum(w_emp[:, i])
        end

        @info "w"
        display(w)
        println()
        @info "w_emp"
        display(w_emp)
        println()
    end

    return hyperedges, clusters
end

function classify_cluster(clu)
    m = countmap(clu)
    k = collect(keys(m))
    v = collect(values(m))
    am = argmax(v)
    if 2 * v[am] > length(clu)
        return (k[am], v[am], length(clu))
    else
        return (0, 0, length(clu))
    end
end
