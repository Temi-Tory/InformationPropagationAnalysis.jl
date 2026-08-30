using Test
using InformationPropagationAnalysis
const IPA = InformationPropagationAnalysis

include("oracle.jl")
using .Oracle

const FIX = joinpath(@__DIR__, "fixtures")

# ─────────────────────────────────────────────────────────────────────────────
# helpers
# ─────────────────────────────────────────────────────────────────────────────

"Run the exact reachability-reliability pipeline, return Dict(node => belief)."
function run_belief(edgelist, node_priors::Dict, edge_probs::Dict, sources)
    out, inc = Oracle.build_indices([Tuple{Int,Int}(e) for e in edgelist])
    out = Dict{Int64,Set{Int64}}(k => Set{Int64}(v) for (k, v) in out)
    inc = Dict{Int64,Set{Int64}}(k => Set{Int64}(v) for (k, v) in inc)
    el  = Vector{Tuple{Int64,Int64}}(edgelist)
    src = Set{Int64}(sources)
    itersets, anc, desc = IPA.Input.find_iteration_sets(el, out, inc)
    fk, jn = IPA.Input.identify_fork_and_join_nodes(out, inc)
    r, u = IPA.new_identify(el, node_priors, edge_probs, src, fk, jn, anc, desc, itersets)
    return IPA.update_beliefs_iterative(el, itersets, out, inc, src, node_priors, edge_probs,
                                        desc, anc, r, jn, fk, u)
end

"CPM longest-path project value (makespan) for scalar Float64 durations."
function run_makespan(edgelist, node_dur::Dict{Int,Float64}, edge_del::Dict{Tuple{Int,Int},Float64})
    el = Vector{Tuple{Int64,Int64}}(edgelist)
    out, inc = Oracle.build_indices([Tuple{Int,Int}(e) for e in edgelist])
    out = Dict{Int64,Set{Int64}}(k => Set{Int64}(v) for (k, v) in out)
    inc = Dict{Int64,Set{Int64}}(k => Set{Int64}(v) for (k, v) in inc)
    itersets, _, _ = IPA.Input.find_iteration_sets(el, out, inc)
    src = Set{Int64}(n for n in keys(out) if !haskey(inc, n) || isempty(inc[n]))
    nv = Dict{Int64,Float64}(Int64(k) => v for (k, v) in node_dur)
    ev = Dict{Tuple{Int64,Int64},Float64}(Tuple{Int64,Int64}(k) => v for (k, v) in edge_del)
    for e in el
        haskey(ev, e) || (ev[e] = 0.0)
    end
    r = IPA.CriticalPath.analyze(itersets, out, inc, src, nv, ev; mode = IPA.CriticalPath.LONGEST_PATH)
    return r.project_value
end

read_csv_rows(path) = [split(strip(l), ',') for l in readlines(path)[2:end] if !isempty(strip(l))]

# ─────────────────────────────────────────────────────────────────────────────

@testset "InformationPropagationAnalysis" begin

@testset "public API surface" begin
    for t in (:Input, :Diamonds, :Probability, :CriticalPath, :Flow)
        @test t in names(IPA)
    end
    for v in (:Interval, :pbox, :new_identify, :update_beliefs_iterative, :critical_path, :analyze_all)
        @test v in names(IPA)
    end
    @test IPA.critical_path === IPA.CriticalPath.analyze
    # curation held: Flow's granular helpers are `public`, not exported
    @test Base.ispublic(IPA.Flow, :birnbaum_importance)
    @test !Base.isexported(IPA.Flow, :birnbaum_importance)
    @test Base.isexported(IPA.Flow, :analyze_all)
end

@testset "Probability — exact reliability vs brute force (synthetic diamonds)" begin
    # single diamond: 1 forks to {2,3}, 4 joins {2,3} — parents share fork ancestor 1
    e1 = [(1,2),(1,3),(2,4),(3,4)]
    np = Dict{Int64,Float64}(n => 0.9 for n in 1:4)
    ep = Dict{Tuple{Int64,Int64},Float64}(e => 0.85 for e in e1)
    got = run_belief(e1, np, ep, [1])
    want = brute_reliability([Tuple{Int,Int}(e) for e in e1],
                             Dict{Int,Float64}(np), Dict{Tuple{Int,Int},Float64}(ep), [1])
    for n in 1:4
        @test got[n] ≈ want[n] atol=1e-10
    end

    # nested / overlapping: inner diamond at join 7 (fork 2), outer at join 8 (fork 1)
    e2 = [(1,2),(1,3),(2,4),(2,5),(3,6),(4,7),(5,7),(6,8),(7,8)]
    np2 = Dict{Int64,Float64}(n => 0.92 for n in 1:8)
    ep2 = Dict{Tuple{Int64,Int64},Float64}(e => 0.8 for e in e2)
    got2 = run_belief(e2, np2, ep2, [1])
    want2 = brute_reliability([Tuple{Int,Int}(e) for e in e2],
                              Dict{Int,Float64}(np2), Dict{Tuple{Int,Int},Float64}(ep2), [1])
    for n in 1:8
        @test got2[n] ≈ want2[n] atol=1e-10
    end
end

@testset "Probability — counterexample-n15 (bug-#1 regression guard)" begin
    dir = joinpath(FIX, "counterexample-n15")
    el, out, inc, src = IPA.Input.read_graph_to_dict(joinpath(dir, "counterexample-n15.EDGES"))
    np = IPA.Input.read_node_priors_from_json(joinpath(dir, "float", "counterexample-n15-nodepriors.json"))
    ep = IPA.Input.read_edge_probabilities_from_json(joinpath(dir, "float", "counterexample-n15-linkprobabilities.json"))
    itersets, anc, desc = IPA.Input.find_iteration_sets(el, out, inc)
    fk, jn = IPA.Input.identify_fork_and_join_nodes(out, inc)
    r, u = IPA.new_identify(el, np, ep, src, fk, jn, anc, desc, itersets)
    bel = IPA.update_beliefs_iterative(el, itersets, out, inc, src, np, ep, desc, anc, r, jn, fk, u)
    # exact value from the diamond-rewrite validation campaign (CUDD-verified,
    # validation/validate_broad.jl → 0 wrong). The pre-rewrite hybrid gave 0.66916.
    @test bel[15] ≈ 0.71858899 atol=1e-6
end

@testset "Probability — power-network, all three value forms" begin
    dir = joinpath(FIX, "power-network")
    el, out, inc, src = IPA.Input.read_graph_to_dict(joinpath(dir, "power-network.EDGES"))
    itersets, anc, desc = IPA.Input.find_iteration_sets(el, out, inc)
    fk, jn = IPA.Input.identify_fork_and_join_nodes(out, inc)

    np_f = IPA.Input.read_node_priors_from_json(joinpath(dir, "float", "power-network-nodepriors.json"))
    ep_f = IPA.Input.read_edge_probabilities_from_json(joinpath(dir, "float", "power-network-linkprobabilities.json"))
    r, u = IPA.new_identify(el, np_f, ep_f, src, fk, jn, anc, desc, itersets)
    bel_f = IPA.update_beliefs_iterative(el, itersets, out, inc, src, np_f, ep_f, desc, anc, r, jn, fk, u)

    # sanity: belief ≤ prior, everything in [0,1]
    for n in keys(bel_f)
        @test 0.0 <= bel_f[n] <= np_f[n] + 1e-12
    end

    # golden master (regression guard) — generated by this validated build; see fixtures/README.md
    goldpath = joinpath(dir, "expected-float-beliefs.csv")
    if isfile(goldpath)
        for row in read_csv_rows(goldpath)
            n = parse(Int, row[1]); v = parse(Float64, row[2])
            @test bel_f[n] ≈ v atol=1e-9
        end
    else
        @warn "no golden belief vector yet — run test/gen_golden.jl"
    end

    # cross-version: the v0.1.0 (5-flat-file) release recorded these same beliefs
    v1 = read_json_dict(joinpath(dir, "v0.1.0-beliefs.json"))
    for (k, v) in v1
        @test bel_f[parse(Int, k)] ≈ Float64(v) atol=1e-12
    end

    # interval priors must bracket the float beliefs
    np_i = IPA.Input.read_node_priors_from_json_interval(joinpath(dir, "interval", "power-network-nodepriors.json"))
    ep_i = IPA.Input.read_edge_probabilities_from_json_interval(joinpath(dir, "interval", "power-network-linkprobabilities.json"))
    ri, ui = IPA.new_identify(el, np_i, ep_i, src, fk, jn, anc, desc, itersets)
    bel_i = IPA.update_beliefs_iterative(el, itersets, out, inc, src, np_i, ep_i, desc, anc, ri, jn, fk, ui)
    for n in keys(bel_f)
        @test bel_i[n].lower <= bel_f[n] + 1e-9
        @test bel_i[n].upper >= bel_f[n] - 1e-9
    end
end

@testset "CriticalPath — PSPLIB j301_1 makespan + interval-split regression" begin
    dir = joinpath(FIX, "psplib-j301_1")
    el, _, _, _ = IPA.Input.read_graph_to_dict(joinpath(dir, "j301_1.EDGES"))
    nd, ed = parse_cpm_inputs(joinpath(dir, "float", "j301_1-cpm-inputs.json"))
    # PSPLIB's own published MPM-Time for this instance is 38.
    @test run_makespan(el, nd, ed) ≈ 38.0 atol=1e-9

    # interval-split must not throw on the degenerate node (the interval_analyze_split fix)
    out, inc = Oracle.build_indices([Tuple{Int,Int}(e) for e in el])
    out = Dict{Int64,Set{Int64}}(k => Set{Int64}(v) for (k, v) in out)
    inc = Dict{Int64,Set{Int64}}(k => Set{Int64}(v) for (k, v) in inc)
    itersets, _, _ = IPA.Input.find_iteration_sets(Vector{Tuple{Int64,Int64}}(el), out, inc)
    src = Set{Int64}(n for n in keys(out) if !haskey(inc, n) || isempty(inc[n]))
    ndi, _ = parse_cpm_inputs_interval(joinpath(dir, "interval", "j301_1-cpm-inputs.json"))
    nv = Dict{Int64,IPA.CriticalPath.ValueInterval}(
        Int64(k) => IPA.CriticalPath.ValueInterval(lo, hi) for (k, (lo, hi)) in ndi)
    ev = Dict{Tuple{Int64,Int64},IPA.CriticalPath.ValueInterval}(
        e => IPA.CriticalPath.ValueInterval(0.0, 0.0) for e in Vector{Tuple{Int64,Int64}}(el))
    local res
    @test (res = IPA.CriticalPath.interval_analyze_split(itersets, out, inc, src, nv, ev;
                    mode = IPA.CriticalPath.LONGEST_PATH)) isa IPA.CriticalPath.IntervalPathResult

    # per-node slack intervals match the recorded oracle bounds
    bounds = Dict(parse(Int, r[1]) => (parse(Float64, r[2]), parse(Float64, r[3]))
                  for r in read_csv_rows(joinpath(dir, "expected-float-bounds.csv")))
    for (n, m) in res.margin
        lo, hi = bounds[n]
        @test m.lo ≈ lo atol=1e-6
        @test m.hi ≈ hi atol=1e-6
    end
end

@testset "Flow — genrmf DIMACS instance: solvers agree, max-flow == min-cut" begin
    dir = joinpath(FIX, "genrmf-dag-small")
    el, out, inc, srcset = IPA.Input.read_graph_to_dict(joinpath(dir, "genrmf_dag_small.EDGES"))
    caps = IPA.Input.read_edge_capacities_from_json(joinpath(dir, "genrmf_dag_small-capacities.json"))
    sinks = [n for n in union(keys(out), keys(inc)) if !haskey(out, n) || isempty(out[n])]
    srcs  = collect(srcset)

    rd = IPA.Flow.solve_max_flow_dinic(el, out, inc, caps, srcs, sinks)
    rk = IPA.Flow.solve_max_flow_edmonds_karp(el, out, inc, caps, srcs, sinks)
    rp = IPA.Flow.solve_max_flow_push_relabel(el, out, inc, caps, srcs, sinks)
    @test rd.max_flow ≈ rk.max_flow atol=1e-6
    @test rd.max_flow ≈ rp.max_flow atol=1e-6
    @test rd.max_flow ≈ rd.mincut_capacity atol=1e-6   # max-flow min-cut
end

@testset "Flow — analyze_all end to end (single-diamond capacity network)" begin
    # 1 --10--> 2 --5--> 4 ,  1 --7--> 3 --8--> 4   ⇒  max s-t flow = 5 + 7 = 12
    el = Tuple{Int64,Int64}[(1,2),(1,3),(2,4),(3,4)]
    out, inc = Oracle.build_indices([Tuple{Int,Int}(e) for e in el])
    out = Dict{Int64,Set{Int64}}(k => Set{Int64}(v) for (k, v) in out)
    inc = Dict{Int64,Set{Int64}}(k => Set{Int64}(v) for (k, v) in inc)
    caps = Dict{Tuple{Int64,Int64},Float64}((1,2)=>10.0, (2,4)=>5.0, (1,3)=>7.0, (3,4)=>8.0)
    kit = IPA.analyze_all(el, out, inc, caps, Int64[1], Int64[4])
    @test kit isa IPA.Flow.FlowCapacityResult
    @test IPA.Flow.solve_max_flow_dinic(el, out, inc, caps, Int64[1], Int64[4]).max_flow ≈ 12.0 atol=1e-9
end

end
