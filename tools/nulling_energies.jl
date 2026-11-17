# Energie impulsow do modelu mieszaniny nulli (tools/nulling_model.py).
# Te same dane co automat (_debased_spcf z zap_verify.jl, kanaly zap_chans16), ale
# linia bazowa liczona OD NOWA z okna B rozlacznego z oknem E_off (W): po odjeciu
# sredniej z tego samego okna E_off == 0 i model szumu sie degeneruje (pmod -debase
# bierze cale off-pulse, wiec W musi byc wyjete z linii bazowej).
#   E_on  = suma (I - baza_B) w oknie on-pulse (n_on binow)
#   E_off = suma (I - baza_B) w oknie W (n_on binow, najdalej od impulsu)
# Etykiety automatu: 0 = zostaje, 1 = RFI (zaps2_rfi), 2 = null (zaps2_null).
# Nic nie zmienia w pazur2. Wyniki: /home/psr/data/OUTPUT/nulling_model/energies/<psr>.csv
#
#     cd ~/software/spats && SPATS_NO_MAIN=1 julia --project=. tools/nulling_energies.jl
using Statistics, Printf
if !isdefined(Main, :SpaTs)
    ENV["SPATS_NO_MAIN"] = "1"
    include(joinpath(@__DIR__, "..", "spats.jl"))
end
isdefined(Main, :_debased_spcf) || Base.include(Main, joinpath(@__DIR__, "..", "modules", "zap_verify.jl"))

const NE_OUT = "/home/psr/data/OUTPUT/nulling_model/"
const NE_PAZUR2 = "/home/psr/data/OUTPUT/pazur2/"
const NE_DATA = "/home/psr/data/new/"

"Pulsary: wszystkie z nullami w zap_summary.csv + `ncontrol` kontrolnych bez nulli (NF powinno wyjść ~0)."
function ne_pulsars(; ncontrol=20)
    nul, ctl = String[], Tuple{String,Float64}[]
    for l in Iterators.drop(eachline(joinpath(NE_PAZUR2, "zap_summary.csv")), 1)
        r = split(l, ",")
        length(r) >= 13 && r[2] in ("ok", "reczny") || continue
        nn = something(tryparse(Int, r[13]), 0)
        s = something(tryparse(Float64, r[6]), 0.0)
        nn > 0 ? push!(nul, r[1]) : (s > 50 && push!(ctl, (r[1], s)))
    end
    # kontrolne: rozlozone po S/N (od slabych do jasnych)
    sort!(ctl, by=last)
    idx = unique(round.(Int, range(1, length(ctl), length=min(ncontrol, length(ctl)))))
    return sort(nul), [ctl[i][1] for i in idx]
end

"Okna: W (E_off) = n_on binow off-pulse najdalej od impulsu, B (linia bazowa) = reszta off-pulse bez W±gap."
function ne_windows(D, bst, bend; gap=10)
    nb = size(D, 2)
    on = bst:bend
    off = _VD.offpulse_bins(D, bst, bend)
    c = (bst + bend) / 2
    dist(j) = min(abs(j - c), nb - abs(j - c))
    srt = sort(off, by=j -> -dist(j))
    length(srt) >= 2 * length(on) + 2gap || return nothing
    W = sort(srt[1:length(on)])
    near = Set(mod1(j, nb) for w in W for j in w-gap:w+gap)
    B = filter(j -> !(j in near), setdiff(off, W))
    length(B) >= max(50, length(on) ÷ 2) || return nothing
    return on, W, B
end

function ne_one(psr)
    p = _VT.read_params(joinpath(NE_PAZUR2, psr, "params.json"))
    bst, bend = Int(p["bin_st"]), Int(p["bin_end"])
    obs = sort(filter(d -> isdir(joinpath(NE_DATA, psr, d)), readdir(joinpath(NE_DATA, psr))))
    tmp = mktempdir()
    try
        D = _debased_spcf(joinpath(NE_DATA, psr, obs[1]) * "/", tmp, bst, bend;
                          chans16=something(get(p, "zap_chans16", nothing), ""))
        N = size(D, 1)
        w = ne_windows(D, bst, bend)
        isnothing(w) && return "okno on-pulse za szerokie na dwa rozłączne okna off-pulse"
        on, W, B = w
        lab = zeros(Int, N)
        lab[_VD.zaps_to_rows(_VD.zaps_rfi_of(p), N)] .= 1
        lab[_VD.zaps_to_rows(_VD.zaps_null_of(p), N)] .= 2
        base = [mean(@view D[i, B]) for i in 1:N]
        Eon = [sum(@view D[i, on]) - length(on) * base[i] for i in 1:N]
        Eoff = [sum(@view D[i, W]) - length(on) * base[i] for i in 1:N]
        good = [lab[i] != 1 && any(!iszero, @view D[i, :]) for i in 1:N]
        σbin = std(vec(D[good, W] .- base[good]))
        mkpath(joinpath(NE_OUT, "energies"))
        open(joinpath(NE_OUT, "energies", "$psr.csv"), "w") do io
            println(io, "# pulsar=$psr")
            println(io, "# N=$N")
            println(io, "# n_on=$(length(on))")
            println(io, "# n_B=$(length(B))")
            println(io, "# sigma_bin=$σbin")
            println(io, "# bin_st=$bst")
            println(io, "# bin_end=$bend")
            println(io, "# W=$(first(W))-$(last(W))")
            println(io, "i,E_on,E_off,label")
            for i in 1:N
                any(!iszero, @view D[i, :]) || continue        # puste impulsy (brak danych)
                println(io, "$(i - 1),$(Eon[i]),$(Eoff[i]),$(lab[i])")
            end
        end
        return ""
    finally
        rm(tmp; recursive=true, force=true)
    end
end

function nulling_energies(; part=nothing, only=nothing, ncontrol=20)
    nul, ctl = ne_pulsars(ncontrol=ncontrol)
    psrs = isnothing(only) ? vcat(nul, ctl) : only
    kp, np_ = isnothing(part) ? (1, 1) : parse.(Int, split(part, "/"))
    mkpath(NE_OUT)
    if kp == 1
        open(f -> foreach(c -> println(f, c), ctl), joinpath(NE_OUT, "control.txt"), "w")
    end
    todo = filter(p -> !isfile(joinpath(NE_OUT, "energies", "$p.csv")), psrs[kp:np_:end])
    println("ENERGIE (część $kp/$np_): $(length(todo)) pulsarów"); flush(stdout)
    for (k, psr) in enumerate(todo)
        msg = try
            ne_one(psr)
        catch e
            e isa InterruptException && rethrow()
            "błąd: " * first(sprint(showerror, e), 200)
        end
        isempty(msg) || open(f -> println(f, "$psr,$(replace(msg, ',' => ';'))"), joinpath(NE_OUT, "skipped_$kp.txt"), "a")
        println("[$k/$(length(todo))] $psr $(isempty(msg) ? "ok" : msg)"); flush(stdout)
    end
end

"Równolegle (jak zap_verify_parallel): n procesów w tle, logi w nulling_model/energie_k.log."
function nulling_energies_parallel(; n=nothing, every=60)
    n = _nproc(n)
    codes = ["include(\"tools/nulling_energies.jl\"); Base.invokelatest(Main.nulling_energies; part=\"$k/$n\")" for k in 1:n]
    cnt() = (d = joinpath(NE_OUT, "energies"); isdir(d) ? length(readdir(d)) : 0)
    _run_parallel(codes, NE_OUT, "energie"; every=every, progress=k -> "pliki energii: $(cnt()) | ")
    println("ENERGIE gotowe: $(cnt()) pulsarów w $(joinpath(NE_OUT, "energies"))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    nulling_energies_parallel()
end
