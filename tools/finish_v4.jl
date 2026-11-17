# Dokonczenie v4 po pelnym przebiegu (wszystko poza P3/dryfem). Uruchom w Julii na kopierniku:
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/finish_v4.jl"); finish_v4()
# Kroki:
#   1. diagnostyka 6 pulsarow, w ktorych metoda reczna jest lepsza (zap_checks_v4/checks.txt);
#   2. kanaly dla pulsarow pominietych w weryfikacji v4: dlugie obserwacje (N > 4000) bez limitu
#      + J0944-1354 (psradd nie laczyl .spCf16: pelny komunikat bledu);
#   3. jesli wyciecie kanalow pomaga (>= 1% po automacie) -> automat ponownie dla tych pulsarow
#      z kanalami wycietymi przed sumowaniem (+ odswiezony pulsar_rfi.spCF);
#   4. raport modelu nulli z oznaczeniem "model nie rozroznia" (z cache, bez dopasowan);
#   5. porownanie automat vs metoda reczna.
# Nic w pazur nie jest zmieniane; w pazur2 tylko pulsary z kroku 3 (jesli jakies).

using Printf
FV_ROOT = "/home/psr/data/OUTPUT/"
FV_MANUAL_BETTER = ["J1632-4621", "J1843+2024", "J1903+2225", "J1905+0616", "J1918+1444", "J2136-1606"]
FV_LONG = ["J1359-6038", "J1917+1737", "J1932+2220", "J1048-5832", "J1141-6545", "J1057-5226"]

function _fv_dedupe_summary(f)
    isfile(f) || return
    L = readlines(f)
    last_ = Dict{String,Int}()
    for (i, l) in enumerate(L[2:end]); last_[split(l, ",")[1]] = i + 1; end
    keep = [1; sort(collect(values(last_)))]
    length(keep) == length(L) && return
    write(f, join(L[keep], "\n") * "\n")
    println("zap_summary.csv: usunięte stare wiersze ($(length(L) - length(keep)))")
end

function _fv_psradd_error(psr; dataroot="/home/psr/data/new/")
    d = joinpath(dataroot, psr)
    obs = sort(filter(x -> isdir(joinpath(d, x)), readdir(d)))
    isempty(obs) && return println("$psr: brak obserwacji")
    fs = sort(filter(x -> endswith(x, ".spCf16"), readdir(joinpath(d, obs[1]); join=true)))
    println("$psr: $(length(fs)) plików .spCf16 w $(obs[1])")
    err = IOBuffer()
    tmp = joinpath(mktempdir(), "t.ar")
    ok = success(pipeline(`psradd -o $tmp $fs`; stdout=devnull, stderr=err))
    println("  psradd: ", ok ? "działa (błąd był chwilowy)" : "BŁĄD:")
    ok || println("  ", join(last(split(String(take!(err)), '\n'; keepempty=false), 6), "\n  "))
    # ktory plik jest zly: kazdy osobno przez psredit
    for f in fs
        success(pipeline(`psredit -c nsubint $f`; stdout=devnull, stderr=devnull)) || println("  nie czyta się: $f")
    end
end

function finish_v4()
    t0 = time()
    println("\n===== 1. DIAGNOSTYKA: metoda ręczna lepsza (6 pulsarów) =====")
    Base.include(Main, "modules/zap_verify.jl")
    Base.include(Main, "tools/zap_checks.jl")
    Base.invokelatest(Main.zap_checks; verify_dir=FV_ROOT * "zap_verify_v4/", out=FV_ROOT * "zap_checks_v4/",
                      psrs=FV_MANUAL_BETTER)

    println("\n===== 2. KANAŁY: pulsary pominięte w weryfikacji v4 =====")
    _fv_psradd_error("J0944-1354")
    vl = FV_ROOT * "zap_verify_v4_long/"
    Base.invokelatest(Main.zap_verify; out=vl, only=vcat(FV_LONG, ["J0944-1354"]), save_rfi=false,
                      max_chan_pulses=10^6)

    println("\n===== 3. AUTOMAT ponownie, jeśli kanały pomagają =====")
    L = SpaTs.chan_first_list(vl)
    if isempty(L)
        println("W żadnym z tych pulsarów wycięcie kanałów nie podnosi S/N o >= 1% po automacie - nic do zmiany.")
    else
        println("Kanały przed sumowaniem: ", join(["$k ($v)" for (k, v) in L], ", "))
        SpaTs.autozap_all(only=collect(keys(L)), chan_first=L, skip_done=false)
        _fv_dedupe_summary(FV_ROOT * "pazur2/zap_summary.csv")
        Base.invokelatest(Main.zap_verify; out=FV_ROOT * "zap_verify_v4_long_po/", only=collect(keys(L)),
                          chans=false, save_rfi=true)
    end

    println("\n===== 4. MODEL NULLI: raport (z gotowych dopasowań) =====")
    try
        run(addenv(`python3 tools/nulling_model.py`, "NULLING_OUT" => FV_ROOT * "nulling_model_v4/"))
    catch e
        println("model nulli: ", e)
    end

    println("\n===== 5. PORÓWNANIE =====")
    SpaTs.compare_zap_methods()
    @printf("\nKONIEC finish_v4 (%.0f min). Wklej: cat %szap_checks_v4/checks.txt %szap_verify_v4_long/summary.txt\n",
            (time() - t0) / 60, FV_ROOT, FV_ROOT)
end
