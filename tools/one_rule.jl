# Test JEDNEJ metody wycinania dla WSZYSTKICH pulsarow (bez mieszania: model w czesci pulsarow, automat w reszcie).
# pazur2 bez zmian; liczy warianty na tych samych danych i porownuje z obecnym v8 i z metoda reczna.
#   A: automat dla wszystkich (bez modelu energii) - jedna metoda, jak v5/v7 bez modelu
#   B: prosta regula dla wszystkich: null = seria >= 5 impulsow bez wykrytej emisji (bez "wygaszen" i pojedynczych
#      pustych impulsow); RFI bez zmian
#   C: B + RFI z calej obserwacji (mediana szumu wszystkich impulsow zamiast 31 sasiadow)
#   ENV["SPATS_NO_MAIN"]="1"; include("spats.jl"); include("tools/one_rule.jl"); one_rule()
# ~20-30 min na 32 procesach, wznawialne. Wynik: final_v7/one_rule/WYNIK.txt
isdefined(Main, :_s8_part) || include(joinpath(@__DIR__, "simple_v8.jl"))
S8_DIR = F7_DIR * "one_rule/"
S8_VARIANTS = [
    ("v8 (obecny): model P>=0.5, bez modelu automat", NamedTuple(), :model_fb, 0.5),
    ("A: automat dla wszystkich", NamedTuple(), :auto, 0.0),
    ("B: null = seria >=5 impulsów bez emisji, dla wszystkich", (use_faded=false, use_null1=false), :auto, 0.0),
    ("C: B + RFI z całej obserwacji", (use_faded=false, use_null1=false, rfi_win=0), :auto, 0.0),
]

function one_rule(; n=32)
    isdefined(Main, :_run_parallel) || Base.include(Main, "modules/zap_verify.jl")
    mkpath(S8_DIR)
    codes = [F7_PRE * "include(\"tools/one_rule.jl\"); Base.invokelatest(_s8_part, $k, $n)" for k in 1:n]
    Base.invokelatest(Main._run_parallel, codes, S8_DIR * "logs", "or"; every=120)
    one_rule_summary()
end

one_rule_summary() = simple_v8_summary(; title="Jedna metoda wycinania dla wszystkich pulsarów (porównanie z v8)", basename_="v8")
