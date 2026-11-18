# Kroki automatu sprawdzane w v7 (test, czy sa potrzebne): nazwa => ustawienia wylaczajace krok.
# Uzywane przez tests/zap_ablation.jl (dane syntetyczne) i tools/final_v7.jl (wszystkie pulsary).
V7_STEPS = [
    ("v6 (wszystkie kroki)", NamedTuple()),
    ("bez pojedynczych śmieci (paski, ACF)", (use_junk=false,)),
    ("bez skoków szumu", (use_jump=false,)),
    ("bez pojedynczych nulli", (use_null1=false,)),
    ("bez wygaszonych odcinków", (use_faded=false,)),
    ("bez zaszumionych odcinków", (use_noisy=false,)),
    ("bez sklejania przerw w nullach", (use_closegaps=false,)),
    ("bez wypełniania krótkich wysp", (min_keep=0,)),
    ("bez bramki S/N", (gate=false,)),
    ("bez bezpiecznika RFI", (rfi_tol=1.0,)),
    ("bez reguły RFI (rms off-pulse)", (rfi=false,)),
    ("najprostszy: RFI rms + nulle + bramka", (use_junk=false, use_jump=false, use_null1=false, use_faded=false,
                                              use_noisy=false, use_closegaps=false, min_keep=0)),
]
_v7_kw(st) = (gate=get(st, :gate, true), rfi=get(st, :rfi, true), rfi_tol=get(st, :rfi_tol, 0.01),
              Base.structdiff(st, (gate=0, rfi=0, rfi_tol=0))...)
