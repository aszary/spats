# Sliding LRFS: dobór długości okna L

Wydzielone 2026-10-09 z `coherent_fold_params.md` (obecna metoda — fold koherentny — z tego nie korzysta;
sLRFS służy tylko jako niezależne odniesienie dla P3(n)). Kod: `modules/p3track.jl` (`sliding_lrfs`, `p3_track`,
`contrast_null`, `good_windows`, `window_length`). Skrypty `~/claude/work/scripts/j1750_sliding_L.jl`,
`p3track_Ltest.jl`, wykresy `~/claude/work/figures/{sliding_L,Ltest}/`, szczegóły w `separations_analysis_log.md`
(2026-10-08).

## 1. Długość okna L w sliding LRFS (J1750-3503)

`j1750_sliding_L.jl`: L = 4 … 512, P3Track (`sliding_lrfs`, `p3_track`, `contrast_null`, `good_windows`).

- Mierzalne jest tylko P3 ≲ L/3 (fmin = 2/L, osłona krawędzi 1/L). Dla P3 ≈ 49 okna L ≤ 64 nic nie widzą.
- L = 128–150: ślad ucieka na krawędź tam, gdzie lokalne P3 rośnie do 60–65.
- L = 196 (= 4·P3, `window_length`): mediana 48.4, prawie bez krawędzi, widoczne zmiany P3 (36–65).
- L = 384–512: zostają dwie średnie (~53 i ~44). Kompromis: dłuższe L = stabilniej, gorsza rozdzielczość czasowa
  (niezależnych okien ~N/L).
- Reguła praktyczna: **L ≈ 3–4 × najdłuższe lokalne P3**.

## 2. L z najdłuższego P3 w średnim LRFS (10 pulsarów)

`p3track_Ltest.jl`: średnie LRFS (Welch, L0 = 256), istotne piki wobec shuffle 99%, L_new = 4·P3_max.

- Działa dla J1825+0004 (reżim P3 32–60 po ~800 widoczny tylko przy L_new).
- U J0034, J1946 i innych najdłuższy pik to okres nulli/burstów, nie dryfu.
- Kryterium energii (widmo E(n) + koherencja fazy κ = |Σ_φF|²/(Σ_φ|F|)²): rozdziela skrajności
  (dryf J1750 κ ≈ 0.05, nulle J0034/J1946 κ ≈ 0.85), ale stały próg κ = 0.5 zawodzi — tło κ zależy od pulsara,
  strefa szara 0.36–0.48 obejmuje też J1825 (AM).
- Dłuższe L zwiększa liczbę przyjętych okien także dla głównego P3 (4·P3 to za mało cykli dla progu shuffle).

