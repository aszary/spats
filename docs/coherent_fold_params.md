# Dobór parametrów: sliding LRFS (L) i `coherent_fold` (lowpass_cutoff)

Stan na 2026-10-08. Skrypty w `~/claude/work/scripts/`, wykresy w `~/claude/work/figures/{sliding_L,Ltest,cohfc}/`,
szczegóły i liczby w `docs/separations_analysis_log.md` (wpisy z 2026-10-08).

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

## 3. `P3FoldViterbi.coherent_fold`: lowpass_cutoff f_c

Demodulacja przy f3 + filtr dolnoprzepustowy f_c to odpowiednik sliding LRFS z L ≈ 1/f_c:
- f_c > rozrzut |1/P3_lok − f3| — inaczej faza nie nadąża za P3;
- f_c ≲ f3/3 — powyżej do fazy przecieka składowa −f3 i modulacja natężenia (nulle).

`coherent_fc_cv.jl`: kroswalidacja po długości (faza z połowy binów, fold na drugiej), miara
ΔR² = R²_cv − null z przetasowanych pulsów. Sam R²_cv rośnie do f_c → f3 (sortowanie pulsów wg podpulsu), stąd null.

| PSR | P3 | optimum f_c (≤ f3/3) | ΔR² opt / 1/300 |
|---|---|---|---|
| J0034-0721 | 6.6 | 0.07 f3 (1/100) | 0.031 / 0.016 |
| J1825+0004 | 14.2 | 0.14 f3 (1/100) | 0.009 / ~0 |
| J0818-3232 | 5.8 | 0.17 f3 (1/35) | 0.069 / 0.011 |
| J1001-5559 | 4.2 | plateau 0.08–0.5 f3 | 0.070 / 0.020 |
| J0959-4809 | 6.0 | 0.33 f3 | 0.026 / ~0 |
| J1750-3503 | 49 | 0.33 f3 (rośnie do f3) | 0.010 / 0.004 |
| J1946+1805 | 19 | 0.33 f3 | 0.032 / ~0 |
| J1626-4537 | 25 | 0.33 f3 | 0.008 / ~0 |
| J1905-0056, J1614+0737 | | brak sygnału | ~0 |

Porównanie foldów (`coherent_fc_compare.jl`): przy 1/300 fold ≈ stałe P3. Przy optimum: J0034 i J1750 — wyraźne pasma
dryfu; J0818, J1001, J1825 — modulacja amplitudy składowych w przeciwfazie (AM); J1946 — modulacja natężenia całego
profilu (nulle), więc wzrost ΔR² nie zawsze oznacza lepszy fold dryfu.

## 4. Maskowanie pulsów o małej amplitudzie (`coherent_mask.jl`)

Maska |s(n)| ≥ q95 z nulla shuffle; P3(n) z regresji fazy ważonej |s|². Skoki P3(n) (J1750 ~100, J0034 przy nullach)
to przejścia |s| przez zero — maska je usuwa: RMS względem sLRFS spada w 6/6 pulsarach (J1750 11.0 → 6.5),
wartości poza zakresem sLRFS znikają (J1946 11% → 0). Wady: maskuje 26–74% pulsów i wycina też odcinki, gdzie lokalne
P3 odchodzi od f3 (tłumienie filtra, J1750 ~800–890). Brzeg `filtfilt` daje osobny artefakt amplitudy (J1001, pierwsze pulsy).

## Wnioski i otwarte sprawy

- Domyślne `lowpass_cutoff = 1/300` w `p3fold_coherent` jest za niskie dla P3 ≲ 15. f_c skalować z f3:
  zakres **0.1–0.33·f3**; stała wartość awaryjna f3/5; najlepiej automatycznie: max ΔR² w [f3/16, f3/3].
- Reguła z rozrzutu śladu sLRFS (f_c = 1.3·q90|Δf|) daje 80–100% optimum w 6/8 pulsarach (J0034 tylko ~50%).
- Do zrobienia: maskowanie tylko dla P3(n), łagodniejszy próg i normalizacja |s| przez tłumienie filtra; odróżnienie
  natężenia od dryfu (κ / widmo energii z progiem kalibrowanym per pulsar); osobne P3 dla trybów (J1825, J0034).
- Kod w `spats.jl` / `modules/` bez zmian — tylko skrypty testowe.
