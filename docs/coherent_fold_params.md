# Fold koherentny (`coherent_fold_agent`): metoda, parametry i test dryfu

**Miejsce startu.** Dokument opisuje obecną metodę (fold koherentny z P3(n)) i bieżący projekt: test hipotezy, że
Song+23 oznaczają dryf także tam, gdzie go nie ma, oraz porównanie z granicą Ė Basu+16. Pełny przebieg: dziennik
`docs/separations_analysis_log.md` (wpisy 2026-10-08 i 2026-10-09). Dobór okna sliding LRFS przeniesiony do
`docs/sliding_lrfs_window.md`.

## 1. Stan i od czego zacząć (2026-10-09)

- **Kod:** `modules/p3fold_viterbi.jl` (`auto_cutoff_agent`, `coherent_fold_agent`, `kernel_fold_agent`,
  `kernel_sigma_agent`), `spats.jl` (`SpaTs.p3fold_coherent_agent` — wykres). Stare `coherent_fold`, `p3fold_coherent`
  bez zmian. Gałąź `claude`.
- **Uruchomienie dla jednego pulsara** (kontener, `psrx`):
  `SpaTs.p3fold_coherent_agent("/home/psr/output/<PSR>[_16]"; datafile=..., plotdir="/home/psr/work/figures/...")`.
  Dla słabych (S/N filtra < 2) lub małego pokrycia P3(n) dodatkowo `threshold_q=0.2, kernel_sigma=0.07`, a gdy
  „nulle” to w rzeczywistości modulacja natężenia — `split_nulls=false`.
- **Partia pulsarów:** `~/claude/work/scripts/highE_agent.jl <lista>` (lista: `katalog plik` w wierszu; fold +
  wariant dla słabych + miara ψ), test nulla: `missing_null.jl` (czyta `~/claude/work/slabe21.txt.tmp`) → wykresy
  `figures/slabe_null/<PSR>.png` (dane / null / dane − null).
- **Wyniki klasyfikacji:** `~/claude/work/coherent_drift_classification.csv` (psr, etykieta Song+23, klasa),
  P–Ṗ: `~/claude/work/scripts/ppdot_coherent_drift.jl` → `figures/ppdot_coherent_drift.png` i `~/output/claude/`.
  Problemy z danymi: `~/claude/work/data_issues.txt`.
- **Ślepy przegląd — ZROBIONY 2026-10-10** (sekcja 6.5): 302 pulsary przeliczone bieżącym kodem z losowymi ID w
  `~/claude/work/blind/raw/` (`scripts/blind_run.jl`, log `logs/blind_run.log`, na końcu „KONIEC”). Werdykty do
  `blind/verdicts.csv` (id,klasa), mapowania `blind/mapping_DO_NOT_READ.csv` nie otwierać przed zapisaniem werdyktów.
  Potem ślepy null dla pochylenie/słabe, na końcu porównanie z poprzednią klasyfikacją, pairshift i Song+23.
- Potem: pozostałe ~120 dryferów Song+23 z Ė ≤ 2·10³² (ślepo), 6 bez danych w `~/output` (J0944-1354, J1402-5021,
  J1605-5257, J1901+0511, J2215+1538, J2346-0609).

## 2. Metoda: wariant C w skrócie (`coherent_fold_agent`)

1. Szablon cyklu dryfu z całej obserwacji (bin FFT przy f3).
2. Każdy impuls rzutowany na szablon → liczba zespolona: kąt = faza w cyklu, długość = wyrazistość.
3. Demodulacja i wygładzenie sąsiednich impulsów, f_c automatyczne (f3/16 … f3/3).
4. Dwa kolejne przebiegi z nośną podążającą za fazą z poprzedniego — zmiany P3 nie są tłumione przez filtr.
5. Z P3(n) wypadają: impulsy o |s| poniżej mediany z danych przetasowanych (`threshold_q`), nulle wykryte z energii
   (`split_nulls`) i brzegi 1/(2 f_c). Fold używa wszystkich impulsów.
6. P3(n) z nachylenia fazy (regresja ważona |s|²) w oknie max(1/(2 f_c), 3·P3, 30 P), osobno w każdym ciągłym odcinku;
   wartości poza [2, 3·P3] odrzucane (pojedynczy punkt −700 w J1919+0134); błędy z
   `n_groups` niezależnych zakresów długości.

## 3. Implementacja (funkcje `_agent`)

- `P3FoldViterbi.auto_cutoff_agent(data, p3, bin_st, bin_end; ybins, grid, nshuffle)` — f_c = max ΔR² (CV po
  długości − null z tasowania) w f3·{1/16, 1/10, 1/8, 1/6, 1/4, 1/3}.
- `P3FoldViterbi.coherent_fold_agent(data, p3, bin_st, bin_end; ybins, lowpass_cutoff=:auto, niter=2, n_groups=4,
  threshold_q=0.5, split_nulls=true)` — fold i P3(n) wariantu C; fold = ŚREDNIA pulsów w binie fazy (nie suma —
  fazy z danych obsadzają biny nierówno, J1539: 31–155 pulsów/bin), zwraca też `counts`, `used`, `nulls`, `amplitude`,
  `lowpass_cutoff`, `cutoff_score`.
- `P3FoldViterbi.kernel_fold_agent(data, phase, sigma; nphase=64)` — fold jądrowy: średnia ważona wszystkich pulsów,
  zawinięty Gauss w fazie o szerokości σ [cykle]; `kernel_sigma_agent(data, p3, bin_st, bin_end, cutoff; grid)` — σ z
  LOO na połówkach binów długości (siatka 0.01–0.3 cyklu, największe σ w 1% od maksimum; LOO na siatce 256 faz, więc
  działa też dla 27 000 pulsów).
- `coherent_fold_agent(...; fold=:kernel, kernel_sigma=:auto, nphase=64)` — DOMYŚLNIE fold jądrowy (`folded` =
  `folded_kernel`, `kernel_sigma`, `sigma_score`); `fold=:bins` → średnia w binach (`folded_bins`, `counts`, `ybins`);
  `kernel_sigma=0.03` — własna szerokość (np. dla drobniejszych szczegółów; LOO wybiera mocne wygładzenie).
- `SpaTs.p3fold_coherent_agent(outdir; datafile, plotdir, name_mod, figtitle, lowpass_cutoff=:auto, threshold_q=0.5,
  split_nulls=true, fold=:kernel, kernel_sigma=:auto, nphase=64)` — odpowiednik `p3fold_coherent` bez okien, wykres `<name_mod>_p3fold_compare.pdf/.png`.
- `coherent_fold` i `p3fold_coherent` bez zmian.

Decyzje projektowe i ich uzasadnienie:
- **Odcinki** — P3(n) bez rozwijania fazy i dopasowania przez przerwę: liczba cykli w przerwie ~P3 jest nieznana
  (J1750: wszystkie skoki > 3 P leżały przy przerwach, wartości obok przerw zależały od założenia do 17.5 P).
- **`threshold_q`** — J1750: 0.5 → pokrycie 72%, 4 przerwy; 0.3 → 81%, 1 przerwa, szersze błędy w słabych odcinkach.
  Na syntetykach 0.3 bez przewagi. Próg wpływa lekko na fold (adaptacyjna nośna interpoluje przez pulsy poniżej progu).
- **`split_nulls`** — energia uśredniona po 5 P, ułamek nulli 2·frac(Ē < 0), epizody ≥ 2 P. Uśrednianie konieczne:
  z pojedynczych pulsów J1750 (S/N energii 1.6, bez nulli) wychodziło fałszywie 18% nulli. Wykrywa: J0034 31%,
  J1750/J0818/J1825 0%. Usuwa piki P3(n) w krótkich nullach (syntetyki: RMS 6.7% → 3.8%); przy skoku fazy po każdym
  nullu pomaga częściowo (15.5% → 12.3%), bo filtr nadal przechodzi przez null.

Testy na danych: f_c J1750 1/147, J0034 1/66, J0818 1/35, J1825 1/142; P3(n) w 60–92% pulsów; 2–10 s na pulsar
(`test_p3fold_coherent_agent*.jl`, `synth_agent_split.jl`).

## 4. Wnioski o metodzie

- Fold: auto f_c zamiast 1/300 — pewny, duży zysk (syntetyki 0.68 → 0.85). Kalman lepiej śledzi fazę, ale fold zyskuje
  mało (+0.01–0.03); na danych rzeczywistych P3(n) z Kalmana bezużyteczne.
- P3(n): wariant C (`coherent_fold_agent`) — mediana błędu ~1.5% na syntetykach, gorzej przy krótkich nullach,
  dużej wędrówce P3 i S/N < 1.5. Ślad sLRFS dokładniejszy tam, gdzie mierzy (~40% pulsów).
- Miary ΔR² / Δr nie nadają się do porównywania metod o różnej sile śledzenia — rozstrzyga benchmark syntetyczny.

## 5. Czy fold pokazuje dane, czy artefakt metody

- **Pochylenie pasm jest w danych:** faza z połowy binów długości + fold drugiej połowy daje ten sam wzór (J0255:
  0.994), lokalne foldy stałym P3 (odcinki ~100 P) pokazują to samo przesunięcie fazy, po odjęciu nulla pochylenie
  zostaje u pewnych dryferów (kontrola J1232).
- **Kontrast zawyżony ~2×:** null = fold jądrowy z przetasowanych pulsów przy WSPÓLNYM szablonie (te same współrzędne
  fazy) daje 25–60% amplitudy foldu, bo faza z dopasowania do szablonu częściowo odtwarza szablon. Kontrastu nie
  traktować jako głębokości modulacji.
- **Falowanie „S” i zygzak przy P3 ≈ 2–3** to w większości artefakt (null odtwarza wzór, korelacja 0.9–0.99).
- **Alias przy P3 ≈ 2–3** (nieokreślony kierunek dryfu) nie jest argumentem przeciw dryfowi — liczy się pochylenie.
- **Złe dane** dają pozorne pochylenie: J0211-8159 — pulsar wędruje w fazie łukiem (zła efemeryda); sprawdzać panel
  pulsów i energii. RFI w pojedynczych pulsach (J1700-3312) trzeba wyzerować i powtórzyć (tam pochylenie zostało).
- J0905-6019, J1742-4616 (S/N < 2): fold płaski, stały fold z modulacją; ψ(φ) płaskie → AM, nie dryf.

## 6. Test hipotezy: Song+23 oznaczają dryf tam, gdzie go nie ma; granica Ė (Basu+16)

### 6.1. Kryteria klasyfikacji (obecne, wzrokowe)

Wykres `*_p3fold_compare.png` w pełnej rozdzielczości: czy plama jasności w kolejnych fazach cyklu P3 przesuwa się
w długości (ukośne pasma) czy tylko zmienia jasność w miejscu (pionowe plamy = AM). Dla podejrzanych — fold − null.
Klasy: `pochylenie` (wyraźne, zostaje po nullu), `umiarkowane`, `slabe_potw` (słabe, zostaje po nullu), `slabe`
(niepewne), `slabe_artefakt` (znika po nullu), `brak_pochylenia`, `zle_dane`.
Pomocniczo miara ψ (`psi_metric.jl`: Δψ wzdłuż składowej przy f3 w 4 odcinkach czasu): stabilne ψ (4/4 znak, CV < 0.5)
to mocne potwierdzenie (pochylenie 48/77), ale jego brak nie wyklucza dryfu (odwrócenia, długie P3, nulle, kilka
składowych: J1750, J1918, J2053). Słabości: subiektywność granic klas, starsze partie oglądane w zmniejszeniu i
starszą wersją kodu, null nie dla wszystkich, ocena NIE była ślepa (w logach z_blk, ψ) → ślepy przegląd (6.5).

### 6.2. Wyniki (przed ślepym przeglądem)

Dryfery Song+23 wg Ė (n = sklasyfikowane z danymi):

| Ė (erg/s) | n | pochylenie (+umiark.) | słabe potw. | słabe | artefakt | brak |
|---|---|---|---|---|---|---|
| ≤ 2·10³² | 132* | 69 | 2 | 4 | 3 | 54 |
| 2·10³²–10³³ (komplet) | 75 | 17 (23%) | 5 | 5 | 8 | 40 |
| > 10³³ (komplet z danymi) | 84 | 2 (2%) | 4 | 6 | 11 | 61 |

\* obciążone: zawiera pulę 104 „pewnych dryferów” wybraną wg pairshift (|z_blk| ≥ 3–5, P3Track). Losowa partia 10
z pozostałych (Ė ≤ 2·10³²): pochylenie 3/10 (J0959-4809, J1840-0809, J1700-3312), J0211-8159 — złe dane.
Wniosek wstępny: dryf zanika stopniowo z Ė, nie urywa się przy 2·10³²; powyżej 10³³ prawie go nie ma. Dryf nad
granicą Basu: J1537-4912 (log Ė 33.45, znany bi-drifter), J1918+1444 (33.71), J1922+1733 (34.60, prawdopodobny);
niepewne J1000-5149, J2043+2740 (~11–12 cykli P3). Wykres P–Ṗ z linią Ė = 2·10³² (Basu, Mitra & Melikidze 2016).

### 6.3. Pulsary analizowane szczegółowo (`cand4_drift.jl`: null, połówki, ψ w odcinkach i oknach, lokalne foldy, pulsy)

- J1537-4912 — bi-drift: główna składowa (163–177°) Δψ −0.19 ± 0.02, słaba (183–193°) +0.17 ± 0.03 (P3Track, stabilne w ćwiartkach; „V” po odjęciu nulla).
- J1922+1733 — pochylone pasmo zostaje po nullu, ψ spójne 3/4.
- J2043+2740, J1000-5149 — niepewne (krótkie obserwacje; J1000: skok fazy ~π między składowymi).
- J1312-5516 — graniczny: ψ niespójne w czasie (pochylenie z jednego odcinka), pairshift −3.3; „slabe”.
- J1016-5345 (P3-only) — słaby/niejednoznaczny.

### 6.4. Kandydaci na dryf wśród P3-only

Lista `~/claude/work/drift_candidates_p3only.txt`. J1016-5345 — słaby/niejednoznaczny: w foldzie koherentnym lekkie
pochylenie (~3°/cykl P3) zostaje po odjęciu nulla z tasowania, ale nachylenie fazy w oknach 128 P jest nieistotne
i zmienia znak (szczegóły w logu 2026-10-09). Pozostałe 9 sprawdzonych P3-only: AM.


### 6.5. Ślepy przegląd (zrobiony 2026-10-10)

302 pulsary pod losowymi ID (`~/claude/work/blind/`), ocena folda bez z_blk, ψ, P3Track i Song+23; potem przeliczenie
po RFI (`scripts/blind_rfi.jl`) i null dla pochylenie/umiarkowane/słabe (`scripts/blind_null.jl`), dopiero na końcu
odsłonięcie mapowania (`scripts/blind_compare.py` → `blind/blind_compare.csv`, log `logs/blind_compare.log`).
Klasa końcowa: słabe + null zostaje → `slabe_potw`, słabe + znika → `slabe_artefakt`, pochylenie/umiarkowane + znika →
`artefakt`, pochylenie + słabnie → `umiarkowane`. 4 złe dane (J dla B6000, B6702, B8146, B9364) bez nulla.

Dryfery Song+23 wg Ė, ślepo (dryf = pochylenie + umiarkowane + słabe potw.; bez złych danych):

| Ė (erg/s) | n | pochylenie | umiark. | słabe potw. | słabe | słabe artef. | artefakt | brak | dryf |
|---|---|---|---|---|---|---|---|---|---|
| ≤ 2·10³² | 130 | 45 | 27 | 4 | 2 | 8 | 9 | 35 | 76 (58%) |
| 2·10³²–10³³ | 75 | 4 | 6 | 2 | 7 | 13 | 6 | 37 | 12 (16%) |
| 10³³–10³⁴ | 64 | 1 | 1 | 2 | 5 | 11 | 2 | 42 | 4 (6%) |
| > 10³⁴ | 18 | 0 | 1 | 0 | 1 | 4 | 0 | 12 | 1 (6%) |

Zgodność z poprzednią (nieślepą) klasyfikacją: z 84 „pochylenie” ślepo dryf 71 (50 pochylenie), 7 artefakt; z 166
„brak” ślepo dryf 14 (głównie umiarkowane/słabe potw. po nullu). Wniosek z 6.2 się utrzymuje: udział dryfu spada
stopniowo z Ė, powyżej 10³³ prawie go nie ma. Pierwszy przedział nadal obciążony pulą wybraną wg pairshift.
Uwaga: 10/18 „artefakt” ma |z_blk| ≥ 3 (3 z ψ 4/4: J0856-6137, J1042-5521, J1819-0925) — null przy wspólnym
szablonie odtwarzający 80–90 % pochylonego wzoru może usuwać też prawdziwy dryf; J1700-3312 ślepo „artefakt”, a w
6.3 pochylenie zostawało po nullu (inne okno pulsów). Te przypadki do sprawdzenia innym nullem (np. bez wspólnego
szablonu). Do zrobienia: P–Ṗ z klasą ślepą.

## 7. Historia testów metody (szczegóły w dzienniku)

### H1. `P3FoldViterbi.coherent_fold`: lowpass_cutoff f_c

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

### H2. Maskowanie pulsów o małej amplitudzie (`coherent_mask.jl`)

Maska |s(n)| ≥ q95 z nulla shuffle; P3(n) z regresji fazy ważonej |s|². Skoki P3(n) (J1750 ~100, J0034 przy nullach)
to przejścia |s| przez zero — maska je usuwa: RMS względem sLRFS spada w 6/6 pulsarach (J1750 11.0 → 6.5),
wartości poza zakresem sLRFS znikają (J1946 11% → 0). Wady: maskuje 26–74% pulsów i wycina też odcinki, gdzie lokalne
P3 odchodzi od f3 (tłumienie filtra, J1750 ~800–890). Brzeg `filtfilt` daje osobny artefakt amplitudy (J1001, pierwsze pulsy).

Poprawka (`coherent_mask2.jl`): adaptacyjna nośna (demodulacja przy fazie z poprzedniego przebiegu, 2 iteracje),
próg = mediana nulla, brzegi 1/(2 f_c) pominięte. Maskuje 10–36% pulsów, pokrycie śladu sLRFS 93–100%, brak zer |s|.
Na wspólnych oknach dokładność A ≈ B ≈ C (maska nie poprawia pomiaru, tylko usuwa złe punkty); C lepsze dla J1750
(RMS 6.6 → 5.8) i J1946. Fold bez zmian. J1750 ~830–900: coherent P3 ≈ 35 vs sLRFS ≈ 65 przy minimum |s| —
hipoteza: odwrócenie dryfu (rzut na szablon wybiera drugą wstęgę). Sprawdzone z subtrack (`j1750_reversal_check.jl`):
pasuje dla 40–110 (D ≈ −0.5 °/P, moc w lustrzanej wstędze), nie pasuje dla 800–900 (D ≈ +0.2); statystycznie
nieistotne (AUC 0.67, null 0.49 ± 0.16). Minimum |s| nie jest wiarygodnym wskaźnikiem odwrócenia.

### H3. Filtr Kalmana i wyłączanie nulli (`kalman_compare.jl`)

Kalman (stan faza + częstotliwość, pomiar z(n) bez filtra, parametry z wiarygodności, gładzenie RTS):
- **Fold**: wzrokowo najostrzejsze pasma (J0034, J1750); na przetasowanych pulsach pasm nie ma, więc są zmierzone,
  nie narzucone. J0818: wzór AM widoczny też po tasowaniu (dla wszystkich metod) — częściowo artefakt.
- **P3(n)**: bezużyteczne w tej postaci — wiarygodność wybiera maksymalny szum procesu, błędy nieskalibrowane.
  Dla P3(n) zostaje wariant C.
- Miary ΔR² i Δr dają sprzeczne rankingi metod → do porównań potrzebny benchmark syntetyczny.

Nulle wg energii (ułamek 2·frac(E<0)): lepszy fold dla dryferów z nullami (J0034, J0818), gorszy przy modulacji
natężenia (J1946).

### H4. Benchmark syntetyczny (`synth_benchmark.jl`, 36 przypadków)

Pasma dryfu ze znaną fazą: P3 ∈ {6, 15, 45}, P3(n) stałe / wędrówka ±20% / skok ×1.3, nulle 0 / 30%, S/N 1.5 / 5.

| metoda | fold (korelacja z prawdą) | faza R | P3(n) RMS / pokrycie |
|---|---|---|---|
| A f_c = 1/300 | 0.68 | 0.55 | — |
| A f_c auto | 0.85 | 0.93 | — |
| C auto | 0.85 | 0.93 | 0.056 / 94% |
| Kalman (−nulle) | 0.86 (0.88) | 0.96 | q_θ = 0: 0.063 / 100%, nieskalibrowany |
| sLRFS | — | — | 0.032 / 39% |

Auto f_c to główny zysk; Kalman lepiej śledzi fazę, ale fold zyskuje mało. Na syntetykach Kalman działa, na danych
nie — generator jest uboższy niż dane (do kalibracji).

### H5. Benchmark syntetyczny v2 (`synth2_bench.jl`, 200 losowych przypadków)

P3 3–50, stałe / wędrówka 5–25% / skok ×0.7–1.4, nulle brak / krótkie / długie (skok fazy po nullu z p 0–1),
S/N 0.8–8, różna geometria podpulsów, 15% AM.

| grupa | mediana błędu P3(n) | pokrycie | faza R |
|---|---|---|---|
| wszystkie | 1.5% | 70% | 0.87 |
| stałe / wędrówka / skok | 0.5% / 2.4% / 1.7% | 69 / 64 / 76% | 0.96 / 0.74 / 0.88 |
| bez nulli / krótkie nulle | 1.1% / 3.1% | 85 / 63% | 0.92 / 0.81 |
| S/N < 1.5 | 2.1% | 61% | 0.75 |
| AM | 3.4% | 53% | 0.76 |

Punkty z błędem > 20%: 4%. Skok P3 (poziomy przed/po w ±5%): 61/78. Błędy σ: w ±2σ 77% (powinno 95%).

## 8. Ewentualne kroki do rozważenia

0. **Przegapiony skok P3** (przegląd wszystkich 200 przypadków): gdy Δf skoku > auto f_c, P3(n) zostaje przy starym
   P3 bez żadnego sygnału błędu (#156: 25.8 zamiast 18.9 przez ~400 P; #153, #110). Najpilniejsze — wynik
   pewny i fałszywy. Wyższe f_c dla P3(n) sprawdzone (`synth2_fcp3.jl`): f3/3 z oknem 3·P3 ≈ obecne (przegapione 11
   zamiast 13 z 78 skoków, gorzej przy stałym P3), f3/2 — 9, ale więcej katastrof; #156 (skok 0.36 f3) nadal źle.
   Do sprawdzenia: nośna startowa z lokalnego P3 (ślad sLRFS) zamiast stałego f3; test zgodności P3(n) ze sLRFS.
1. **Okno P3(n) niezależne od f_c** — ograniczyć okno regresji i brzegi do ~2–3·P3. Teraz przy auto f_c = f3/16
   (stałe P3) okno i brzegi to 8·P3; przy częstych nullach pokrycie spada do zera (syntetyk #161). Odwrotnie przy
   małym P3 i f_c = f3/3: okno 1.5·P3 ≈ 5 P i P3(n) szarpie się z puls na puls (#78, #92, #199) — częściowo
   naprawione (okno ≥ 3·P3); J1041/J1056 (P3 ≈ 4) nadal poszarpane w krótkich odcinkach. Sprawdzone
   (`smallp3_window.jl`, `window_minpulses.jl`): okno ≥ 30 P zmniejsza błąd dla P3 < 8 z 1.6% do 0.85% bez wpływu na
   większe P3 i skoki — wprowadzone; `threshold_q` 0.3 daje +4% pokrycia (J1041/J1056: +20%) przy podobnym błędzie —
   domyślnie zostaje 0.5, 0.3 jako opcja dla słabych pulsarów / małego P3.
2. **Osobne f_c dla P3(n)** — auto f_c optymalizuje fold; przy wędrówce P3 ≳ 20% potrzeba f_c ≈ f3/3, a dobór
   wybierał mniej w 6/56 przypadków (faza gubiona, #58). Np. f_c dla P3(n) ≥ rozrzut P3 z pierwszego przebiegu lub ze
   śladu sLRFS.
3. **Kalibracja błędów σ** — niedoszacowane ~1.5× (w ±2σ 77%); empiryczny współczynnik lub inna metoda.
4. **Osobna demodulacja w odcinkach między nullami** — żeby skok fazy po nullu nie przechodził przez filtr; sens
   tylko dla odcinków ≳ 4·P3.
5. **AM z głęboką modulacją** — fazy niskiej jasności wykrywane jako nulle (syntetyk #20: 40%); opisać lub
   odróżnić okresowe nulle od przypadkowych.
6. **Generator syntetyków** — dopasować do danych (LRFS, S/N, rozrzut fazy J0034/J1750), bo Kalman działa na
   syntetykach, a na danych nie.
6b. **Piki P3(n) na początku odcinków** (J1137 ~23, J1655 ~19, J1056 ~9). Odrzucanie punktów z jednostronnym oknem
   sprawdzone i odrzucone — zbyt duża utrata pokrycia (J1056 69 → 25%), piki tylko częściowo znikają.
6c. **ybins dla foldu koherentnego** — w foldzie koherentnym ybins to rozdzielczość w fazie cyklu, niezależna od P3.
   p3_ybins z params niejednorodne (find_ybins, find_ybins_old, ręczne). Sprawdzone (`ybins_kernel.jl`): CV (LOO,
   połówki binów) daje ybins 4–16, prawie jak wyrocznia (syntetyki 0.761 vs 0.766; 2·P3: 0.677). Fold JĄDROWY
   (zawinięty Gauss w fazie, σ z CV) lepszy w 40/40 syntetyków (0.813) i ma wyższy LOO na 7/7 pulsarach —
   WPROWADZONY jako domyślny (siatka σ do 0.3 cyklu, własne σ opcjonalnie).
6d. **Czy fold pokazuje dane, czy artefakt metody** (J0255 i 5 innych; `j0255_check.jl`, `null_subtract.jl`).
   Pochylenie pasm jest w danych: faza z połowy binów długości + fold drugiej połowy daje ten sam wzór (0.994),
   lokalne foldy stałym P3 (odcinki 100 P) pokazują to samo przesunięcie fazy, po odjęciu nulla pochylenie zostaje
   u wszystkich dryferów. Kontrast jest zawyżony ~2×: null (tasowane pulsy, ten sam szablon) daje 25–60% amplitudy
   foldu, bo faza z dopasowania do szablonu częściowo odtwarza szablon. Kontrastu nie traktować jako głębokości
   modulacji. Przy P3 ≈ 2–3 kierunek/tempo dryfu niejednoznaczne (alias). Ewentualne rozszerzenie
   `coherent_fold_agent` o te testy (opcjonalnie): fold-null z tasowania (wspólny szablon) i fold dane − null,
   nadwyżka amplitudy nad nullem, fold z połówek binów jako kontrola; koszt ~10 dodatkowych przebiegów na pulsar.
6e. **Fallback do stałego P3** — przy S/N ≲ 2 fold koherentny bywa płaski (σ_k na granicy siatki), a fold stałym P3
   pokazuje modulację (J0905-6019, J1742-4616). Porównać LOO foldu koherentnego i stałego P3, wybrać lepszy.
6f. **Próg |s| zawyżony przy modulacji natężenia** — J1703-3241: null z tasowania (36.6) > mediana |s| danych (23.1),
   pokrycie P3(n) 2% przy S/N 11.7. Tasowanie przenosi moc fluktuacji energii do pasma f3. Rozważyć null z pulsami
   znormalizowanymi do energii (lub z odjętą składową energii) przy liczeniu progu.
7. **Inne** — odróżnienie natężenia od dryfu (κ / widmo energii z progiem per pulsar); osobne P3 dla trybów
   (J1825, J0034); reguła f_c z rozrzutu śladu sLRFS (80–100% optimum w 6/8 pulsarach) jako alternatywa dla CV.
