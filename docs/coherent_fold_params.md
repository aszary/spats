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

Poprawka (`coherent_mask2.jl`): adaptacyjna nośna (demodulacja przy fazie z poprzedniego przebiegu, 2 iteracje),
próg = mediana nulla, brzegi 1/(2 f_c) pominięte. Maskuje 10–36% pulsów, pokrycie śladu sLRFS 93–100%, brak zer |s|.
Na wspólnych oknach dokładność A ≈ B ≈ C (maska nie poprawia pomiaru, tylko usuwa złe punkty); C lepsze dla J1750
(RMS 6.6 → 5.8) i J1946. Fold bez zmian. J1750 ~830–900: coherent P3 ≈ 35 vs sLRFS ≈ 65 przy minimum |s| —
hipoteza: odwrócenie dryfu (rzut na szablon wybiera drugą wstęgę). Sprawdzone z subtrack (`j1750_reversal_check.jl`):
pasuje dla 40–110 (D ≈ −0.5 °/P, moc w lustrzanej wstędze), nie pasuje dla 800–900 (D ≈ +0.2); statystycznie
nieistotne (AUC 0.67, null 0.49 ± 0.16). Minimum |s| nie jest wiarygodnym wskaźnikiem odwrócenia.

## 5. Filtr Kalmana i wyłączanie nulli (`kalman_compare.jl`)

Kalman (stan faza + częstotliwość, pomiar z(n) bez filtra, parametry z wiarygodności, gładzenie RTS):
- **Fold**: wzrokowo najostrzejsze pasma (J0034, J1750); na przetasowanych pulsach pasm nie ma, więc są zmierzone,
  nie narzucone. J0818: wzór AM widoczny też po tasowaniu (dla wszystkich metod) — częściowo artefakt.
- **P3(n)**: bezużyteczne w tej postaci — wiarygodność wybiera maksymalny szum procesu, błędy nieskalibrowane.
  Dla P3(n) zostaje wariant C.
- Miary ΔR² i Δr dają sprzeczne rankingi metod → do porównań potrzebny benchmark syntetyczny.

Nulle wg energii (ułamek 2·frac(E<0)): lepszy fold dla dryferów z nullami (J0034, J0818), gorszy przy modulacji
natężenia (J1946).

## 6. Benchmark syntetyczny (`synth_benchmark.jl`, 36 przypadków)

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

## Wariant C w skrócie (stan obecny `coherent_fold_agent`)

1. Szablon cyklu dryfu z całej obserwacji (bin FFT przy f3).
2. Każdy impuls rzutowany na szablon → liczba zespolona: kąt = faza w cyklu, długość = wyrazistość.
3. Demodulacja i wygładzenie sąsiednich impulsów, f_c automatyczne (f3/16 … f3/3).
4. Dwa kolejne przebiegi z nośną podążającą za fazą z poprzedniego — zmiany P3 nie są tłumione przez filtr.
5. Z P3(n) wypadają: impulsy o |s| poniżej mediany z danych przetasowanych (`threshold_q`), nulle wykryte z energii
   (`split_nulls`) i brzegi 1/(2 f_c). Fold używa wszystkich impulsów.
6. P3(n) z nachylenia fazy (regresja ważona |s|²) w oknie max(1/(2 f_c), 3·P3, 30 P), osobno w każdym ciągłym odcinku;
   wartości poza [2, 3·P3] odrzucane (pojedynczy punkt −700 w J1919+0134); błędy z
   `n_groups` niezależnych zakresów długości.

## Implementacja (funkcje `_agent`)

- `P3FoldViterbi.auto_cutoff_agent(data, p3, bin_st, bin_end; ybins, grid, nshuffle)` — f_c = max ΔR² (CV po
  długości − null z tasowania) w f3·{1/16, 1/10, 1/8, 1/6, 1/4, 1/3}.
- `P3FoldViterbi.coherent_fold_agent(data, p3, bin_st, bin_end; ybins, lowpass_cutoff=:auto, niter=2, n_groups=4,
  threshold_q=0.5, split_nulls=true)` — fold i P3(n) wariantu C; fold = ŚREDNIA pulsów w binie fazy (nie suma —
  fazy z danych obsadzają biny nierówno, J1539: 31–155 pulsów/bin), zwraca też `counts`, `used`, `nulls`, `amplitude`,
  `lowpass_cutoff`, `cutoff_score`.
- `SpaTs.p3fold_coherent_agent(outdir; datafile, plotdir, name_mod, figtitle, lowpass_cutoff=:auto, threshold_q=0.5,
  split_nulls=true)` — odpowiednik `p3fold_coherent` bez okien, wykres `<name_mod>_p3fold_compare.pdf/.png`.
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

## Benchmark syntetyczny v2 (`synth2_bench.jl`, 200 losowych przypadków)

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

## Wnioski

- Fold: auto f_c zamiast 1/300 — pewny, duży zysk (syntetyki 0.68 → 0.85). Kalman lepiej śledzi fazę, ale fold zyskuje
  mało (+0.01–0.03); na danych rzeczywistych P3(n) z Kalmana bezużyteczne.
- P3(n): wariant C (`coherent_fold_agent`) — mediana błędu ~1.5% na syntetykach, gorzej przy krótkich nullach,
  dużej wędrówce P3 i S/N < 1.5. Ślad sLRFS dokładniejszy tam, gdzie mierzy (~40% pulsów).
- Miary ΔR² / Δr nie nadają się do porównywania metod o różnej sile śledzenia — rozstrzyga benchmark syntetyczny.

## Ewentualne kroki do rozważenia

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
   kandydat do wprowadzenia; σ_k często na granicy siatki 0.15 → rozszerzyć.
7. **Inne** — odróżnienie natężenia od dryfu (κ / widmo energii z progiem per pulsar); osobne P3 dla trybów
   (J1825, J0034); reguła f_c z rozrzutu śladu sLRFS (80–100% optimum w 6/8 pulsarach) jako alternatywa dla CV.
