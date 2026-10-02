# Analiza separations — dziennik i przewodnik

Wyznaczanie Δseparation (offset składowych profilu między 1023 a 1523 MHz) dla listy
`input/separations_todo.csv`. Wyniki lądują w `input/separations.csv`.
Powiązane dokumenty:
[`profile_narrowing_summary.md`](profile_narrowing_summary.md) — zwięzłe podsumowanie całej
analizy (dane, metoda, wyniki); ten plik to dziennik chronologiczny stojący za nim.
[`component_position_methods.md`](component_position_methods.md).

Plik jest dostępny także jako `~/claude/work/NOTES.md` (symlink) — zgodnie z konwencją
prowadzenia dziennika z `CLAUDE.md`.

---

# START TUTAJ (dryf vs P3-only) — stan na 2026-10-02

Prace nad metodą rozstrzygania, czy pulsar dryfuje, czy jest tylko P3-only (etykiety Song+23), 533 pulsary:

| metoda | opis | wpisy w tym dzienniku | wynik |
|---|---|---|---|
| travel (T_cv, f_trav, ρ) | [`travel_test_method.md`](travel_test_method.md) | 2026-09-22 … 2026-09-30 | niezadowalająca (czułość, f_trav, zgodność z oceną wzrokową) |
| P3Track (sliding LRFS → fold z kompensacją P₃ → faza szablonu) | [`p3track_method.md`](p3track_method.md) — §0, §8 sprawy otwarte, §9 co nie zadziałało | 2026-10-01 (cd. 4) … 2026-10-02 | wobec Song+23 zmienia bardzo niewiele |

Następny krok (2026-10-02): nowa sesja szuka lepszej metody — najpierw przeczytać §0, §8, §9 obu dokumentów.
Wyniki P3Track: `~/output/claude/p3track_batch/` (CSV v1–v4b, wykresy v4); skrypty `~/claude/work/scripts/`.

---

# START TUTAJ (separations) — stan na 2026-09-17

## Gdzie co jest

| co | gdzie |
|---|---|
| wyniki | `~/claude/software/spats/input/separations.csv` (91 wierszy) |
| lista zadań | `~/claude/software/spats/input/separations_todo.csv` (106 pulsarów) |
| dane wejściowe | `~/output/claude/<PSR>_16/` (host) = `/home/psr/output/<PSR>_16/` (kontener) |
| kod | `~/claude/software/spats` = `/home/psr/software/spats` — **bind mount, ten sam katalog** |
| ten dziennik | `docs/separations_analysis_log.md` w repo (symlink: `~/claude/work/NOTES.md`) |
| skrypty pomocnicze | `~/claude/work/scripts/` (poza repo) |
| logi | `~/claude/work/logs/` |
| diagnoza p3-foldów | `~/claude/work/p3fold_diagnoza.md` + `scan_all.csv` |
| backupy | `~/claude/work/backup_p3fold/`, `separations_przed_poprawka.csv` |

Repo: branch `claude`, wszystko wypchnięte na `origin/claude`. `master` nietknięty.

## Stan: LISTA UKOŃCZONA

**106 pulsarów = 91 zapisanych + 15 odrzuconych.** Nic nie czeka na policzenie.

## Jak policzyć nowego pulsara (workflow)

```bash
# 1. kryterium jakości + podgląd wyniku w jednym przebiegu (nic nie zapisuje)
#    edytuj liste `targets` w skrypcie, potem:
psrx julia --project=/home/psr/software/spats /home/psr/work/scripts/batch_auto.jl

# 2. porównaj z todo.csv, obejrzyj "KRYTERIUM: ... odrzucam [...]"

# 3. zapis (psr brany z nazwy katalogu; keep = maska pulsów do zachowania)
psrx julia --project=. -e 'include("modules/data.jl"); k = trues(N); k[i]=false;
  Data.analyse_p3folds_16_new_agent("/home/psr/output/<PSR>_16", "norefine"; n_comp=NC, keep=k)'

# 4. commit + push (pierwszy push w sesji czasem pada - wtedy powtórz sam push)
cd ~/claude/software/spats && git add input/separations.csv && git commit -m "..." \
  && (git push origin claude || psrx git push origin claude)
```

`psr=""` = tryb podglądu (nie zapisuje do `separations.csv`). Bez tego argumentu nazwa
pulsara jest brana z katalogu i wynik **jest zapisywany**.

## Pułapki, które kosztowały najwięcej czasu

1. **Nie ufaj własnej liście kroków — licz wiersze w pliku.** Dwa razy pomyliłem się
   w rachunku postępu; raz cały pulsar (J1919+1745) wypadł z przetwarzania niezauważony,
   bo skrypt kontrolny miał błąd w pętli.
2. **`include("spats.jl")` uruchamia `SpaTs.main()`** — sprawdź, co jest tam
   odkomentowane. Do pojedynczej funkcji używaj `include("modules/data.jl")`.
3. **Ocena wykresów „na oko" jest zawodna i wolna.** Twarde kryterium liczbowe
   (`check_order.jl` / `batch_auto.jl`) wychwyciło błędy, których nie widziałem, i naprawiło
   4 pulsary z przeciwnym znakiem.
4. **Próg absolutny w kryterium jakości był zły** — realny offset bywa duży (J1843-0459
   ma 11 binów i to sygnał). Kryterium jest adaptacyjne: `median + 4·MAD`.
5. **Kryterium adaptacyjne nie wykryje awarii systematycznej** (gdy złe są wszystkie pulsy).
   Sygnałem jest wtedy sama **mediana |Δmu|**: zdrowe pulsary 0.2-4.9 bina, zepsute 16-23.
6. **`component_offsets` sortuje komponenty po pozycji** (`GaussianFit.jl:258`), więc G1 to
   zawsze komponent lewy. Nie ma problemu „odwróconej numeracji" — raz błędnie to zgłosiłem.

## Zadania otwarte

1. **J1527-5552** — zapisany wpis jest **błędny** (artefakt duplikatów p3-foldu w `norefine`).
   `refine` daje +1.175 ± 0.173 przy χ² < 4 zamiast −0.857 ± 0.694. Wymaga decyzji:
   czy przejść na `refine` dla dotkniętych pulsarów (niespójność metody w pliku!).
2. **Ostrzeżenie w kodzie**, gdy `p3_ybins > round(P3)` przy całkowitym P3 — 10 pulsarów
   dotkniętych, mechanizm w `p3fold_diagnoza.md`.
3. **Zweryfikować 3 wyniki rozbieżne z `todo.csv`**: J1824-0127 (+0.960 vs −0.310),
   J1825-0935 (+0.513 vs −0.273), J1919+0134 (−0.599 vs −0.972).
4. **Zdecydować, co z `todo.csv`** — trzy jego wartości okazały się błędne (patrz niżej),
   więc plik nie jest wiarygodnym punktem odniesienia.

## Najważniejszy wynik merytoryczny

**J1733-3716: wartość −12.208 w `todo.csv` to artefakt jednego złego dopasowania.**
Poprawnie **−1.827 ± 0.166** (χ²=0.5). To był największy odstający punkt na diagramie P-Ṗ
z offsetami (`Plot.ppdot_offsets`) — komentarz w `modules/plot.jl` opisywał go jako
„flagged in the input as needing a redo". Podobnie **J1714-1054**: −2.241 → −0.944 ± 0.073.
Oraz **J1757-2421, J1803-3329, J1807-0847** miały w `todo.csv` przeciwny znak.

---

# DZIENNIK CHRONOLOGICZNY


## 2026-09-11 — przegląd pulsarów z input/separations_todo.csv

Workflow: `Data.analyse_p3folds_16_new_agent` (faza 1, `psr=""`) generuje wykresy per-pulse do
przeglądu w `outdir`, potem po ocenie `keep` — faza 2 zapisuje wynik do `input/separations.csv`.
Wszystkie wyniki commitowane i pushowane na branch `claude`.

### Zrobione (21):
J0134-2937, J0151-0635, J0304+1932, J0459-0210, J0729-1836 (pulse 10 skip — zły fit high-freq),
J0738-4042, J0818-3232, J0820-4114, J0837+0610, J0907-5157 (pulse 3 skip — zamiana etykiet G2/G3),
J1001-5559, J1001-5939, J1017-5621 (pulse 1 skip — zamiana etykiet G2/G3), J1041-1942,
J1057-5226, J1110-5637, J1112-6926, J1119-7936, J1123-4844, J1137-6700, J1224-6407.

### POMINIĘTE — wymaga ręcznego przeglądu:
- **J1231-4609** (n_comp=2): systemowy problem, nie pojedynczy zły pulse.
  - Pulse 2 i 3 mają identyczne dane (duplikat w p3-foldzie).
  - Niestabilna identyfikacja komponentów w profilu High: ostry pik profilu jest zawsze przy
    bin~538, ale fit "High G1" czasem łapie ten ostry pik (pulsy 1,4,8), a czasem szeroki,
    słabszy garb przy bin~522 (pulsy 5,6,7) — podczas gdy "Low" konsekwentnie trzyma się
    ostrego piku jako G1. To winduje χ² do 60-108 (dla porównania: inne pulsary miały χ²<30).
  - `todo.csv` miał to jako `dsep=0.006, sigma=0.0, grade=8` (już wcześniej uznane za
    nieistotne/niepewne przez oryginalną ręczną analizę).
  - Nic nie zapisane do `separations.csv`. Do rozważenia: fit z dopasowaniem komponentów po
    pozycji (nie po kolejności/amplitudzie) zamiast obecnego `GaussianFit`/`component_offsets`.

- **J1232-4742** (n_comp=2): ten sam typ problemu co J1231-4609, jeszcze bardziej ekstremalny
  (χ² = 265-657!). "High G1" konsekwentnie przy bin~440-445, "Low G1" przy bin~453-457 — stały
  rozjazd ~10-15 binów (zbyt duży na realną ewolucję z częstotliwością; inne pulsary miały <3
  biny). `todo.csv`: `dsep=-4.882, sigma=4.6, grade=10` (niespójne z moim wynikiem -6.14 ± 2.46).
  Nic nie zapisane do `separations.csv`.

- **J1257-1027** (n_comp=2): ten sam typ problemu co J1231-4609/J1232-4742, ale z jasnym,
  potwierdzonym wzorcem: **co trzeci pulse jest identycznym duplikatem poprzedniego**
  (pulse3≡pulse2, pulse6≡pulse5, pulse9≡pulse8 — te same wartości fitu co do cyfry).
  To wygląda na systematyczny artefakt generowania p3-foldu (P3 niedzielący się równo przez
  liczbę sub-integracji?), nie przypadek. χ²(G1) = 95-113. `todo.csv`: `dsep=0.784, sigma=7.1,
  grade=10, poszerzenie=1`. Nic nie zapisane do `separations.csv`.
  **TODO dla użytkownika**: rozważyć zbadanie `Data.p3fold_psrdata`/procesu generowania
  `pulsar_*.debase.p3fold_norefine` — dwa pulsary z rzędu (J1231-4609, J1257-1027) pokazują
  wzorce duplikacji, może to systemowy problem dotykający więcej pulsarów z listy.

- **J1302-6350** (n_comp=2): inny typ problemu — bardzo niskie S/N, chaotyczny szum na całej
  szerokości okna (100-860 binów) w KAŻDYM sprawdzonym pulsie (1,2,5,10,15,20 z 20 — rozłożone po
  całym zakresie), brak czystego profilu pulsara widocznego na oko. Wynik końcowy: same NaN
  (prawdopodobnie fit się nie zbiega sensownie / dzielenie przez błąd bliski zeru).
  `todo.csv`: `dsep=-8.868, sigma=6.4, grade=10` (wcześniej uznane za wiarygodne — inna metoda
  lub inne dane wejściowe niż `pulsar_*.debase.p3fold_norefine` użyte tutaj?). Nic nie zapisane
  do `separations.csv`. Wymaga sprawdzenia czy dane wejściowe (`p["bin_st"]`/`bin_end`, plik
  p3fold) są poprawne dla tego pulsara.

- **J1306-6617** (n_comp=2): ten sam typ niestabilności identyfikacji komponentów co
  J1017-5621/J1257-1027. χ²(G1) = 402-415! Potwierdzone wizualnie: w pulsie 2 "Low G2" przeskoczyło
  na drugą stronę ostrego piku (bin~565) podczas gdy "High G2" zostało przy słabym garbie (~456).
  `todo.csv`: `dsep=-2.505, sigma=17.0, grade=10` — mój wynik +5.53 ma nawet inny znak.
  Nic nie zapisane do `separations.csv`.

- **J1328-4921** (n_comp=3): ekstremalna niestabilność etykiet — χ²(G3) = 2165! "Low G3" przy
  bin~501-509, "High G3" przy bin~558-575 (rozjazd ~50-75 binów). `todo.csv`: `dsep=1.543,
  sigma=2.8, grade=10`. Nic nie zapisane do `separations.csv`.

- **J1402-5021** (katalog na dysku: `J1402-5124_16`, stara nazwa — ATNF go przemianował,
  patrz `ppdot` output wcześniej; n_comp=3): niestabilność etykiet we WSZYSTKICH trzech
  komponentach jednocześnie, χ² = 56-104. "Low G3" skacze bin~497↔554 między pulsami.
  `todo.csv`: `dsep=-0.463, sigma=2.3, grade=9`. Nic nie zapisane do `separations.csv`.

- **J1414-6802**: PRAWDOPODOBNIE dryfujące podpulsy (subpulse drift) — ostry pik płynnie
  przesuwa się w fazie pulsa-do-pulsa (bin 504→519 przez 9 pulsów), a nowy komponent
  pojawia się z drugiej strony okna w pulsach 9-10. `analyse_p3folds4` zakłada nieruchome
  komponenty, więc χ² astronomicznie wysokie (1297!) może odzwierciedlać złamanie tego
  założenia astrofizycznego, nie błąd fitu. Zapisane mimo to na życzenie użytkownika (duży
  końcowy błąd ±0.2° i tak to sygnalizuje). `todo.csv`: `dsep=-0.358, sigma=6.6, grade=10`.
  **Uwaga**: wynik w `separations.csv` (dsep=-0.3506) dla tego konkretnego pulsara może nie
  mieć takiej samej interpretacji fizycznej jak dla pulsarów bez dryfu — warto rozważyć
  osobną metodę śledzenia pasma dryfu, jeśli ten pulsar ma znaczenie dla dalszej analizy.

### Zapisane z zastrzeżeniem:
- **J1502-6128**: zachowane tylko pulsy 2,3 (n=2; pulse 1 — zamiana etykiet G1/G2 między Low
  a High, pulse 4 — Low G2 złapało wąski artefakt przy bin 536 zamiast garbu przy 496).
  Wynik −1.746 ± 0.951 vs `todo.csv` −1.095 — zgodne w granicach błędu.
- **J1527-5552** (n_comp=3): zachowane 17 z 20 pulsów (pominięte 1, 19, 20 — zamiana etykiet
  High G2/G3; w 19/20 High G2 przy bin~547 podczas gdy Low G2 przy ~502). Dodatkowo duplikaty
  danych: 2≡3, 4≡5, 17≡18, 19≡20 (ten sam wzorzec co J1231-4609/J1257-1027 — zawyżają n
  i zaniżają błędy). Wynik −0.857 ± 0.694 ma **przeciwny znak** niż `todo.csv` (+0.836),
  χ²(G3)=86 wciąż wysokie. **Do weryfikacji ręcznej.**

## 2026-09-14 — partia 7 pulsarów (J1528-4109 … J1555-3134)

Zmiana w kodzie: `Plot.analyse_p3folds4_agent` dostał guard na `fit.converged` przed
`print_fit_summary` — nieudany fit zostawia `nothing` w polach podsumowania, co rzucało
`MethodError(Float64, (nothing,))` i **przerywało cały przebieg**. Ten sam brak zabezpieczenia
ma oryginalny `analyse_p3folds4` (interaktywny); `analyse_average_offset` guard ma.

### Zapisane (4):
- **J1539-6322** −1.484 ± 0.214 (n=13; todo.csv −1.405) ✓
- **J1543+0929** −2.160 ± 0.281 (n=9; todo.csv −2.160) ✓ idealna zgodność
- **J1548-4821** −4.048 ± 0.739 (n=5; todo.csv −4.047) ✓ idealna zgodność, ale dane bardzo
  zaszumione — w pulsach 3,5 "High G2" szeroki i przesunięty (~522-531) vs "Low G2" (~541-545)
- **J1555-3134** −1.108 ± 0.067 (n=21; todo.csv −1.131) ✓ najczystsze dane w tej partii

### Pominięte (3):
- **J1528-4109** (n_comp=3): wynik NaN. Komponenty rozrzucone poza profilem (Low G3≈473,
  High G2≈566, High G3≈539 przy profilu siedzącym w ~505-525). `todo.csv`: −0.660, sigma 3.4.
- **J1535-4415** (n_comp=2): **wszystkie parzyste pulsy (2,4,…,20) zawierają Inf/NaN** i fit
  się nie zbiega — kolejny regularny artefakt p3-foldu (por. duplikaty w J1231-4609/J1257-1027/
  J1527-5552). Z pozostałych 10 pulsów: −1.84 ± 9.01, χ²=382 — bezużyteczne.
  `todo.csv`: −4.904, sigma 4.2.
- **J1536-3602** (n_comp=3): χ² = 1658/1364/323. High G3 skacze bin 487↔445 (poza profil).
  `todo.csv`: −0.867, sigma 15.5.

### (nieaktualne — lista ukończona 2026-09-16, patrz sekcja końcowa)

## 2026-09-16 — ZBADANY problem z p3-foldami

Pełna diagnoza: **`~/claude/work/p3fold_diagnoza.md`**, dane: `~/claude/work/scan_all.csv`.

Przyczyna: gdy `P3` jest dokładnie całkowite a `p3_ybins > P3`, `pfold -p3fold_norefine`
nadpróbkowuje fazę → `p3_ybins − P3` wierszy to kopie sąsiadów; przy `p3_ybins = 2×P3`
co drugi wiersz jest pusty (stały) → `normalize_per_pulse` robi z niego NaN.

Skan 106 pulsarów: **95 czystych**, 6 z duplikatami, 4 z martwymi wierszami, 1 bez katalogu.
Z dotkniętych tylko **J1527-5552** trafił do `separations.csv` — i jego wpis jest **błędny**
(norefine −0.857 ± 0.694 vs refine +1.175 ± 0.173; todo.csv +0.836, χ² spada z 86 do <4).

**Wcześniejsze notatki o duplikatach wymagają korekty**: to nie jest "systemowy artefakt
dotykający wiele pulsarów" (jak sugerowałem przy J1231-4609/J1257-1027), tylko wąski,
deterministyczny efekt konfiguracji P3/ybins w 10 z 106 przypadków. Pulsary uznane
wcześniej za czyste pozostają czyste.

**Próba naprawy przez `p3_ybins ≤ round(P3)` — NIEUDANA.** Przetestowane na J1527-5552:
artefakt znika, ale χ² rośnie i wynik się pogarsza. Najlepszy jest `refine` na
oryginalnych plikach. Duplikaty były objawem; przyczyną jest niezmierzone P3
(`p3_error` 2–18). Stan przywrócony z backupu, pozostałe 9 pulsarów nietknięte.
Szczegóły i tabela porównawcza: `p3fold_diagnoza.md`.

## 2026-09-16 — partia: 8 zaległych + 6 nowych

### Zapisane bez zastrzeżeń (11)
J1557-4258, J1559-4438, J1609-4616, J1638-3815, J1650-1654, J1651-7642, J1703-4442,
J1707-4417 (przegląd z 2026-09-14, wszystkie pulsy zachowane) oraz J1722-3207
(+0.0001 ± 0.1433), J1727-2739 (−1.5294 ± 0.1806), J1741-0840 (−0.4815 ± 0.0787).

### Zapisane po odrzuceniu wadliwego pulsu (2) — DUŻA zmiana wyniku
Do znalezienia złych pulsów użyty nowy skrypt `scripts/check_order.jl`: sprawdza
kolejność `mu` komponentów w każdym pulsie i zgodność Low/High. Pewniejsze niż
oglądanie wykresów.

- **J1714-1054**: pulse 2 odrzucony (fit High wstawił oba komponenty na lewy pik:
  mu = 495, 492). Wynik: −2.277 ± 1.172 (χ²=303) → **−0.944 ± 0.073 (χ²=0.91)**.
  `todo.csv` −2.241 odpowiada wersji ZE złym pulsem.
- **J1733-3716**: pulse 1 odrzucony (High: mu = 465, 462 — oba na lewym piku).
  Wynik: −12.208 ± 8.463 (χ²=4066) → **−1.827 ± 0.166 (χ²=0.51)**.
  `todo.csv` −12.208 zgadza się co do cyfry z wersją ZE złym pulsem. To ten sam
  pulsar, który w komentarzu w `modules/plot.jl` (skala kolorów `ppdot`) jest opisany
  jako *"the largest offsets reach ~12 deg (J1733-3716, flagged in the input as
  needing a redo)"* — ekstremalny outlier okazał się artefaktem jednego złego fitu.

### Zapisane z zastrzeżeniem (1)
- **J1720-2933**: pulsar z **dryfem podpulsów** — `mu` przesuwa się monotonicznie przez
  serię (507→505→503→501→498→496→495→493, równolegle 533→…→513) i wraca w pulsie 9
  do 509/534, czyli pełen cykl P3 (P3=2.454, ybins=9). W 6 z 9 pulsów numeracja G1/G2
  jest przez to odwrócona. Low/High zgodne w **100%** pulsów, więc Δsep = −0.0372 ± 0.0745
  jest poprawne (todo.csv −0.070). Zapisane na decyzję użytkownika.
  Ten sam przypadek co [J1414-6802] — oba wymagałyby metody śledzącej pasmo dryfu.

  **KOREKTA (2026-09-16).** Napisałem wyżej, że kolumny lon 1/lon 2/sep „mieszają różne
  fazy dryfu i nie należy ich używać". To był **błędny wniosek**. `GaussianFit.component_offsets`
  (`modules/GaussianFit.jl:258-259`) **sortuje komponenty po `mu`** przed sparowaniem
  Low↔High, więc G1 zawsze oznacza komponent lewy — kolejność, w jakiej fit je znalazł,
  nie ma znaczenia. Realny problem J1720-2933 to dryf **pozycji** (lewy komponent wędruje
  w zakresie 493-509 binów), co podbija χ² longitude, ale średnia pozycja pozostaje
  sensowna. Kolumny lon/sep są użyteczne, tylko z dużym rozrzutem.
  Potwierdzenie: wszystkie 71 zapisanych wierszy ma **dodatnią** separację, co przy
  odwróconej numeracji byłoby niemożliwe.

## 2026-09-16 — partia 11 pulsarów (J1746+2245 … J1822+0705)

Nowe, twarde kryterium w `scripts/check_order.jl`: dla każdego pulsu sortuje komponenty
po pozycji i liczy `max |mu_low − mu_high|`. Powyżej **10 binów** to zawsze zły fit
(realne offsety w tej analizie to < 3 biny przy nbin=1024). Zastąpiło oglądanie wykresów —
szybsze i wychwytuje dokładnie ten tryb awarii, który psuł wyniki.

### Zapisane bez odrzuceń (6)
J1746+2245 (+0.391 ± 0.421), J1801-2920 (−0.419 ± 0.053), J1810-5338 (−0.665 ± 0.083),
J1811-4930 (−0.786 ± 0.175), J1822+0705 (−0.939 ± 0.169) — wszystkie zgodne z `todo.csv`.
- **J1806-1154** (−0.446 ± 0.178, n=20): wynik zgadza się z `todo.csv` co do cyfry, ale
  ten pulsar ma **duplikaty w p3-foldzie** (P3=12 całkowite, ybins=20 → 8 duplikatów).
  Efektywnie ~12 niezależnych faz, więc **błąd jest zaniżony** o czynnik ~√(20/12) ≈ 1.3.
  `todo.csv` ma ten sam problem, stąd zgodność.

### Zapisane po odrzuceniu złych pulsów (4) — wszystkie miały PRZECIWNY ZNAK przed poprawką
- **J1757-2421**: odrzucone [1,8,9,10,12,22,23] → n=16.
  −2.185 ± 1.508 (χ²=90) → **+1.637 ± 0.253** (χ²=0.9). `todo.csv` +1.700 ✓
- **J1803-3329**: odrzucone [12,13,14,15,16,17,19] → n=12.
  +1.006 ± 1.099 (χ²=623) → **−0.441 ± 0.204** (χ²=13). `todo.csv` −0.466 ✓
- **J1807-0847**: odrzucone [3,4,5,6] → n=3.
  −1.232 ± 1.579 (χ²=344) → **+0.431 ± 0.122** (χ²=0.1). `todo.csv` +0.421 ✓
  UWAGA: zostały tylko 3 pulsy.
- **J1808-3249**: odrzucone [5,6,8,9,10,12,13,15,16,17] → n=10.
  Wynik NaN → **−0.058 ± 0.145** (χ²=1-7). `todo.csv` −0.349 (rozbieżność ~2σ).

### Pominięte (1)
- **J1819+1305** (n_comp=3): po odrzuceniu 14 z 20 pulsów zostaje n=6, a χ² longitude
  wciąż 110–583. `todo.csv` ma `dsep=0.007, sigma=0.0`, czyli oryginalna analiza też
  uznała to za nieistotne. Nic nie zapisane.

### Znaleziona luka w `_offset_summary` (nie naprawiona)
J1808-3249 i J1819+1305 dawały NaN, bo **pojedynczy** pulse miał `offset_err = NaN`
(osobliwa macierz kowariancji przy skrajnie złym dopasowaniu). `_offset_summary`
filtruje tylko przypadek, gdy WSZYSTKIE błędy są zerowe
(`!all(offset_data[c].err .== 0.0)`), więc jeden NaN zatruwa całą średnią ważoną
przez `w = 1/err²`. Propozycja: odfiltrować niefinite `err` przed liczeniem wag.

## 2026-09-16 — partia 10 pulsarów (J1823+0550 … J1847-0402)

### Poprawione kryterium w `scripts/check_order.jl`: ABSOLUTNE → ADAPTACYJNE
Próg stały (10 binów) okazał się zły: J1843-0459 miał odrzucone **wszystkie 5** pulsów,
mimo że wynik (+2.038 vs todo +1.978) i χ² (0.1-0.2) były znakomite — jego realny offset
to 2-4°, czyli 6-11 binów. Próg absolutny mylił sygnał z błędem.
Nowe kryterium odrzuca pulsy ODSTAJĄCE od reszty danego pulsara:
`d > median(d) + 4 * 1.4826 * MAD(d)`, plus bezpiecznik 60 binów.
Efekt: J1843-0459 5/5 zachowanych, J1834-0426 21/21 (było 11 odrzuceń), J1842-0359 12/12.

### Zapisane (10)
- bez odrzuceń: J1834-0426 (+1.313), J1834-1202 (−0.905), J1840-0809 (+0.441),
  J1842-0359 (−3.788), J1843-0459 (+2.038), J1847-0402 (+1.074) — zgodne z `todo.csv`
- **J1823+0550** (+0.844 ± 0.124, n=3/6): 3 wiersze p3-foldu martwe (P3=3, ybins=6).
  Wynik mimo to idealnie zgodny z todo.csv (+0.844).
- **J1827-0750** (−1.384 ± 0.762, n=3/10): 5 wierszy martwych (P3=5, ybins=10) + 2 odrzucone.
  Wynik zgodny z todo.csv (−1.383), ale jeden fit zwrócił `A = 0.0000` i NaN w
  niepewnościach — statystyka bardzo słaba.
- **J1824-0127** (+0.960 ± 0.401) — `todo.csv` −0.310, rozbieżność ~3.2σ
- **J1825-0935** (+0.513 ± 0.304) — `todo.csv` −0.273, rozbieżność ~2.6σ

  Dla obu ostatnich: kryterium jakości NIE wskazuje żadnych wadliwych pulsów, numeracja
  spójna, fity poprawne. Komponenty są jednak słabo rozdzielone (separacja 2.98° i 2.61°,
  czyli ~8 binów), co przy takich odległościach czyni dopasowanie wrażliwym. Precedens
  J1733-3716 i J1714-1054 pokazuje, że rozbieżność z `todo.csv` bywa winą starej wartości,
  nie nowej — ale tutaj przyczyna nie jest ustalona. **Do weryfikacji ręcznej.**

## 2026-09-16 — partia 11 pulsarów (J1850+0026 … J1919+0134)

Workflow zautomatyzowany: `scripts/batch_auto.jl` liczy kryterium jakości i od razu
podgląd wyniku (`psr=""`) w jednym przebiegu Julii. Wykresy przeglądowe i tak powstają,
ale decyzja o odrzuceniu pulsów jest podejmowana na podstawie liczb, nie oglądania.

### Zapisane — wszystkie 11
Dziesięć zgadza się z `todo.csv` znakomicie (często co do trzeciego miejsca), χ² niskie:
J1850+0026 (−0.792, odrzucony pulse 1), J1852-0635 (−3.482, odrzucony 8),
J1900-2600 (+0.120), J1901+0716 (+0.013, odrzucone 2 i 14), J1901-0906 (−0.286),
J1909+1102 (−0.197), J1910+0728 (−0.885), J1912+2104 (+0.215), J1914+0219 (−0.172),
J1914+1122 (+0.040).

- **J1919+0134** (−0.599 ± 0.060, n=12/13, odrzucony pulse 11): `todo.csv` ma −0.972.
  Rozbieżność ~2σ względem błędu starej wartości. Bez odrzucenia pulse 11 wychodzi
  −0.485 ± 0.492 przy **χ²(G2) = 350**, czyli ten pulse na pewno jest zły — ale żaden
  z wariantów nie daje −0.972, więc oryginalna analiza odrzucała inny zestaw pulsów
  (interaktywne 's'). Zostawiam wersję z χ²=0.8.

## 2026-09-16 — ZAKOŃCZENIE listy (ostatnie 9)

### Zapisane (6)
J1921+2003 (−1.267 ± 0.229), J1922+1733 (−1.988 ± 0.624), J1933+1304 (−0.142 ± 0.032),
J2037+1942 (+0.015 ± 0.115), J2046+1540 (−0.113 ± 0.055, odrzucony pulse 13),
J2053-7200 (−1.767 ± 0.199) — wszystkie zgodne z `todo.csv`, bez odrzuceń poza J2046+1540.

### Odrzucone (3) — ten sam, nowy tryb awarii
- **J1932+1059** (n_comp=3): mediana |Δmu| = 17.4 bina, χ²(G2)=109. `todo.csv` +0.525.
- **J1402-5021** (katalog `J1402-5124_16`, n_comp=3): mediana |Δmu| = 16.1 bina,
  χ² = 56-104. `todo.csv` −0.463. (Odrzucony już wcześniej, potwierdzone ponownie.)
- **J1819+1305** (n_comp=3): mediana |Δmu| = 23.4 bina, wynik NaN. `todo.csv` +0.007.

**Ograniczenie kryterium adaptacyjnego.** W tych trzech przypadkach rozjazd Low↔High jest
**systematyczny** — dotyczy wszystkich pulsów, nie pojedynczych. Kryterium
`d > median(d) + 4·MAD` szuka outlierów względem reszty danego pulsara, więc gdy złe są
wszystkie pulsy, nie wskazuje niczego do odrzucenia. Sygnałem jest wtedy sama **mediana**
|Δmu|: zdrowe pulsary mają 0.2-4.9 bina, te trzy mają 16-23 biny (5.7-8.2°), podczas gdy
`todo.csv` podaje dla nich offsety < 1°. Warto dodać do skryptu ostrzeżenie, gdy
mediana |Δmu| przekracza np. 10 binów.

### Stan danych po tym epizodzie
- `~/output/claude/J1527-5552_16/` — przywrócony do oryginału (params.json + 4 p3foldy).
- Backup 10 pulsarów zachowany w `~/claude/work/backup_p3fold/` (50 plików, 16 MB).
- `input/separations.csv` — bez zmian; wpis J1527-5552 **nadal jest tym błędnym**
  (−0.857 z norefine). Do poprawienia, gdy zapadnie decyzja o metodzie.

---

# PODSUMOWANIE listy separations_todo.csv (2026-09-16; poprawki z 17.09 niżej)

**106 pulsarów: 91 zapisanych w `input/separations.csv`, 15 odrzuconych.**

### Odrzucone (15) i przyczyny
- **Systematyczny rozjazd Low↔High we wszystkich pulsach** (mediana |Δmu| 16-23 biny przy
  realnych offsetach < 1°): J1932+1059, J1402-5021, J1819+1305, J1328-4921, J1536-3602.
- **Niestabilne etykietowanie / bardzo wysokie χ²**: J1231-4609, J1232-4742, J1257-1027,
  J1306-6617, J1625-4048, J1627-5936.
- **Martwe wiersze p3-foldu lub NaN**: J1528-4109, J1535-4415, J1822-4209.
- **Zbyt niskie S/N**: J1302-6350.

### Wyniki wymagające weryfikacji ręcznej (zapisane, ale oznaczone)
- **J1527-5552** — wpis pochodzi z `norefine` i jest **błędny** (artefakt duplikatów
  p3-foldu); `refine` daje +1.175 ± 0.173 przy χ² < 4 zamiast −0.857 ± 0.694.
- **J1824-0127** (+0.960 vs todo −0.310) i **J1825-0935** (+0.513 vs todo −0.273) —
  rozbieżność 2.6-3.2σ bez wykrywalnej przyczyny; komponenty słabo rozdzielone (~8 binów).
- **J1919+0134** (−0.599 vs todo −0.972) — moja wersja ma χ²=0.8, stara nie odtwarza się
  żadnym zestawem odrzuceń.
- **J1414-6802**, **J1720-2933** — pulsary z dryfem podpulsów; Δsep poprawne, ale χ²
  longitude wysokie z natury zjawiska.
- **J1806-1154** — błąd zaniżony (~√(20/12)) przez duplikaty w p3-foldzie.
- **J1823+0550** (n=3/6), **J1827-0750** (n=3/10), **J1807-0847** (n=3/7) — bardzo słaba
  statystyka po odrzuceniach.

### Wartości w todo.csv, które okazały się błędne
- **J1733-3716**: todo −12.208 (największy outlier na diagramie P-Ṗ) to artefakt jednego
  złego fitu. Poprawnie: **−1.827 ± 0.166** (χ²=0.5). Zgodne z komentarzem w kodzie
  (`modules/plot.jl`: "flagged in the input as needing a redo").
- **J1714-1054**: todo −2.241 → **−0.944 ± 0.073** (χ² z 303 na 0.91).
- **J1757-2421, J1803-3329, J1807-0847**: wszystkie miały w todo przeciwny znak niż
  poprawny wynik po odrzuceniu wadliwych fitów.

### Narzędzia powstałe przy okazji (`~/claude/work/scripts/`)
- `check_order.jl` — kryterium jakości per pulse (adaptacyjne, `median + 4·MAD`).
- `batch_auto.jl` — kryterium + podgląd wyniku w jednym przebiegu Julii.
- `scan_all.jl`, `diag_*.jl` — diagnostyka p3-foldów.

### Otwarte zadania dla użytkownika
1. Poprawić wpis **J1527-5552** (decyzja: `norefine` czy `refine` jako metoda).
2. ~~Luka w `_offset_summary`~~ — **NAPRAWIONE 2026-09-17**, patrz sekcja niżej.
3. Rozważyć ostrzeżenie w kodzie, gdy `p3_ybins > round(P3)` przy całkowitym P3
   (10 pulsarów dotkniętych, patrz `p3fold_diagnoza.md`).
4. Zweryfikować 4 wyniki rozbieżne z todo.csv (lista wyżej).


---

## 2026-09-17 — NAPRAWIONA luka w `_offset_summary` (modules/plot.jl)

### Problem
Wagi liczone są jako `w = 1/err^2`. Pojedynczy pomiar z `err = NaN` (osobliwa macierz
kowariancji w `GaussianFit`) albo `err = 0` (waga `Inf`) zamieniał **każdą** sumę poniżej
w NaN, bo w arytmetyce zmiennoprzecinkowej dowolne działanie z NaN daje NaN.
Efekt: jeden zły pulse z dwudziestu unieważniał całego pulsara.
Stary filtr `!all(err .== 0.0)` tego nie łapał — `NaN == 0.0` jest `false`,
więc komponent przechodził dalej.

### Poprawka
Filtr per komponent zamiast per pulsar:
```julia
usable(c) = findall(e -> isfinite(e) && e > 0, offset_data[c].err)
filter!(c -> !isempty(usable(c)), comps)
```
a w pętli `off/lon/err` brane są tylko z indeksów `usable(c)`; `n` to liczba faktycznie
użytych pomiarów. Gdy coś odpadnie, wypisywany jest komunikat
`G<c>: dropped N of M measurements (err = 0, NaN or Inf)`.

### Weryfikacja — brak regresji
Przeliczone 5 kontrolnych pulsarów o różnym `n_comp` i `n` (J0134-2937, J1555-3134,
J1834-0426 z n_comp=4, J1919+1745, J1909+1102): wszystkie wartości **identyczne co do
ostatniej cyfry** z zapisanymi, zero komunikatów o odrzuceniach.
Dodatkowo: żaden z 91 zapisanych wierszy nie zawiera NaN/Inf — gdyby `err` miał NaN
lub 0, wynik byłby NaN i nie zostałby zapisany. Poprawka nie zmienia żadnego
istniejącego wyniku. Backup sprzed zmiany: `~/claude/work/separations_przed_poprawka.csv`.

### Wpływ na dwa pulsary, które dawały NaN
- **J1808-3249** — zapisany wynik (−0.058 ± 0.145, n=10) odtwarza się identycznie:
  wadliwy pulse był już wykluczony kryterium jakości, więc poprawka nic nie zmienia.
- **J1819+1305** — po poprawce daje liczbę zamiast NaN (+1.756 ± 1.263 bez odrzuceń),
  ale χ² = 76-600. Ma **drugi, niezależny problem**: systematyczny rozjazd Low↔High
  (mediana |Δmu| = 23.4 bina = 8.2°, przy `todo.csv` +0.007). NaN był objawem
  towarzyszącym, nie przyczyną. **Pozostaje odrzucony.**

Wniosek: poprawka usuwa realną pułapkę na przyszłość (i dotyczy też interaktywnego
`analyse_p3folds4` oraz `analyse_average_offset`, bo `_offset_summary` jest wspólne),
ale nie odblokowuje żadnego z 15 odrzuconych pulsarów.

---

## 2026-09-17 — nowy diagram P-Ṗ: frakcyjne zwężenie ΔW/W

### Odzyskanie poprawionych Δsep
`_offset_summary` liczyło Δsep i ΔW/W, ale zapisywało do `separations.csv` wyłącznie `sep`
(mid-band). Poprawione Δsep z całej kampanii istniały więc tylko w logach.
Wyciągnięte skryptowo z `~/claude/work/logs/*.log`: dla każdego bloku
`separation G..-G.. =` → `Δ separation =` → `<psr> written to` pobrana para (sep, Δsep),
przypisana do pulsara z linii `written to`. **Walidacja: `sep` z logu zgadza się z wierszem
w `separations.csv` co do 1e-4° dla 91/91** — czyli sparowany został dokładnie ten przebieg,
który zapisał wiersz. Kopia wyciągu: `~/claude/work/separations_final.csv`.

21 wartości różni się od `offsets.csv` o >0.05°, w tym potwierdzone wcześniej artefakty:
J1733-3716 −12.208 → −1.827, J1714-1054 −2.241 → −0.944, J1527-5552 +0.836 → −0.857,
zmiany znaku J1824-0127 (−0.310 → +0.960) i J1825-0935 (−0.273 → +0.513).
`offsets.csv` NIE był ruszany (per-komponentowych offsetów nie da się odtworzyć z logów).

### Zmiany w repo (gałąź claude)
- `input/separations.csv` — dwie nowe kolumny `dsep,dsep err`, wypełnione dla 91 wierszy.
  Backup sprzed zmiany: `~/claude/work/separations_przed_dsep.csv`.
- `_offset_summary` — writer zapisuje teraz także `dsep`/`dsep err`, więc na przyszłość
  ta wielkość nie ginie.
- `_read_separations(file; grades=offsets.csv)` — nowy reader, zwraca `dsep/sep` (bezwymiarowe)
  w tym samym NamedTuple co `_read_offsets`; oceny jakości dołączane po nazwie z `offsets.csv`
  (w `separations.csv` ich nie ma).
- `_ppdot` — `offsets` przyjmuje teraz także gotowy `Dict`, nie tylko ścieżkę; nowe kwargi
  `offset_norm` (`:symlog` / `:linear`), `offset_label`, `offset_legend`; romb w legendzie
  tylko gdy faktycznie jest pulsar 1-komponentowy; na skali liniowej colorbar dostaje tick
  w zerze, a drabinka ticków startuje od zera (inaczej 0.01 nachodziło na 0).
- `Plot.ppdot_separations(outdir)` — nowa funkcja. Domyślnie `highlight=nothing`:
  wszystkie 91 pulsarów jest w `pulsars.txt`, więc magentowa warstwa 533 punktów tylko
  zasłaniała symbole. Żeby ją przywrócić: `highlight=".../input/pulsars.txt"`.

### Wynik
`~/output/claude/ppdot_separations.png` / `.pdf`.
90 z 91 narysowanych (J1714-1054 ma ocenę < 6 w `offsets.csv` — ocena dotyczy starego,
błędnego fitu, do rewizji). Skala koloru liniowa ±0.2, 4 pomiary poza zakresem.
Statystyka ΔW/W: mediana −0.045, zakres −0.309…+0.322, 66/91 ujemnych (zwykły RFM).
Skrajne: J1824-0127 +0.322 ± 0.136, J1430-6623 −0.309 ± 0.031, J1922+1733 −0.266 ± 0.087.

### Brak regresji
`ppdot_offsets` i `ppdot` przeliczone: 159/161 pulsarów, zakres ±2.2, 15 nasyconych —
identycznie jak przed zmianą (wartości zgodne z komentarzem w kodzie).

### Otwarte
- Ocena jakości (`ocena` w `offsets.csv`) opisuje stare fity; J1714-1054 wypada przez to
  z wykresu, choć jego poprawiony wynik ma χ² = 0.91. Warto przenieść oceny do
  `separations.csv` albo zrewidować je dla 21 poprawionych pulsarów.

### Analiza zależności ΔW/W od parametrów pulsara (2026-09-17)

Próbka jak na wykresie: 90 pulsarów (ocena ≥ 6), P i Ṗ z `psrcat.db` (replika `read_psrcat`
w pythonie, zwalidowana liczbą rekordów 2930). Tabela: `~/claude/work/ppdot_sep_table.json`,
skrypt wykresu: `scripts/frac_vs_params.py` → `~/output/claude/frac_vs_params.png`.

**Wagi 1/σ² są tu nieużywalne.** χ²/dof wokół średniej ważonej = 64 → rozrzut próbki (σ = 0.096)
jest 8× większy niż mediana błędu pomiarowego (0.013). Wszystkie fity ważone dają absurdalne
istotności (35-83σ przy r ≈ 0). Używane: Spearman + nachylenie MNK z błędem bootstrapowym.

| parametr | ρ_S | p | nachylenie ΔW/W na dekadę |
|---|---|---|---|
| log Ėdot | +0.25 | 0.018 | +0.016 ± 0.010 |
| log τ_c | −0.25 | 0.019 | −0.029 ± 0.012 |
| log Ṗ | +0.18 | 0.094 | +0.024 ± 0.011 |
| log B_d | +0.08 | 0.44 | +0.031 ± 0.019 |
| log P | −0.10 | 0.33 | −0.007 ± 0.029 |

Ėdot/τ_c/Ṗ to ten sam trend (te wielkości są współliniowe). Brak zależności od P i B_d.
Odporny na podpróbki: grade≥8 (ρ = ±0.28, p ≈ 0.01), bez 4 podejrzanych wyników (p ≈ 0.02),
|f|/err ≥ 3 (p = 0.02-0.08), usunięcie najstarszego J1548-4821 (ρ = −0.225, p = 0.034).
Po korekcie na 7 testowanych parametrów — **poszlaka, nie wynik**. Bonferroni ×7 daje 0.124,
ale jest zbyt ostry: te parametry to funkcje P i Ṗ, dwie pierwsze składowe główne biorą 88%
wariancji (N_eff ≈ 5.2). Dokładny rachunek to permutacja ΔW/W ze statystyką max |ρ| po
wszystkich siedmiu, 200 000 prób: **p_globalne = 0.075**. Z drugiej strony Ė wskazał już
poprzedni raport, więc jako hipoteza postawiona z góry test nie wymaga korekty i zostaje
p = 0.018 — tyle że wszystkie 90 pulsarów nowej próbki było w starej (105 wielokomponentowych,
pokrycie 100%), więc to nie jest niezależne potwierdzenie. Uczciwy przedział: 0.018–0.075.
Moc: przy prawdziwym ρ = 0.25 i N = 90 to 66% na p < 0.05 i 27% na 3σ; na rozstrzygnięcie
trzeba 165 i 285 pulsarów.

Wielkość efektu: mediana ΔW/W −0.062 (dolny kwartyl Ėdot) vs −0.035 (górny), różnica 2.7 p.p.,
czyli ~28% rozrzutu próbki. Poszerzenia (ΔW/W > 0) 13% vs 30% (Fisher p = 0.28).

**Kontrola normalizacji:** |Δsep| w stopniach silnie koreluje z szerokością profilu
(ρ = +0.62), |ΔW/W| już nie (ρ = −0.14, p = 0.21) — dzielenie przez `sep` faktycznie usuwa
zależność od szerokości, czyli frakcja jest właściwą wielkością.

**Wniosek:** dominującego efektu nie ma. Rozrzut jest 8× większy od błędów, więc rządzi nim
coś spoza płaszczyzny P-Ṗ — najbardziej naturalnie geometria (α, β). Następny krok:
skorelować ΔW/W z geometrią RVM (katalogi `*_rvm` w `~/output/claude/`).

## 2026-09-17 — inwentaryzacja ~/output/claude: co da się odzyskać

Użytkownik potwierdził: pozostałe pulsary z próbki Song et al. (2023) BYŁY próbowane,
ale wyniki nie zostały zapisane. Stąd brak śladu w `offsets.csv` (162 wiersze przy 533 w
`pulsars.txt`).

### Stan danych
534 katalogi `<psr>_16` = cała próbka macierzysta. Zawartość (odczyt, nic nie ruszane):

| plik | zmierzone (91) | próba b/w (71) | reszta (372) |
|---|---|---|---|
| `pulsar_low.debase.p3fold_refine` | 91 | 71 | 356 |
| `pulsar_high.debase.p3fold_refine` | 91 | 71 | 348 |
| `params.json` | 91 | 71 | 371 |

**510 z 533 ma komplet low+high p3foldu** (refine i norefine). `params.json` trzyma
`p3`, `p3_error`, `p3_ybins`, `bin_st`, `bin_end`, `nbin`, `nsubint`. P3-foldy są ASCII
(nagłówek + wiersze `sub chan bin wartość`), czytelne bez PSRCHIVE.

Wniosek kosztorysowy: **etap kosztowny (dedyspersja, debasing, LRFS/2DFS, p3-fold) jest
policzony dla całej próbki.** Ponowny pomiar offsetów = wyłącznie dopasowanie gaussów na
gotowych p3-foldach, bez dotykania `~/data`.

### Wąskie gardło: n_comp
`params.json` nie zawiera liczby składowych, a `analyse_p3folds_16_new` wymaga jej jako
argumentu. Test wykonalności (`scripts/ncomp_probe.py`, wynik `ncomp_probe.csv`):
prosty licznik pików (find_peaks, próg 5σ nad rms poza oknem, prominencja 3σ) na profilu
scalonym z p3-foldu, walidowany na 91 pulsarach o znanym `n_comp`:

- **zgodność dokładna 63%** (57/91) — za mało na przebieg wsadowy
- główny tryb błędu: **19 z 75 dwuskładnikowych uznane za jednoskładnikowe** — składowe
  zlane (rzędu 8 binów), piku nie ma, ale dopasowanie dwóch gaussów sobie radzi
- poprawnie „≥2 składowe": 77%

Licznik pików to zły przyrząd. Właściwy test: dopasować n = 1, 2, 3 i wybrać po BIC —
do zrobienia i zwalidowania na tych samych 91.

### Rozbicie 443 pulsarów bez wyniku ΔW/W
| status automatu | liczba |
|---|---|
| ok (piki zgodne low/high) | 260 |
| rozjazd liczby pików low vs high | 131 |
| za niski S/N | 27 |
| brak p3-foldu | 14 |
| katalog szczątkowy (sam `params.json`) | 11 |

Czyli **391 ma sprawne p3-foldy** i nadaje się do ponownej próby.

Uboczne: flaga „rozjazd" trafia 9 z 15 pulsarów odrzuconych w tej kampanii, przy 33 ze 91
zapisanych — sygnał jest, ale nie rozdziela czysto, więc nie nadaje się na samodzielne
kryterium odrzucania.

### Selektor n_comp po BIC — WYNIK NEGATYWNY (2026-09-17)

`scripts/ncomp_bic.py` — dla każdego katalogu czyta p3-fold low i high, sumuje po binach
P3, dopasowuje n = 1..4 gaussów (wielostart: piki, równy podział, podział przesunięty)
i zapisuje RSS(n). Wybór n robiony offline, żeby dało się skanować kryterium bez
ponownego dopasowywania. Dane: `ncomp_rss_known.csv` (105), `ncomp_rss_all.csv` (534).

**BIC z prawdziwym szumem nie działa w ogóle.** Profil scalony ma ogromne S/N (suma
p3_ybins wierszy po ~1000 pulsach), więc rms off-pulse jest o rzędy mniejszy niż realna,
niegaussowska struktura profilu — kryterium zawsze wybiera n = 4. Użyta postać bezskalowa,
z wariancją estymowaną z residuów: `BIC = N ln(RSS/N) + lam·k·ln N`.

Trafność na 105 pulsarach o znanym n_comp (`separations_todo.csv`):

| metoda | trafność |
|---|---|
| **stała n = 2 (linia bazowa)** | **79%** |
| licznik pików (`ncomp_probe.py`) | 61% |
| BIC, lam = 8, min(low, high) | 58% |
| BIC przycięte do [2,3] | 68% |
| zgoda piki+BIC, inaczej 2 — **w próbie** | 85% |
| to samo, **uczciwa walidacja 2-fold** | **79%** |

85% było artefaktem strojenia na tych samych danych. Po podziale na połowy (lam wybierane
na treningu, ocena na teście) zostaje 79%, czyli dokładnie linia bazowa. McNemar konsensus
vs stała 2: poprawia 9, psuje 5, **p = 0.42**. Żadnej poprawy.

Rozkład prawdy: 83× n=2, 21× n=3, 1× n=4. Błędy reguły konsensusu: (3→2) 12×, (2→3) 5×,
(4→3) 1× — dominuje niedoszacowanie.

**Dlaczego to nie mogło zadziałać.** `n_comp` to liczba składowych, które da się ŚLEDZIĆ
w p3-foldzie, a profil scalony wyrzuca dokładnie tę informację (dryf), która ją definiuje.
Do tego realne profile mają skrzydła i asymetrie, które kryterium chce opisać dodatkowymi
gaussami, a składowe zlane (rzędu 8 binów) nie dają osobnego minimum. Negatywny wynik
dotyczy więc profilu scalonego, nie zagadnienia w ogóle — cechą fizycznie właściwą jest
struktura 2-D p3-foldu (liczba ścieżek dryfu w płaszczyźnie faza–P3). Niesprawdzone.

**Wniosek operacyjny: selektor jest niepotrzebny.** Puścić cały wsad z `n_comp=2`:
79% trafień od razu, a błędy są wykrywalne kryteriami jakości z poprzedniej kampanii
(χ², mediana |Δμ|). Kolejka do ręcznego przeglądu ~82 z 391 zamiast 391.

Flaga „sprawdź, czy nie 3-składnikowy" (zgoda piki+BIC na n=3) w walidacji 2-fold:
precyzja 53%, czułość 45%. Za słaba na decyzję, wystarczająca do posortowania kolejki.

### Lista wsadowa gotowa: `~/claude/work/batch_todo.csv` (2026-09-17)

Pełny przebieg `ncomp_bic.py rss all` na 534 katalogach (21 min) → `ncomp_rss_all.csv`.
Status: 471 ok, 38 za niski S/N, 14 bez p3-foldu, 11 katalogów szczątkowych.

Po odjęciu 91 pulsarów, które już mają ΔW/W: **380 pulsarów do przeliczenia**
(70 próbowanych bez wyniku + 310 nigdy nietkniętych). Kolumny: `psr`, `n_comp` (wszędzie 2,
patrz wynik negatywny wyżej), `flaga3` (11 pulsarów do sprawdzenia pod kątem n=3),
`snr_low`, `snr_high`, `grupa`. Posortowane malejąco po S/N — najmocniejsze profile idą
pierwsze, więc przerwany wsad zostawia najlepszy możliwy podzbiór.

S/N profilu scalonego (low): mediana 63, kwartyle [30, 151], min 6. Progu S/N nie nakładam:
przy rozrzucie populacyjnym 0.096 pomiar z błędem 0.08 nadal niesie 59% wagi, a słabe
przypadki i tak odsieją kryteria jakości (χ², mediana |Δμ|).

Oczekiwany plon, ostrożnie: 380 × 0.79 (trafiony n_comp) × 0.7–0.86 (skuteczność
dopasowania; górny kraniec z poprzedniej kampanii, która była WYSELEKCJONOWANA) ≈ 210–260
nowych pomiarów. Razem z obecnymi 90 daje to N ≈ 300–350, czyli okolice sufitu 349
i moc ~90% na 3σ przy rho = 0.25.

### Test wsadowy na 10 pulsarach (2026-09-17)

`scripts/batch_run.jl` — czyta `batch_todo.csv`, dla każdego pulsara liczy kryterium
jakości per puls (mediana + 4·MAD na max|Δμ|), stosuje maskę `keep` i woła
**`Data.Plot.analyse_p3folds4_agent`** (wariant nieinteraktywny; zwykły `analyse_p3folds4`
czeka na klawisz i zawiesiłby `psrx`). Zapis kierowany przez `separations=` do pliku
roboczego, obrazki przez `outdir=` do `/home/psr/work/review/<psr>/` — repo i katalogi
danych nietknięte.

Trzy bramki: `MED_MAX = 10` binów (systematyczny rozjazd Low↔High we wszystkich pulsach —
otwarte zadanie z 16.09, teraz zaimplementowane), `MIN_PULSES = 5`, oraz kontrola
separacji po fakcie.

**Kontrola separacji dodana po pierwszym przebiegu testowym.** J0837-4135 dał
W = 0.100° ± 0.103° — `n_comp=2` narzucone pulsarowi jednoskładnikowemu dopasowało dwa
gaussy w to samo miejsce, a `_offset_summary` policzył z tego ΔW/W = +0.095 ± 0.306.
Progi `MIN_SEP = 1.0°` i `MIN_SEP_SN = 5`: w 91 zweryfikowanych pulsarach minimum to
W = 2.61° przy 8.8σ, więc nie odrzucą niczego, co już przeszło weryfikację. Wiersz
odrzuconego pulsara jest usuwany z pliku wyjściowego (`drop_row`).

Wynik na pierwszej dziesiątce (posortowanej po S/N):

| status | liczba | pulsary |
|---|---|---|
| ok | 6 | J0820-1350, J1001-5507, J1243-6423, J1534-5334, J1903+0135, J1921+2153 |
| rozjazd systematyczny | 2 | J1645-0317 (mediana 13.8 bina), J1327-6222 (12.9) |
| składowe zlane | 1 | J0837-4135 |
| za mało pulsów | 1 | J2048-1616 (3 z 6 po odrzuceniach) |

**Plon 60%** — powyżej dolnego krańca prognozy (0.79 × 0.7 ≈ 55%), ale to próbka 10
pulsarów o najwyższym S/N, więc raczej górne oszacowanie niż typowe. Czas: ~5 min na
10 pulsarów w jednym procesie Julii, czyli ~3 h na całe 380.

Wyniki: `~/claude/work/separations_batch.csv` (6 wierszy), logi `logs/batch_test10*.log`,
obrazki przeglądowe `~/claude/work/review/` (4.7 MB dla 7 pulsarów → ~250 MB dla 380).

Do rozważenia przed pełnym przebiegiem: `bad_pulses` i `analyse_p3folds4_agent` dopasowują
gaussy dwukrotnie (raz na kryterium, raz na wynik) — pełny przebieg da się skrócić o połowę,
jeśli kryterium będzie liczone wewnątrz funkcji agentowej.

### Optymalizacja batch_run.jl — i korekta kosztorysu (2026-09-17)

**Domniemane wąskie gardło nie istniało.** Podejrzenie padło na podwójne dopasowanie gaussów
(raz w kryterium, raz w funkcji agentowej). Pomiar (`scripts/timing_probe2.jl`) pokazał coś
innego:

| etap | czas na pulsara |
|---|---|
| `include` + `using` | 9.7 s (raz na proces) |
| wczytanie p3-foldów | 0.02 s |
| kryterium jakości (wszystkie fity) | **0.02–0.06 s** |
| `analyse_p3folds4_agent` | **1.8–6.8 s** |

Po rozgrzaniu JIT dopasowania są praktycznie darmowe — cały koszt to matplotlib, ~0.3 s na
puls. Usunięcie podwójnego fitu dałoby zero.

**Co zostało zrobione:**
- `Plot.analyse_p3folds4_agent` dostał kwarg `plots=true`; `plots=false` pomija rysowanie
  (blok `if plots` wokół figury). Zmiana w repo, 10 linii.
- `batch_run.jl`: jedno wczytanie zamiast dwóch, przebieg z `plots=false`, a wykresy
  dorysowywane **tylko dla pulsarów idących do przeglądu** — gdy kryterium odrzuciło jakiś
  puls, gdy max χ²/dof offsetu > `CHI_REVIEW = 5`, gdy wynik odpadł na kontroli separacji,
  albo gdy pulsar ma `flaga3` w `batch_todo.csv`. Nowy status `ok_do_przegladu`.
- 4. argument `plots` wymusza rysowanie dla wszystkich (regeneracja obrazków, benchmark).

**Walidacja:** `separations_batch.csv` po optymalizacji jest **identyczny co do bajtu** z
wersją sprzed niej, w obu trybach rysowania.

**Zysk mniejszy, niż się wydawało:** 41 s vs 49 s na dziesiątce (powtarzalnie), czyli ~16%.
Obrazków 2.9 MB zamiast 4.7 MB, 3 katalogi do obejrzenia zamiast 7 — i to jest realna
korzyść: kolejka przeglądu zawiera tylko przypadki, które faktycznie tego wymagają.

**KOREKTA:** wcześniejszy wpis mówił „~5 min na 10 pulsarów, ~3 h na całe 380". To było
oszacowanie, którego nie zmierzyłem, i było zawyżone. Zmierzone: 41 s na 10 pulsarów, z czego
~15 s to jednorazowy narzut (`include` + JIT), więc koszt krańcowy to ~2.6 s na pulsara.
**Pełny przebieg 380 pulsarów: ~17 minut**, nie 3 godziny.

## 2026-09-17 — pełny przebieg 380 pulsarów: NIEZALEŻNE POTWIERDZENIE TRENDU

`batch_run.jl 1 380` — 380 pulsarów w ~20 min. Wynik: `~/claude/work/separations_batch_full.csv`,
log `logs/batch_full.log`, obrazki `~/claude/work/review/` (123 katalogi, 138 MB).

| status | liczba |
|---|---|
| rozjazd systematyczny Low↔High | 184 |
| ok, do przeglądu | 102 |
| ok, czyste | 47 |
| za mało pulsów | 26 |
| składowe zlane (n_comp=2 na jednoskładnikowym) | 21 |

**Plon 149 z 380 = 39%** — poniżej prognozy 55%, bo dominującym trybem awarii okazał się
rozjazd systematyczny (48% przypadków), rosnący wraz ze spadkiem S/N wzdłuż posortowanej listy.

### Wynik merytoryczny

Nowe 149 pulsarów jest **rozłączne** ze starymi 91, więc po raz pierwszy da się sprawdzić
trend z Ė na niezależnej próbce. Ė było wskazane przez poprzedni raport, więc jako hipoteza
postawiona z góry test nie wymaga korekty na liczbę porównań.

| próbka | N | ρ_S | p | przy ustalonym log P |
|---|---|---|---|---|
| stare (zweryfikowane) | 91 | +0.257 | 0.014 | +0.259, p=0.013 |
| **nowe czyste, rozłączne** | **46** | **+0.350** | **0.017** | +0.273, p=0.066 |
| nowe czyste bez ekstremum | 45 | +0.418 | 0.0043 | +0.341, p=0.022 |
| **stare + nowe czyste** | **137** | **+0.286** | **0.0007** | +0.259, p=0.0023 |
| stare + nowe czyste, bez ekstremum | 136 | +0.305 | 0.0003 | +0.282, p=0.0009 |

Z 0.018 na jednej próbce zrobiło się **0.0007 przy N = 137**, z niezależnym potwierdzeniem
po drodze. To już nie jest poszlaka.

### Zastrzeżenia, bez których ta liczba jest myląca

1. **Żaden z 149 nowych pomiarów nie był oglądany przez człowieka.** 46 „czystych" przeszło
   wyłącznie bramki automatyczne.
2. **Wynik zależy od odrzucenia 102 ze 149.** Z wszystkimi: ρ = −0.025, p = 0.76 — trend
   znika. To zachowanie oczekiwane (mediana σ(ΔW/W) w kolejce do przeglądu to 0.046 wobec
   0.029 dla czystych, i są tam wartości rzędu ΔW/W = −1.45), ale **to jest największy wybór
   analityczny w całej tej analizie i musi zostać zweryfikowany ręcznie.**
3. Brak dowodu, że cięcie jakościowe jest obciążone względem Ė: mediana log Ė 31.93 (czyste)
   vs 32.21 (do przeglądu), KS p = 0.41. Test ma jednak umiarkowaną moc.
4. Błędy nowych są 2× większe od starych (0.029 vs 0.013), ale wciąż 3× mniejsze od rozrzutu
   populacyjnego (0.096), więc niosą ~92% wagi — zgodnie z wcześniejszym rachunkiem.
5. **Bramki są niekompletne.** J1717-3425 dostał status „czysty" z ΔW/W = −0.74 ± 0.11, przy
   maksimum 0.32 w 91 zweryfikowanych. Dodana bramka `MAX_FRAC = 0.4` kieruje takie przypadki
   do przeglądu. Ekstremum **osłabia** trend, nie napędza go (bez niego p = 0.0043 zamiast
   0.017 na próbce niezależnej).

### Następny krok
Przejrzeć 123 katalogi w `~/claude/work/review/`. Dopóki to nie jest zrobione,
`separations_batch_full.csv` NIE jest scalany do `input/separations.csv`.

### Przegląd 123 katalogów kolejki — i KOREKTA wcześniejszego wniosku (2026-09-17)

Metoda: zamiast oglądać ~1500 obrazków per-puls, `scripts/review_sheets.py` parsuje
`print_fit_summary` z `logs/batch_full.log` i rysuje jeden panel na pulsara — ślad μ każdej
składowej pulse po pulsie, osobno 1023 i 1523 MHz, składowe posortowane po μ (tak jak robi
`component_offsets`). Arkusze po 12: `~/output/claude/review_sheets/sheet_01..11.png`.
Werdykty: `~/claude/work/review_verdicts.csv` (123 wiersze z uzasadnieniem).

Uwaga metodyczna: pierwsza wersja metryki liczyła „skrzyżowania" jako niezgodność kolejności
surowego wydruku fitu z posortowaną — to artefakt, bo `component_offsets` i tak sortuje.
Zastąpione liczbą skoków śladu po posortowaniu.

**Werdykt: 57 zachowanych, 66 odrzuconych** ze 123. Razem z 47 auto-czystymi daje to
**104 nowe pomiary** (z 149 przed przeglądem).

Dominujące powody odrzucenia: bistabilne etykietowanie (ślad przeskakuje między dwiema
pozycjami), pasmo high rozrzucone przy stabilnym low, ślady zbiegające się lub pokryte
(n_comp=2 na jednej składowej), trwały rozjazd low↔high w jednej składowej.

### KOREKTA: niezależne potwierdzenie nie przeżyło przeglądu

Wcześniejszy wpis raportował na 46 auto-czystych ρ = +0.350, p = 0.017 jako niezależne
potwierdzenie trendu z Ė. **Po dołożeniu 57 pulsarów zachowanych w przeglądzie sygnał
w próbce niezależnej słabnie poniżej istotności:**

| próbka | N | ρ_S | p | przy ustalonym log P |
|---|---|---|---|---|
| stare (zweryfikowane) | 91 | +0.257 | 0.014 | +0.259, p=0.013 |
| nowe, tylko auto-czyste (poprzedni wpis) | 46 | +0.350 | 0.017 | +0.273, p=0.066 |
| **nowe po pełnym przeglądzie** | **103** | **+0.158** | **0.110** | +0.089, p=0.37 |
| wszystko po przeglądzie | 194 | +0.207 | 0.0037 | +0.175, p=0.015 |

Czyli p = 0.0007, które raportowałem na próbce „stare + auto-czyste", schodzi do **0.0037**
przy N = 194, a **niezależne potwierdzenie przestaje być istotne (p = 0.11)**. Te 46
auto-czystych było podpróbką wybraną kryteriami (zero odrzuconych pulsów, niskie χ²), która
dała mocniejszą korelację niż pełna zweryfikowana próbka. Nie umiem rozstrzygnąć, czy to
przypadek, czy cięcie auto-czyste selekcjonuje coś, co korelację zawyża.

### Czego to nie podważa

Populacyjny charakter zwężenia jest odporny: **125 z 195 pulsarów zwęża profil (64%),
p = 1.0 × 10⁻⁴**, mediana ΔW/W = −0.024.

### Różnica rozkładów stare vs nowe

Nowe: mediana −0.009, IQR [−0.042, +0.027], 57% zwężeń. Stare: −0.045, IQR [−0.088, +0.003],
73%. KS p < 0.001. Nowe pulsary zwężają się słabiej. Prawdopodobna przyczyna: stare 91
pochodzą z `separations_todo.csv`, listy wyselekcjonowanej z przypadków, w których offset
dawał się zmierzyć i był istotny — czyli z obciążeniem w stronę dużych |ΔW/W|.

### Luka w bramkach
`MAX_FRAC` dodany PO przebiegu, więc nie zadziałał: **J1717-3425 (ΔW/W = −0.742)** przeszedł
jako auto-czysty i nie trafił do kolejki, mimo że maksimum w 91 zweryfikowanych to 0.32.
Nie był oglądany. Bez niego: nowe p = 0.064, wszystko p = 0.0020. Do przejrzenia ręcznie.

Scalona tabela: `~/claude/work/separations_final_merged.csv` (195 wierszy, kolumna `zrodlo`:
stare / auto / przeglad). Nadal NIE scalone do `input/separations.csv`.

## 2026-09-18 — Porównanie `separations_maciej.csv` z istniejącymi pomiarami

Cel: wczytać plik Macieja (50 pulsarów) do `input/` i skonfrontować z `separations.csv`,
`separations_merged.csv`, `separations_todo.csv`.

Skrypty: `work/scripts/cmp_maciej.py` (tabela porównań), `work/scripts/profile_check.py`
(profile średnie). Log: `work/logs/cmp_maciej.log`. Wykres: `~/output/claude/cmp_maciej_profile.png`.

### Pokrycie
Wszystkie 50 pulsarów Macieja są w `separations_todo.csv` (i wszystkie mają puste `zrobione`).
47 jest w `separations.csv`, 48 w `merged`. Nowych pulsarów brak.
Brak w `separations.csv`: J1627-5936, J1819+1305, J1932+1059.

### Zgodność
41/47 zgodnych z `separations.csv` w granicach 3σ, w tym 32 identyczne co do 1e-4 stopnia.
`ncomp` zgodne wszędzie.

Rozbieżne (>3σ), zawsze przez JEDEN skrajny komponent, nie przez przesunięcie całego profilu:
J1733-3716 (12.1σ), J1901+0716 (10.1σ), J1757-2421 (7.8σ), J1714-1054 (5.9σ),
J1803-3329 (4.2σ), J1808-3249 (3.9σ). Dodatkowo vs `merged`: J1819+1305 (11.9σ).

### Rozstrzygnięcie z zapisanych fitów per-pulse (`<PSR>_16/component_offsets.txt`)
Odtworzenie procedury `_offset_summary` z surowych mu wszystkich pulsów:

| pulsar | rekonstrukcja | separations.csv | Maciej | merged |
|---|---|---|---|---|
| J1757-2421 | 18.652 ± 0.132 | 18.7415 | 20.6813 | 18.7415 |
| J1803-3329 | 4.075 ± 0.082 | 4.1425 | 4.7241 | 4.1425 |
| J1808-3249 | 9.591 ± 0.112 | 9.6328 | 11.1026 | 9.6328 |
| J1819+1305 | 17.403 ± 0.134 | — | 18.5010 | 13.3531 |

Wniosek: dla trzech pierwszych `separations.csv` odtwarza się z surowych fitów (drobne różnice
= maska `keep`), wartości Macieja nie. **J1819+1305 to osobny problem: wartość w `merged`
(13.3531) nie wychodzi z zapisanych fitów dla ŻADNEJ pary komponentów** — lony to
169.74 / 175.30 / 187.14. 13.35 leży na skraju zakresu osiągalnego tylko przy bardzo agresywnym
odrzuceniu pulsów. Do sprawdzenia, skąd pochodzi.

### Profile średnie (pdv -FTt na pulsar.low/high)
- J1733-3716: drugi komponent ma szczyt 201.5 (low) / 200.4 (high) → `separations.csv` (201.87)
  trafia w szczyt, Maciej (196.69) siedzi na zboczu narastającym.
- J1714-1054: drugi komponent słaby (0.31 low, 0.17 high), szczyty 185.6 / 184.6 → nierozstrzygalne,
  obie wartości (184.86 / 185.47) w obrębie komponentu.
- J1901+0716: sporny komponent to szerokie ramię bez własnego maksimum → nierozstrzygalne z profilu
  średniego, brak `component_offsets.txt`.

Otwarte: (1) skąd 13.3531 dla J1819+1305 w merged; (2) czy wciągać dane Macieja jako osobne
`zrodlo` — na razie NIE scalone, plik leży jako `input/separations_maciej.csv`.

---

## 2026-09-22 — Test „travel": dryf vs modulacja amplitudowa bez użycia P3

**Cel.** Nowa metoda rozróżniania dryfu od P3-only, niezależna od kryterium Song et al. 2023
(offset centroidy mocy w 2DFS względem osi 1/P2 = 0). Motywacja: centroida mocy jest obciążona
przez stochastyczną zmienność kształtu pulsu, która piętrzy moc wzdłuż 1/P2 = 0.

**Metoda.** Dwuwymiarowa autokorelacja fluktuacji, K(Δ,τ) = Σ δI(n,φ)·δI(n+τ,φ+Δ), i jej część
antysymetryczna A(Δ,τ) = K(Δ,τ) − K(−Δ,τ) — „czy długość φ wyprzedza φ+Δ, czy odwrotnie".
Modulacja amplitudowa jest separowalna, δI = a(φ)w(n), więc K = C_a(Δ)·C_w(τ), a C_w jest
**dokładnie** parzysta w τ (przeindeksowanie sumy, nie stacjonarność). Stąd A ≡ 0 tożsamościowo:
dla dowolnego w(n) — wędrujące P3, nulling, brak okresowości — i dla a(φ) zmieniającego znak,
czyli dla antyfazowej AM, która w `phase_modulation3` czyta 25σ fałszywego dryfu.

Formalnie ten sam kanał informacji co asymetria 2DFS, ale estymator rzutuje część symetryczną
do zera algebraicznie zamiast ją uśredniać w centroidzie, i całkuje po całej płaszczyźnie (Δ,τ)
zamiast po jednym binie f3.

**Implementacja.** `modules/travel.jl` (moduł `Travel`), `Plot.travel`, wrapper
`SpaTs.travel_test(outdir; ...)` z parametrem `datafile` (obsługuje też `pulsar_high_debase.txt`
z katalogów `_16`). Null: surogat separowalny rank-1 + szum bootstrapowany z ciągłych pasków
off-pulse z losowym przesunięciem cyklicznym w czasie. Kontrola off-pulse wbudowana.

**Walidacja** (`work/scripts/travel_check.jl`, `travel_stress.jl`, `Travel.selftest()`):

| test | wynik |
|---|---|
| FFT vs suma wprost | zgodność 4e-16 |
| pole separowalne (wędrujące P3 + nulling + skok znaku) | max\|A\| = 1.1e-16 × K(0,0) |
| wzór bieżący syntetyczny | P2 = 18.0 (prawda 18), P3 = 12.0 (prawda 12), rank1 = 1.00 |
| J0820-1350 (B0818−13, znany dryfer) | >999σ, rank1 = 0.999, P3 = 4.8 vs katalog 4.78 |
| J1907+0731 (P3-only) | 1.7σ (p = 0.07), rank1 = 0.09, bloki ≈ 0.00 |
| J2053-7200 (wobble P3 opróżnia bin LRFS) | **110σ** — `phase_modulation` daje tam 0.5σ |
| J1750-3503 (dryf odwracający kierunek) | 163σ, bloki pokazują oba znaki wprost |
| J1110-5637 (silny koherentny dryf) | 81σ |

Kontrola off-pulse przeszła wszędzie (|σ| ≤ 2.1).

**Uwagi.** (1) Przy bardzo silnym dryfie (T − ⟨T⟩)/σ eksploduje do 10⁵ — sensowną liczbą jest
p-value ograniczone przez `nreal`; printout i rysunek zakrywają to jako „>999σ". (2) P2 i P3
z przejść przez zero modów SVD są **poglądowe** (relacja P/2 jest ścisła tylko dla czystej
sinusoidy; harmoniczne i zanik koherencji przesuwają przejście o kilka–20%). (3) Projekcja
blokowa liczona leave-one-out — projekcja na mapę globalną zawierałaby człon własny i dawała
~1/√n_bloków dla samego szumu (początkowo dawało to mylące 0.48 dla J1907+0731).

**Otwarte:** batch po 114/115 pulsarach P3-only (katalogi `<PSR>_16/`, plik
`pulsar_high_debase.txt`); kontrola pozytywna na 418 dryferach z `input/drift_pulsars_P3.txt`;
analiza populacyjna (π₀ z rozkładu p-value) zamiast 115 niezależnych progów 3σ.

### 2026-09-22 (cd.) — kierunek degradacji: `T_inc` i projekcja matched

**Problem.** Test travel w wersji podstawowej nie może zdegradować etykiety `drift` do P3-only:
niska `T` to brak dowodu, nie dowód braku. Dodatkowo globalna mapa A ma ślepy punkt — dryfer
odwracający kierunek **symetrycznie** kasuje się w sumie i czyta jako brak ruchu.

**Dodane w `Travel.travel_test`:**

1. `T_inc = Σ_b Σ A_b²` — suma niekoherentna po blokach impulsów. Zamyka ślepy punkt reversera.
   Zweryfikowane w `selftest`: dla zbalansowanego reversera **T/T_inc = 3.4e-6**. Płaci za to
   wyższym progiem szumu (nb bloków szumu zamiast jednego uśrednienia), więc dla stałego dryfu
   jest mniej czuła — raportować obie, przy braniu lepszej doliczyć karę za trials.
2. `drift_template(max_lag, max_dphi, p2, p3; npulses, non)` + pola `frac`, `frac_err`,
   `frac_sig`, `frac_limit` — projekcja matched na mapę, jaką dałby dryf o **deklarowanym** P2:
   `frac = ⟨A, szablon⟩ / (K(0,0)·‖szablon‖²)` = ułamek mocy modulacji siedzący w dryfie o tej
   geometrii. To pozwala hipotezę dryfu **odrzucić**, nie tylko nie potwierdzić.

**Błąd wyłapany przez selftest:** pierwsza wersja szablonu dawała `frac` = 0.767 zamiast 1 dla
syntetycznego dryfu o znanym P2/P3. Przyczyna: korelacja liniowa sumuje po (N−τ)(M−|Δ|) parach
przy opóźnieniu (Δ,τ) wobec NM w zerze, więc mapa jest przycięta trójkątnym taperem nawet dla
idealnie koherentnego wzoru. Taper jest czystą geometrią, znaną dokładnie — po wstawieniu go do
szablonu `frac` = 1.000. Bez tego projekcja zaniżała o 23% przy N=600, M=40.

**Wyniki projekcji matched:**

| przypadek | frac | limit 3σ |
|---|---|---|
| J0820-1350, szablon z właściwym P2 = −13.7 | +0.118 | — |
| J0820-1350, odwrócony znak P2 = +13.7 | −0.118 | — |
| J0820-1350, błędny P2 = +45 | −0.042 | — |
| J1907+0731 (P3-only), wmówiony dryf P2 = +20 | 0.0004 ± 0.0001 | **0.0007** |

Czyli zakres dynamiczny prawdziwy dryfer / P3-only ≈ **300×**. `frac` prawdziwego dryfera to 0.118,
nie 1 — szablon jest pojedynczą sinusoidą, a realna modulacja ma harmoniczne i traci koherencję,
więc `frac` jest dolnym ograniczeniem. **Próg degradacji (`demote_frac`, domyślnie 0.05) jest
placeholderem** — trzeba go skalibrować przebiegiem po `input/drift_pulsars_P3.txt`.

**Otwarty problem do rozstrzygnięcia:** J1907+0731 (kontrola negatywna) daje `T_inc` = 4.8–5.0σ przy
`T` = 1.7σ. Kontrola off-pulse dla obu statystyk czysta (T_off = −0.3σ, **T_inc_off = 0.2σ**), więc
to nie jest asymetria szumu ani zły null. Najbardziej prawdopodobna diagnoza: null używa surogatu
**rank-1** (separowalnego), a jeśli pole jest rank ≥ 2 z niezależnymi modami czasowymi, surogat
zaniża wariancję i zawyża istotność. Wskazówka zgodna: projekcje blokowe ≈ 0.00, czyli mapy bloków
**nie zgadzają się ze sobą** — to sygnatura losowego uporządkowania między niezależnymi modami,
a nie dryfu. Możliwa poprawka: surogat rank-r z niezależnie przesuwanymi w czasie modami.
Do czasu rozstrzygnięcia: `T_inc` bez zgodności blokowej nie jest kandydatem na dryf.

Werdykt degradacji jest zablokowany, gdy `T_inc` > 3σ — sprawdzone, dla J1907+0731 nie drukuje się.

### 2026-09-22 (cd. 2) — iloraz R: próg bez kalibracji na skażonej próbce

**Problem podniesiony przez AS.** Kalibracja progu degradacji na 418 pulsarach `drift` jest
cyrkularna — to jest właśnie próbka podejrzana o skażenie. Kierunek biasu: ucząc się na mieszance
dowiesz się, że „dryfer może mieć frac ≈ 0", próg spadnie i nic nie zdegradujesz. Nie psuje to
kontroli fałszywych pozytywów, tylko kasuje moc testu.

**Rozwiązanie: normalizacja z tego samego pulsara.** Model dryfu rozkłada się na dwie połowy
o **równych współczynnikach**:

```
cos(2π(Δ/P₂ − τ/P₃)) = cos(2πΔ/P₂)cos(2πτ/P₃) + sin(2πΔ/P₂)sin(2πτ/P₃)
                       └─ parzysta w Δ ─┘        └─ nieparzysta w Δ ─┘
```

Nieparzystą mierzy `A` (antisym_map), parzystą nowa `E` (`sym_map`). Modulacja amplitudowa, będąc
separowalną, wkłada wszystko w parzystą i nic w nieparzystą. Stąd `R = frac_odd / frac_even` = 1
dla dryfu, 0 dla AM — a siła modulacji, harmoniczne, zanik koherencji i taper mnożą obie połowy
tak samo i **kasują się w ilorazie**. Zero liczb z zewnątrz.

**Weryfikacja syntetyczna** (`selftest`): czysty dryf R = 1.000. Fala bieżąca + stojąca o tych
samych okresach: R = 0.607 przy przewidzianym analitycznie 0.600 (fala stojąca to pół bieżącej
w przód + pół w tył, więc R = (a_f²−a_b²)/(a_f²+a_b²)).

**Weryfikacja na danych — to jest właściwy wynik:**

| pulsar | frac_odd | frac_even | R | klasa |
|---|---|---|---|---|
| J0820-1350 | 0.1180 | 0.1067 | **1.106 ± 0.002** | dryfer |
| J1110-5637 | 0.0062 | 0.0075 | **0.827 ± 0.034** | dryfer |
| J2053-7200 | 0.0047 | 0.0036 | **1.305 ± 0.025** | dryfer |
| J1750-3503 | 0.0114 | 0.0179 | **0.640 ± 0.003** | dryfer (reverser) |
| J1907+0731 | 0.0004 | 0.0002 | NaN | P3-only |

Sedno: **`frac_odd` rozciąga się na 25× między dryferami (0.0047–0.118), a R tylko na 2×
(0.64–1.31)**. Na progu bezwzględnym J1110-5637 i J2053-7200 wyglądałyby jak „prawie zero" obok
J0820-1350 i zostałyby błędnie zdegradowane. R je poprawnie trzyma przy 1. To jest dokładnie ta
zmienność, której nie da się wyuczyć ze skażonej próbki — i którą R usuwa bez próbki.

Próg `demote_R` = 0.3 leży w luce między obserwowanym pasmem dryferów a zerem z konstrukcji dla AM.

**Zastrzeżenia (udokumentowane w kodzie):**
- R **nie jest ścisłym ułamkiem i może przekroczyć 1** (tu 1.11 i 1.31). Mianownik zbiera też
  parzystą projekcję modulacji nie-wędrującej, a ta nie ma ustalonego znaku. Pierwsza wersja testu
  syntetycznego dała R = 1.32 dla dryfu + rampy monotonicznej. Odczyt jest porządkowy, nie
  dosłowny. **Patologia siedzi przy górnym końcu, a degradacja rozgrywa się przy dolnym**, gdzie
  do zera dąży licznik i iloraz zachowuje się dobrze.
- **Reverser zaniża R** (J1750-3503: 0.64, najniżej z czwórki) — globalna mapa częściowo się kasuje.
  Przed odczytem R sprawdzić `T_inc` i `block_proj`, a przy obu znakach ograniczyć się do jednego
  epizodu przez `pulse_st`/`pulse_end`.
- J1907+0731 daje NaN, bo `frac_even` nie jest zmierzone na 3σ — przy wmówionym P2 = 20 nie ma
  koherentnej modulacji, więc nie ma mianownika ani hipotezy dryfu do odrzucenia. Guard zadziałał.

**Błąd wyłapany przy okazji:** zmiana sygnatury `_travel_maps` na 4-krotkę zepsuła rozpakowanie
w selfteście (test reversera czytał 2000 zamiast 3.4e-6). Selftest to złapał.

**Następny krok:** batch po 533 (418 + 115) z tabelą, rozkład R, sprawdzenie bimodalności.
Etykiety Song et al. wchodzą wtedy jako **zbiór testowy**, nie treningowy.

### 2026-09-23 — pełny przebieg 533, przejście na ρ, kontrola okien

Opis metody: `docs/travel_test_method.md` (zaktualizowany do stanu na dziś).

**Zmiany metodologiczne.** `R = odd/even` zastąpione przez `ρ = √2·odd/hypot(odd,even)` — ten sam
kąt, ale sinus zamiast tangensa, więc bez bieguna, bez NaN-ów i **bez dwóch progów jakości**, które
byłyby potrzebne tylko do tłumienia blow-upu. Geometria: |P₂| skanowane po projekcji parzystej,
P₃ brane z LRFS (swobodny fit P₃ ucieka w róg dużych P₂/P₃ — J2053-7200 dopasowało 63 zamiast 3.06).

**Dane.** Pierwszy przebieg mieszał pełne pasmo (85 pulsarów) z podpasmem `low` = 3/16 (430) —
różnica czułości 2.25×, niedopuszczalna. Pełne pasmo odtworzone z `pulsar.spCf16` i **zostaje w
katalogach `_16`** jako `pulsar.full` / `pulsar_full.debase.gg` / `pulsar_full_debase.txt`
(~72 MB/pulsar, ~38 GB, QNAP ma 1.8 TB wolnego).

**Wynik.** 533 policzone, 515 bez błędu, 3 odrzucone przez kontrolę off-pulse, 387 z detekcją.
W reżimie rozdzielczym (P₂fit ≤ M/2): **mediana ρ dryferów = 1.011 przy n = 172**, wobec
przewidywania teorii dokładnie 1, bez ani jednego parametru swobodnego. p3only: 0.324 (n=7).
Sama geometria też rozróżnia: mediana P₂fit/M = 0.39 (drift) wobec 1.39 (p3only).

Pozorna zależność ρ od S/N okazała się **paradoksem Simpsona** — w obu grupach geometrii ρ jest
płaskie, S/N steruje tylko przynależnością do grupy.

**Kontrola okien on-pulse** (`check_onpulse.jl`, profile średnie z `pdv -t -F -T`): M medianowo
116 binów wobec zasięgu emisji 3σ = 56, szczyt wewnątrz okna u 512/515. Okna są hojne, ale nie
błędne. Syntetyk potwierdza, że **poszerzanie okna szkodzi** (rozcieńcza sygnał, zawyża P₂ — przy
M = 250 zamiast 100 P₂fit wychodzi 309 zamiast 150) i że większy zakres Δ też nie pomaga.
Informacja o okresie w długości pochodzi wyłącznie z obszaru, który świeci.

**Kluczowe zastrzeżenie — test stabilności okna.** ρ liczone przy 1.0/1.5/2.0 × W₃σ: prawdziwe
dryfery są odporne (J0034-0721: 1.042/1.067/1.048), ale obiekty o niskim ρ **chwieją się nawet
40-krotnie** (J1239+2453: 0.011/0.374/0.461; J2048-1616: 0.021/0.003/0.481). Więc walidacja
populacyjna stoi, ale **żadna lista kandydatów nie jest wiarygodna bez testu stabilności** — w tym
J1810-5338, najczystszy kandydat do promocji, który go nie przeszedł (0.566/0.438/0.770).

**Otwarte:** (1) test stabilności w pipeline i przeliczenie kandydatów; (2) niewyjaśnione ρ > 1 u
25 dobrych dryferów (J1519-6106: 1.315 ± 0.008 przy rank1 = 0.97) — model przewiduje ρ ≤ 1;
(3) null rank-1 wobec pola rank ≥ 2 (J1907+0731); (4) P₂ do zamiany na stopnie.

### 2026-09-23 (cd.) — przestawienie hierarchii: T/T_inc jako produkt główny, ρ jako dodatek

Po serii testów wyszło, że myliłem dwa różne pytania:

- **A: czy jest jakikolwiek ruch?** — odpowiada `T`/`T_inc`, bez żadnych założeń (tożsamość A ≡ 0).
  To jest binarne pytanie, które stawia klasyfikacja Song et al.
- **B: czy modulacja to zasadniczo sztywna translacja?** — odpowiada ρ, ale wymaga mierzalnej
  geometrii.

ρ postawiłem jako produkt główny i to był błąd hierarchii. Wynik na pełnej próbce dla pytania A:
detekcja ≥5σ u **330/405 dryferów (82%)** i **57/107 P3-only (53%)**. Tych 53% nie wolno jednak
ogłaszać przed naprawą nullu rank-1, bo część może być artefaktem modelu zerowego.

**Priorytet numer jeden to teraz surogat rank-1.** Dotąd był punktem na liście otwartych; skoro
produktem głównym jest detekcja, to od niego zależy wszystko. Objaw wzorcowy: J1907+0731, T_inc =
4.8σ przy czystych kontrolach off-pulse i projekcjach blokowych ≈ 0.

**Wariant bez szablonu — sprawdzony i odrzucony.** Stosunek norm √(‖A‖²/‖E‖²) miał usunąć fit P₂.
Na syntetyku działa i jest całkowicie odporny na okno (0.819/0.829/0.829/0.829 tam, gdzie wersja
szablonowa spada 0.872 → 0.653), ma też analityczną korektę na próbkowanie τ (przy P₃ = 2.05 surowe
0.407 → po korekcie 1.017). Ale na danych realnych rozdzielczość spada z 2.9× na **1.25×**, a
J0034-0721 (podręcznikowy dryfer) czyta **0.203**. Powód: ‖E‖² zbiera całe tło mapy, a ‖A‖² nie, bo
A znika przy Δ = 0. **To jest odpowiedź na pytanie, po co w ogóle fit P₂: szablon jest jedyną
rzeczą, która wycina część związaną z dryfem i odrzuca resztę.** Syntetyk tego nie pokazał, bo
pojedyncza sinusoida daje A i E identyczną strukturę.

**Mechanizm zależności od okna — rozpracowany.** Okno nie psuje ρ bezpośrednio (przy P₂ ustalonym
na sztywno ρ jest płaskie). Okno przeciąga **fit P₂** (150 → 309 przy M = 100 → 250), a zawyżone P₂
kaleczy asymetrycznie: przy Δ → 0 cos → 1, a sin → 0, więc kanał nieparzysty traci pokrycie z
sygnałem. Naprawa zweryfikowana: `max_dphi` z W₃σ zamiast z M usztywnia P₂fit na 150.2 we wszystkich
szerokościach okna. **Niewdrożone w pełnym przebiegu.**

**Fizyka: P₂ > W jest dopuszczalne** (mało iskier albo linia widzenia prostopadła do pierścienia).
Wtedy widać najwyżej jeden podpuls naraz i **P₂ nie jest mierzalne, tylko ograniczone od dołu** —
ucieczka fitu na kraniec siatki jest poprawną odpowiedzią, a nie usterką. Wartości P₂fit > M to
limity. W tym reżimie T/T_inc działają, ale z ρ nie wolno robić degradacji.

**Wędrujące P₃ podnosi ρ powyżej 1** (syntetyk: 1.387 przy P₃ wahającym się 3–7, T = 8530σ). To
pierwsza hipoteza tłumacząca ogon ρ > 1.1 u 25 realnych dryferów.

Dokumentacja `docs/travel_test_method.md` przepisana wokół nowej hierarchii; §11 to lista ograniczeń
uporządkowana wg ważności, §12 zawiera sześć wycofanych wniosków.

### 2026-09-23 (cd. 2) — dudnienie udaje dryf; to problem klasyfikacji, nie kalibracji

**Demonstracja.** Dwie nakładające się składowe, każda z własnym P₃, **nic się nie przemieszcza**:

| P₃ składowych | dudnienie | T | T_inc | spójność blok. |
|---|---|---|---|---|
| 7.0 / 7.4 | 129 P | **37.9σ** | **73.3σ** | **0.95** |
| 7.0 / 9.0 | 32 P | −0.1σ | 9.7σ | 0.76 |
| 7.0 / 13.0 | 15 P | 0.9σ | 1.0σ | 0.20 |
| rank-1 (ścisłe H₀) | — | −0.3σ | −1.4σ | 0.19 |

Groźne są bliskie okresy — długie dudnienie nie zdąży się uśrednić. **Reguła „T_inc bez zgodności
blokowej nie jest kandydatem" jest niewystarczająca**: najgorszy przypadek ma zgodność 0.95, bo
dudnienie dwóch ściśle okresowych sygnałów jest deterministyczne i powtarza się z bloku na blok.

**Próba naprawy nullu (droga 1) — odrzucona.** Surogat rank-r z randomizacją faz Fouriera
(`rank_r_modes`, `phase_randomize!`, `surrogate_rank`). Fałszywki znikły (74σ → −1.2σ), ale
**prawdziwy dryf spadł z 19112σ na 1.9σ**. Powód jest nieusuwalny: dryf *jest* relacją fazową
między dwoma modami w kwadraturze, więc randomizacja faz losuje dokładnie to, co stanowi sygnał.
Poszerzając null tak, by objął dudnienie, obejmuje się nim również dryf. Domyślny `surrogate_rank`
przywrócony na 1; opcja zostaje w kodzie jako zapis sprawdzonego wariantu.

**Przeformułowanie.** To nie był błąd kalibracji. Null rank-1 jest poprawny (dla pola rank-1 daje
−0.3σ), liczba 74σ jest prawdziwa — w danych naprawdę jest uporządkowanie czasowe. Błędny był mój
krok **„jest uporządkowanie ⇒ jest dryf"**. Rozdzielenie musi nastąpić **po** detekcji.

**Właściwy dyskryminator (droga 2) — działa.** Dryf ma wyprzedzenie trwałe, dudnienie odwraca znak
co pół okresu dudnienia:

| | 4 bloki (250 P) | 10 bloków (100 P) |
|---|---|---|
| dudnienie (129 P) | +0.98 +0.97 +0.97 +0.95 | −0.12 −0.68 −0.99 −0.97 −0.05 +0.93 −0.99 −0.99 +0.90 +0.61 |
| prawdziwy dryf | +1.00 ×4 | **+1.00 ×10** |

**Blokada:** `nblocks` jest przycinane do `N ÷ (4·max_lag)`, co przy realnych danych daje maksymalnie
4 bloki — dokładnie reżim, w którym dudnienie udaje spójność. **W pełnym przebiegu ten test nie miał
szans zadziałać.** Wdrożenie skanu wymaga poluzowania przycięcia, a to krótszego `max_lag` dla
krótkich bloków. Do rozwiązania.

Dokumentacja: §7 przepisany (7.1 konstrukcja, 7.2 dudnienie, 7.3 odrzucona naprawa, 7.4 skan po
długości bloku), §11 pkt 1 przeformułowany, §12 dopisany siódmy wycofany wniosek.

### 2026-09-24 — pełny przebieg v2: `max_dphi` z W₃σ + skan spójności

**Cel.** Wdrożyć w pełnym przebiegu obie zaległe zmiany (§8.3 i §7.4 opisu metody) i dopiero
wtedy ocenić frakcję detekcji wśród P3-only.

**Polecenia.** `~/claude/work/scripts/travel_batch_full.jl` z `max_dphi = clamp(W₃σ ÷ 2, 1, M ÷ 2)`
(W₃σ z `onpulse_check.csv`, jest dla wszystkich 515) → `~/output/claude/travel_batch_v2.csv`
(52 min, pilot 6 pulsarów: `travel_v2_pilot.csv`). Podsumowanie:
`python3 ~/claude/work/scripts/travel_v2_summary.py` → `~/claude/work/logs/travel_v2_summary.log`.
Stary `travel_batch_full.csv` zostaje do porównania.

**Kontrole.** 18 błędów (te same co wcześniej), 2 odrzucone off-pulse. Model zerowy nadal
skalibrowany: σ off-pulse T 1.03, T_inc 0.91 (było 0.91/0.85), mediany ≈ 0. `max_dphi` medianowo
0.49 × M/2.

**Pytanie A — detekcja ≥5σ:**

| | stary (M/2) | nowy (W₃σ/2) | zmiana |
|---|---|---|---|
| drift | 330/405 (81.5%) | **368/406 (90.6%)** | +37, −0 |
| p3only | 57/107 (53.3%) | **79/107 (73.8%)** | +22, −0 |

Istotność rośnie medianowo 1.8×. Nowe detekcje pochodzą z pulsarów o oknie najmocniej
przewymiarowanym (`max_dphi` 0.35 × M/2): stary zakres Δ dokładał same biny szumowe i rozcieńczał T.

**Skan spójności — rozstrzyga, czym są te detekcje.** `min(cons)` wśród detekcji, wg siły detekcji:

| siła | drift: n, mediana | p3only: n, mediana |
|---|---|---|
| 5–20σ | 54, 0.086 | 31, 0.018 |
| 20–100σ | 101, 0.231 | 23, 0.035 |
| >100σ | 213, **0.419** | 25, **0.088** |

Przy porównywalnej sile detekcji P3-only mają spójność kilkakrotnie niższą, a ~⅓ z nich ma
`min(cons) < 0`. Poniżej 0.2 jest 71/79 detekcji p3only wobec 129/368 dryferów. Czyli: **74% P3-only
nie jest czystą modulacją amplitudową, ale ich uporządkowanie czasowe w zdecydowanej większości
się nie odtwarza** — nie ma trwałego wyprzedzania jak u dryferów. Przejście przez zero (sygnatura
dudnienia) przy detekcji ≥20σ: 0/314 drift, 1/48 p3only (J1001-5939: 0.34 → −0.05).

**ρ (P₂fit ≤ M/2, z detekcją):** drift n = 200, mediana **0.990** (kw. 0.733–1.082), było 1.011
przy n = 172; p3only n = 9, mediana 0.202.

**Otwarte.** (1) Ciąg „detekcja → spójność” potrzebuje progu albo modelu odniesienia dla `min(cons)`
zależnego od S/N — dotąd tylko rozkład; kandydat: ten sam skan na surogatach rank-1 i na
syntetycznym dryfie o danym S/N. (2) Frakcji 74% nie wolno podawać bez drugiego wiersza tabeli
(spójność). (3) Reżim P₂/M ≤ 0.5 wciąż liczony względem zadeklarowanego M, nie W₃σ.

### 2026-09-24 (cd.) — odniesienie dla min(cons) → null nie jest skalibrowany dla zmienności nieseparowalnej

**Cel.** Wyznaczyć odniesienie dla `min(cons)` zależne od siły detekcji (§11 pkt 1b).

**Polecenia** (`~/claude/work/scripts/`, wyniki w `~/output/claude/`):
- `travel_cons_ref.jl` → `travel_cons_ref.csv`: 30 losowych dryferów z sig ≥ 100 + biały szum
  (k·σ_off, k = 0…16) oraz syntetyk (sztywny dryf / kierunek losowy w epizodach 30 P / rank-1).
  Analiza: `travel_cons_ref_summary.py` → `travel_cons_ref.png`, log `travel_cons_ref_summary.log`.
- `travel_null_highsnr.jl` → `travel_null_highsnr.csv` (rank-1 przy amp 0.5–5, 10 ziaren) oraz
  `travel_modstrength.csv`: `mod = var(on)/var(off) − 1` po highpass dla wszystkich pulsarów.
- `travel_null_nonsep.jl` → `travel_null_nonsep.csv`: pola **bez uporządkowania czasowego** (impulsy
  niezależne), nieseparowalne: jitter, losowe podpulsy, losowe podpulsy + AM z P₃.

**1. `min(cons)` silnie zależy od S/N nawet dla idealnego dryfu.** Syntetyczny sztywny dryf:
0.07 przy 16σ, 0.20 przy 59σ, 0.57 przy 452σ, ~0.93 przy >1000σ. Zdegradowane dryfery: p10 < 0
aż do 100σ. Poniżej ~50σ `min(cons)` nie rozróżnia niczego.

**2. GŁÓWNY WYNIK: null rank-1 nie jest skalibrowany dla nieseparowalnej zmienności impuls-do-impulsu.**
Rank-1 jest w porządku do mod = 336 (0/40 realizacji ≥ 5σ). Ale pola bez żadnego ruchu:

| pole | mod ≈ 1 | mod ≈ 5 | mod ≈ 20 | mod ≈ 60–90 |
|---|---|---|---|---|
| jitter | 8–17σ | 32–82σ | 173–337σ | — |
| losowe podpulsy | 23–59σ | 146–254σ | 514–1071σ | 1900–4200σ |
| losowe podpulsy + AM | 39–73σ | 189–397σ | 601–1174σ | 2600–5000σ |

Kontrola off-pulse wszystko przepuszcza. `min(cons)` ≈ 0 lub ujemne w każdym przypadku.

Mechanizm: tożsamość §2 mówi, że **wartość oczekiwana** A dla takiego pola jest zero, ale T = ΣA²
zbiera też **wariancję** A. Surogat (wiodący mod SVD + szum z off-pulse) nie zawiera nieseparowalnej
zmienności on-pulse, więc ta wariancja nie wchodzi do rozkładu zerowego, a T rośnie liniowo z mod.
Obwiednia dla losowych podpulsów w tym syntetyku: sig ≈ 60·mod.

**Konsekwencje.**
- Detekcja T/T_inc ≥ 5σ **nie jest** dowodem uporządkowania czasowego przy silnej modulacji. Frakcje
  91% / 74% z przebiegu v2 **nie mają interpretacji** „nie jest czystą modulacją amplitudową”.
- Poniżej obwiedni (σ ≤ 60·mod): 59/78 detekcji p3only, 176/361 drift (orientacyjnie, bo obwiednia jest
  z jednej geometrii syntetyku).
- Skan spójności **jest** odporny na ten efekt (fluktuacje A z niezależnych impulsów nie odtwarzają się
  między blokami), więc to on, a nie T, mierzy trwałe uporządkowanie.

**Otwarte — wymaga decyzji.** Naprawa nullu, kandydaci:
(a) statystyka krzyżowa między blokami `T_cv = Σ_{b≠b'} ⟨A_b, A_b'⟩`: wartość oczekiwana 0 dla każdego
pola z fluktuacjami niezależnymi między blokami, więc nie wymaga modelu zmienności; wariancja z
jackknife po blokach. Traci reversera (jak T).
(b) surogat z permutacją kolejności impulsów: A części separowalnej zostaje ≡ 0 przy każdej kolejności,
a niezależna zmienność on-pulse jest zachowana. Ale niszczy też korelacje symetryczne w czasie.

### 2026-09-24 (cd. 2) — naprawa (a): statystyki krzyżowe między blokami, przebieg v3

**Konstrukcja** (`Travel.crossblock_test`, pola `cv_*` w wyniku `travel_test`).
`T_cv = Σ_{b≠b'} ⟨A_b, A_b'⟩ = ‖ΣA_b‖² − Σ‖A_b‖²` — diagonala, w której siedzi obciążenie od wariancji,
wycięta wprost. Null z **randomizacji znaków bloków**: brak uporządkowania = symetria względem odwrócenia
czasu, która daje A_b → −A_b. Wariancja analitycznie `2Σ_{b≠b'}G²`, p z 10⁵ losowań (dokładnie dla
2^(B−1) ≤ 10⁵). Bez surogatów i bez modelu zmienności. Cena: istotność ograniczona przez liczbę bloków,
z_max = √(B(B−1)/2) (22.3 przy B = 32). Dodatkowo `T_adj = Σ_b ⟨A_b, A_{b+1}⟩` (sąsiednie bloki), dodatni
dla ruchu dłuższego od bloku niezależnie od kierunku.
Błąd po drodze: `2^(B−1)` przepełnia Int64 przy B = 64 → nieskończona enumeracja; naprawione
(`B − 1 ≤ floor(log2(nflip))`).

**Walidacja** (`~/claude/work/scripts/travel_cv_validate.jl` → `travel_cv_validate{,_adj}.csv`,
9 typów × 9 amplitud × 4 ziarna):

| przypadek | T (null rank-1) | T_cv, nb = 32 | T_adj, nb = 64 |
|---|---|---|---|
| rank-1, jitter, losowe podpulsy, podpulsy+AM | do 1574σ | skalibrowane: p<0.05 w 2.2%, max z = 3.15 (n = 432 z dudnieniami) | max z 1.98 |
| dudnienie 7.0/7.4 | do 905σ | z ≈ −0.7 (przy nb = 8: z ≈ 5 — bloki = okres dudnienia) | **z ≈ 5.5** |
| dudnienie 7.0/9.0 | do 152σ | z ≈ 0 | z ≈ −6 (nb = 32: +4.4) |
| sztywny dryf | 8.9σ przy amp 0.2 | z = 4.7; od amp 0.3 z = 14–22 | 1.2–5.6 |
| dryf, P₃ ±20% | 13σ przy amp 0.3 | z = 12 | — |
| reverser (epizody 30 P / 100 P) | do 10⁵σ | **nie widzi** (z ≈ 0) | **z = 5–7** |

T_adj nie odróżnia reversera od dudnienia — lokalnie to ta sama rzecz (§7.2). Nie wnosi nic ponad T_cv na
danych (niżej), zostaje jako diagnostyka. Progi: z_cv(nb = 32) ≥ 5, z_adj(nb = 64) ≥ 5.

**Przebieg v3** (`travel_batch_full.jl` → `~/output/claude/travel_batch_v3.csv`, 53 min; podsumowanie
`travel_v3_summary.py` → `~/claude/work/logs/travel_v3_summary.log`). 18 błędów, 2 odrzucone off-pulse.

| | n | T/T_inc ≥ 5 (v2) | **T_cv ≥ 5, nb = 32** | nb = 64 | max(32, 64) |
|---|---|---|---|---|---|
| drift | 406 | 368 (91%) | **271 (67%)** | 238 (59%) | 278 (68%) |
| p3only | 107 | 79 (74%) | **12 (11%)** | 14 (13%) | 17 (16%) |

(Kolumna nb = 32 z fallbackiem na największy nb ≥ 16 dla 30 pulsarów; w wariantach nb = 64 i max bez.)
Tabela krzyżowa: p3only 67 detekcji tylko przez T, 0 tylko przez T_cv; drift 102 tylko T, 9 tylko T_cv.
z_cv rośnie z siłą T u dryferów (46% → 86% detekcji), u p3only nie (10–21% w każdym przedziale T).
T_adj: 57 detekcji drift, 0 p3only, wszystkie pokrywają się z T_cv.
ρ dla detekcji T_cv (P₂fit ≤ M/2): drift n = 167, mediana **1.000** (kw. 0.845–1.080); p3only n = 3, 0.217.

**Wniosek.** Po wyjęciu obciążenia od zmienności nieseparowalnej **trwałe uporządkowanie czasowe ma 2/3
dryferów i ~11% P3-only** (7–16% zależnie od podziału). Wcześniejsze 74% było w większości artefaktem
modelu zerowego. 12 P3-only z detekcją T_cv: J0837+0610, J1057-5226, J1048-5832, J1633-4453, J1701-3130,
J1810-5338, J1632-4621, J1121-5444, J1555-0515, J1816-5643, J1722-3207, J1130-6807.

**Otwarte.**
1. Dryfery z silnym T, ale bez T_cv (40 z T ≥ 100) mają dłuższe P₃ (mediana 12 wobec 6): bloki 32 P
   z lag ≤ 8 słabo pokrywają długi P₃. Kilka przechodzi przy nb = 64. Do sprawdzenia: lag_b niezależny od
   L÷4 albo bloki dobierane do P₃. Część (J0738-4042: T = 1.2·10⁵σ, z_cv 1.6; J1430-6623: z_cv < 0) może
   naprawdę nie mieć trwałego uporządkowania.
2. Reverser o losowych epizodach: T_cv go nie widzi, T_adj myli go z dudnieniem. Nieusuwalne bez modelu.
3. Test stabilności okna dla 12 P3-only z detekcją.

### 2026-09-29 — f_trav (siła trwałego dryfu) i wykresy P-Pdot

**Cel.** Wykres P-Pdot z wielkością mówiącą, jak silny jest dryf, niezależną od S/N (z_cv, T i
min(cons) od niego zależą, więc się nie nadają).

**f_trav** (`Travel.crossblock_test`, pola `f_trav`, `f_trav_err`, `f_trav_ksnr` w `travel_test`): ułamek mocy
fluktuacji w trwałym dryfie. Mapy blokowe dzielone przez liczbę par (L−τ)(M−Δ) i przez wariancję
fluktuacji z odjętym szumem off-pulse (k); moc trwała tylko z iloczynów między blokami (jak T_cv),
`f = sgn(m)√|m|/k`. Dla fali bieżącej komórka to 2 sin sin, RMS = 1. Błąd: jackknife po blokach.

**Walidacja** (`travel_ftrav_validate.jl` → `travel_ftrav_validate.csv`):
- mieszanina q·dryf + (1−q)·AM: f = 0.00 / 0.21 / 0.44 / 0.67 / 0.91 dla q = 0 … 1, stałe w zakresie
  36σ–34 000σ i prawie niezależne od nb; jitter, losowe podpulsy: |f| < 0.03 (stary T do 1207σ);
- 30 silnych dryferów + szum: f(k)/f(0) = 1.00 / 1.00 / 1.01 / 1.01 przy starej sile 2884σ → 42σ;
- **ograniczenie:** przy k/σ²_szumu ≲ 0.02 mianownik jest różnicą prawie równych wariancji i f jest
  zawyżone (syntetyk: 1.39 zamiast 0.91 przy k_snr ≈ 0.01).

**Skala na danych.** Mianownik zbiera całą zmienność impuls-do-impulsu (wahania energii, zmiany kształtu,
nulling — k_snr sięga setek), więc realne dryfery mają f ≈ 0.02–0.34, nie ~0.9. Przebieg v4 (`travel_batch_v4.csv`,
detekcje identyczne z v3): drift z detekcją T_cv f mediana **0.049** (kw. 0.030–0.090, max 0.344, n = 249);
p3only 0.018 (0.014–0.036, n = 11).

**Podejrzane: k_snr < 0 u 23 z 283 detekcji** (do −0.97): wariancja on-pulse mniejsza niż szacunek szumu
z off-pulse'u, czyli off-pulse zawyżony (emisja lub artefakty bazowe w paskach). T_cv tego nie dotyczy
(bez modelu szumu); f dla nich niezdefiniowane. Tego samego szacunku szumu używał stary null T —
do sprawdzenia, czy te 23 to te same pulsary, które mają dziwne kontrole off-pulse.

**Wykresy** (`Plot.ppdot_travel`, skrypt `~/claude/work/scripts/travel_ppdot.jl`):
`~/output/claude/ppdot_travel_ftrav.{pdf,png}` i `ppdot_travel_rho.{pdf,png}`. Kształt = etykieta Song+23,
kolor = f_trav (log 0.01–0.3) lub |ρ| (0–1.3; tylko P₂fit ≤ M/2), jasnoszare = trwałe uporządkowanie bez
wartości, puste = brak detekcji T_cv. `_ppdot` dostał trzy opcjonalne kwargi (`overlay`,
`population_color`, `legend_loc`); domyślne zachowanie bez zmian.

### 2026-09-29 (cd.) — zależność od Ė, τ_c, B

**Polecenia.** `~/claude/work/scripts/psrcat_dump.jl` → `~/claude/work/psrcat_ppdot.csv` (P, Ṗ z `input/psrcat.db`,
ta sama funkcja co wykresy); `travel_vs_ppdot.py` → `logs/travel_vs_ppdot.log`; `travel_vs_edot_fig.py` →
`logs/travel_vs_edot_fig.log`, `~/claude/work/figures/travel_vs_edot.png`. Ė = 4π²IṖ/P³ (I = 10⁴⁵), 513 pulsarów.

**1. Częstość trwałego uporządkowania maleje z Ė** (etykieta drift): 71% / 78% / 63% / 48% / 44% w przedziałach
Ė = 10²⁹⁻³¹ / 10³¹⁻³² / 10³²⁻³³ / 10³³⁻³⁴ / >10³⁴ erg/s; Mann-Whitney p = 9·10⁻⁵ (w drift), 3.5·10⁻⁸ (całość).
P3-only płasko ~10% (mała statystyka). **Nie wynika z S/N:** trend jest w obu połowach k_snr i starego T,
najsilniej w połowie o wyższym S/N (89% → 38%); logistyczna det ~ Ė + k_snr + T: Ė z = −4.6. Siła modulacji
k_snr nie koreluje z Ė (ρ_S = +0.02).

**2. f_trav maleje z Ė**: drift z detekcją n = 249, ρ_S = −0.30 [−0.41, −0.19]; po kontroli log P₃, log k_snr,
log T (proxy S/N) −0.31. Mediany w przedziałach Ė spadają ~0.06 → ~0.03. τ_c: +0.28 (to samo, bo τ_c i Ė są
silnie zależne w próbce). **B nie gra** (częściowa +0.03); po kontroli zostają P (+0.30) i Ė, Ṗ nie.
P₃ nie koreluje z Ė w tej próbce (+0.08), więc to nie efekt tłumienia przy P₃ → 2.

**3. |ρ| nie zależy od niczego** (|ρ_S| ≤ 0.02 dla Ė, τ_c, B, P, Ṗ; n = 167). Tam, gdzie dryf jest, ma charakter
sztywnej translacji niezależnie od Ė; z Ė zmienia się to, *jak często* jest trwały i *jaką część* zmienności
pulsów stanowi.

**Zastrzeżenia.** Selekcja próbki Song+23 (etykieta drift sama zależy od Ė: 90% → 62% udziału); ρ i f_trav tylko
dla detekcji; f_trav ma w mianowniku całą zmienność impuls-do-impulsu, więc „mniejsze f przy wysokim Ė” może
znaczyć zarówno słabszy dryf, jak i silniejszą niezależną zmienność (k_snr kontrolowane, ale nie jest tym samym).
Zgodne jakościowo z Basu et al. (2016) — dryf w pulsarach o Ė ≲ 10³² — ale tu jako trend ciągły, bez progu.

### 2026-09-30 — opis metody przepisany wokół T_cv / f_trav

`docs/travel_test_method.md`: nowy §0 (podsumowanie prostymi słowami + interpretacja fizyczna jako hipotezy),
§5 przepisany (5.1 dlaczego nie T, 5.2 T_cv, 5.3 walidacja, 5.4 f_trav, 5.5 wyniki v3/v4, 5.6 zależność od Ė,
5.7 kontrole), §4, §7.1, §7.4, §11, §12 (wycofany wniosek 8: „74% P3-only”), §13 zaktualizowane.
Wykresy w `docs/figures/` (~620 kB): `travel_null_calibration.png`, `travel_ftrav_validation.png`
(skrypt `~/claude/work/scripts/travel_doc_figures.py`), `ppdot_travel_{ftrav,rho}.png`, `travel_vs_edot.png`.

### 2026-10-01 — p3fold_coherent dla 10 pulsarów z najbardziej widocznym dryfem

**Kryteria wyboru** (ustalone przed wyborem; `~/claude/work/logs/p3fold_selection.log`): z_cv ≥ 5,
0.8 ≤ |ρ| ≤ 1.2, P₂fit/M ≤ 0.5, P₃ ≥ 2.5, N/P₃ ≥ 30; ranking wg **V = f_trav × k_snr** (moc trwałego dryfu
w pojedynczym impulsie względem szumu). Kryteria spełnia 88 pulsarów.

**Top 10** (V; matched-filter SNR z `p3fold_coherent`): J2139+2242 (37.9; 107), J0820-1350 (23.0; 90),
J1428-5530 (15.3; 35), J1932+1059 (14.7; 23), J1059-5742 (12.5; 23), J1041-1942 (8.6; 11), J0255-5304 (8.4; 33),
J1703-1846 (6.9; 21), J0934-5249 (5.6; 33), J0034-0721 (3.6; 20).

**Kod.** `SpaTs.p3fold_coherent` dostał kwargi `datafile` (dla `_16`: `pulsar_full_debase.txt`), `plotdir`,
`name_mod`; `p3_ybins` zaokrąglane (J1720-0212 ma 5.8 w params.json). Domyślne zachowanie bez zmian.
Skrypt `~/claude/work/scripts/p3fold_top10.jl`; wykresy `~/claude/work/figures/p3fold/<PSR>_coherent_p3fold_compare.{png,pdf}`,
SNR w `~/claude/work/p3fold_top10_snr.csv`.

### 2026-10-01 (cd.) — p3fold_coherent dla 10 P3-only z najbardziej prawdopodobnym dryfem

**Kryterium:** ranking wg z_cv (B = 32) wśród P3-only — prawdopodobieństwo trwałego uporządkowania; geometrii
(ρ, P₂) nie wymagano, bo dla P3-only jest przeważnie niemierzalna (`~/claude/work/logs/p3fold_p3only_selection.log`).
Wszystkie z katalogów `_16` (`pulsar_full_debase.txt`).

| PSR | z_cv | P₃ | f_trav | SNR złożenia |
|---|---|---|---|---|
| J0837+0610 | 12.3 | 2.17 | 0.044 | 8.6 |
| J1057-5226 | 11.7 | 8.51 | 0.004 | 9.6 |
| J1048-5832 | 9.9 | 17.4 | 0.011 | 14.5 |
| J1633-4453 | 9.2 | 16.5 | 0.016 | 3.5 |
| J1701-3130 | 9.1 | 28.0 | 0.038 | 2.8 |
| J1810-5338 | 8.2 | 4.81 | 0.016 | 3.2 |
| J1632-4621 | 7.6 | 15.1 | 0.012 | 3.9 |
| J1121-5444 | 7.1 | 31.6 | 0.037 | 3.3 |
| J1555-0515 | 6.4 | 2.33 | 0.035 | 5.0 |
| J1816-5643 | 6.2 | 19.7 | — (k_snr < 0) | 1.7 |

**Obserwacje (obejrzane 3 najsilniejsze).** Brak pochyłych pasów jak u dryferów. Widać modulację jasności
**poszczególnych składowych** z P₃ (J1057-5226: tylko środkowa składowa ~bin 95; J1048-5832: dwie składowe
~100 i ~125; J0837+0610: P₃ ≈ 2.17, blisko Nyquista). Trwałe uporządkowanie z T_cv pochodzi więc najpewniej ze
**stałego przesunięcia fazy modulacji między składowymi** (jedna składowa systematycznie wyprzedza drugą),
a nie z przesuwania się podpulsu — to degeneracja z §11 pkt 12 opisu metody. Do sprawdzenia: profil fazy
modulacji (faza LRFS przy f₃) w funkcji długości — liniowy = dryf, schodkowy = opóźnienie między składowymi.
Pozostałe 7 złożeń ma SNR 1.7–5, czyli są w dużej mierze szumowe.

Skrypt `p3fold_top10.jl` przyjmuje teraz argumenty `[lista] [katalog wyjściowy] [plik SNR]`; wykresy
`~/claude/work/figures/p3fold_p3only/`, SNR `~/claude/work/p3fold_p3only_snr.csv`.

### 2026-10-01 (cd. 2) — p3fold_coherent dla 10 P3-only bez trwałego ruchu (kontrola)

**Kryterium:** P3-only z wyraźną modulacją (k_snr ≥ 1, N/P₃ ≥ 30; 35 z 107), ranking wg **najniższego z_cv**
(B = 32). Sam najniższy z_cv wybrałby najsłabsze, zaszumione pulsary — tu brak ruchu jest pomiarem, nie
brakiem sygnału (`~/claude/work/logs/p3fold_p3only_worst_selection.log`).

| PSR | z_cv | stary T | P₃ | k_snr | SNR złożenia |
|---|---|---|---|---|---|
| J1401-6357 | −1.29 | **4652σ** | 2.21 | 283 | 27.0 |
| J1146-6030 | −1.25 | 616σ | 10.95 | 11.3 | 7.2 |
| J1901+1306 | −0.97 | −1σ | 11.83 | 1.4 | 1.3 |
| J1001-5939 | −0.33 | 103σ | 2.09 | 15.5 | 20.5 |
| J0955-5304 | −0.04 | 26σ | 3.54 | 1.6 | 3.0 |
| J1825+0004 | 0.03 | 150σ | 14.22 | 1.5 | 3.1 |
| J1143-5158 | 0.25 | 11σ | 5.03 | 1.6 | 2.4 |
| J1757-2421 | 0.30 | 0σ | 24.38 | 4.1 | 1.9 |
| J1603-2531 | 0.43 | **1709σ** | 48.62 | 27.8 | 23.9 |
| J2307+2225 | 0.63 | 18σ | 3.48 | 2.0 | 3.1 |

**Obserwacje (obejrzane 3 z wysokim SNR: J1401-6357, J1603-2531, J1001-5939).** Czysta modulacja jasności całej
składowej: poziome pasy, bez nachylenia i bez przesunięcia fazy między częściami profilu. To dokładnie
przypadek, który stary T fałszywie wykrywał na tysiącach σ (J1401: 4652σ, J1603: 1709σ), a T_cv poprawnie
zeruje. Kontrast z „najlepszymi” P3-only (wpis wyżej), gdzie modulowane są poszczególne składowe z możliwym
opóźnieniem fazy, jest subtelny wizualnie — rozstrzygnie profil fazy modulacji w funkcji długości.

Wykresy `~/claude/work/figures/p3fold_p3only_worst/`, SNR `~/claude/work/p3fold_p3only_worst_snr.csv`.

### 2026-10-01 (cd. 3) — J1825+0004: dryf widoczny, ale nietrwały (przypadek graniczny)

Użytkownik widzi dryf w złożeniu J1825+0004, który trafił do kontroli „bez ruchu” (z_cv = 0.03).
**Etykieta Song+23: P3-only** (`p3only_pulsars_P3.txt`: 14.1(9)). Skrypty `~/claude/work/scripts/j1825_{blocks,lrfs}.jl`.

- v4: T = 93σ, **T_inc = 150σ > T**, spójność blokowa −0.21 (4 bloki), z_cv ≈ 0 przy każdym B — sygnatura
  zmiany kierunku, na którą T_cv jest ślepy.
- **Jasność skacze ~2× ok. impulsu 715.** Moc map A w blokach po 715 jest ~10× większa, a znak wiodącego wzoru
  zmienia się co ~65 P (+0.93, −0.99, +0.47, −0.84, −0.17); przed 715 znak słaby, ale stały (+0.05…+0.13).
- **LRFS impulsy 1–700:** P₃ ≈ 14.9 w wiodącej składowej, faza przy f₃ monotoniczna wzdłuż składowej
  (−1.2 → +1.9 rad na binach 18–30, P₂ ≈ 24 biny) — **klasyczny dryf**. T_cv na 1–700: z = 3.5–4.0 (B = 8–32),
  poniżej progu (słaby sygnał, pulsar ciemny).
- **Impulsy 716–1040:** zmienność zdominowana przez wolne fluktuacje (65–325 P), P₃ ≈ 14 nie dominuje; T_cv ≈ 0.
- Normalizacja map per blok (równe wagi) nie zmienia wyniku (całość −0.8…+0.4; 1–700: +3.5…+4.0), więc to nie
  efekt ważenia, tylko rzeczywista niestacjonarność.

**Wniosek.** Dryf jest, ale tylko w części obserwacji i słaby; na całej obserwacji nie jest trwały. Oko ma rację,
T_cv też (w swoim sensie). Ograniczenie metody: **dryf obecny tylko w odcinku obserwacji jest rozmywany**.
Kandydat do batcha: T_cv osobno w połowach/odcinkach obserwacji jako diagnostyka niestacjonarności.

### 2026-10-01 (cd. 4) — nowa metoda, krok 1: sliding LRFS i ślad P₃(t) (`P3Track`)

Cel: zastąpić test travel. Krok 1 (polecenie użytkownika): sliding LRFS jak Fig. 4 w Szary+2022 (J1750-3503),
z jak najkrótszym oknem, do wykrywania odcinków ze stabilnym P₃.

Kod: `modules/p3track.jl` (`sliding_lrfs`, `p3_track`, `contrast_null`, `good_windows`), `Plot.sliding_lrfs`.
Skrypt `~/claude/work/scripts/p3track_test.jl`, log `~/claude/work/logs/p3track_test.log`,
wykresy `~/claude/work/figures/p3track/<PSR>_sliding_lrfs_L<L>.png` (L = 16…256, stride 1).

- Okno z taperem Hann (periodyczny) + zero-padding ×8; P₃ z dopasowania Gaussa do najwyższego **lokalnego**
  maksimum wewnątrz (2/L, 0.5) — argmax globalny łapał czerwony szum od nulli (J0034-0721) i wolnych zmian.
- Pik bliżej niż 1/L od fmin = 2/L → `edge` (wyciek DC). Mierzalne praktycznie **P₃ ≲ L/3**.
- Jakość okna: S/N względem off-pulse'u bezużyteczne (J0820: mediana 3000–12000). Zamiast tego **kontrast**
  (pik / mediana widma) wobec progu 99% z **lokalnego tasowania impulsów w obrębie okna**. Tasowanie globalne
  zawyżało próg w ciemnych odcinkach (J1825+0004, L = 64: 3.79 globalnie vs 2.38 z impulsów 1–700).
- Wyniki (udział dobrych okien, mediana P₃): J0820-1350 L=16 76% (4.79), L=32 95%; J0151-0635 L=64 94% (14.25);
  J1825+0004 L=64 53% (14.5) — dobre okna dokładnie w 1–680, część jasna po ~715 odrzucona;
  J0034-0721 L=32 25% (6.64), tylko w seriach bez nulli; J1750-3503 L=128 36%, L=256 65% (48.7).
- J1750-3503 L=128 jakościowo zgodne z Fig. 4b (minimum P₃ ≈ 45 ok. impulsu 230, ~37–42 w 500–600),
  ale P₃ ≈ 40–60 przy L=128 to tylko 2–3 cykle — na granicy.

Otwarte: kryterium „stabilnego P₃” (rozrzut estymatora przy krótkim L, np. J0820 L=16: P₃ 4.4–5.3 — szum
estymatora czy realna zmienność?); wybór L względem P₃.

### 2026-10-01 (cd. 5) — P3Track krok 2: odcinki ciągłego P₃ i grupy do foldowania

Ustalenia z użytkownikiem: L = max(16, 4·P₃) (`window_length`); monotoniczna zmiana P₃ jest OK (kompensacja
przy foldowaniu); impulsy grupować po P₃ i foldować grupy osobno.

Kod: `p3_segments` (ciągłość śladu: sąsiednie dobre okna, przerwa ≤ L/2, |Δf₃| ≤ 0.25/L; zakres impulsów =
środki okien ± L/2, przycięte w połowie odstępu do sąsiedniego odcinka), `p3_groups` (łączenie po medianie f₃
z tolerancją 1/L ≈ rozdzielczość okna Hann), `merge_sections` (stykające się odcinki tej samej grupy).
Skrypt `~/claude/work/scripts/p3track_segments.jl`, log `~/claude/work/logs/p3track_segments.log`,
wykresy `~/claude/work/figures/p3track/*_seg.png`.

| pulsar | L | impulsy w odcinkach | grupy (P₃, impulsy) |
|---|---|---|---|
| J0820-1350 | 19 | 98% | 1 odcinek, P₃ 4.77 |
| J0151-0635 | 58 | 100% | 14.3 (999); 7.46 (40) = 2. harmoniczna |
| J1825+0004 | 57 | 65% | 14.5 (1–696); część jasna po ~715 odrzucona |
| J0034-0721 | 48 | 66% | 6.52 (7 serii między nullami, 599); 10.1 (940–1031, tryb A/przejście) |
| J0034-0721 | 26 | 48% | 6.6 (10 serii); tryb A (P₃ ~ 12) poza zasięgiem L/3 |
| J1750-3503 | 196 | 87% | 2 grupy, P₃ 38–68; przy N = 1031 i P₃ ~ 50 metoda na granicy |

- tol = 0.5/L dzieliło szum estymatora (J1825: odcinek 33 P przy 11.7 obok 14.5, Δf = 0.95/L) → tol = 1/L.
- Bez scalania J0820 przy L = 19 rozpadało się na 13–20 odcinków o tym samym P₃ (rozrzut estymatora > 0.25/L).

Otwarte: rozpoznawanie harmonicznych (grupa z f₃ ≈ 2·f₃ grupy głównej to ten sam reżim); L wg P₃ z params
odcina dłuższe P₃ innych trybów (J0034 tryb A przy L = 26); minimalna liczba impulsów grupy do foldowania.

### 2026-10-01 (cd. 6) — P3Track krok 3: harmoniczne, min. 5 cykli, fold z kompensacją zmiennego P₃

Decyzje użytkownika: harmoniczne dołączać do grupy głównej; grupa ≥ 5·P₃ impulsów; L z params.json (test L
z najdłuższego P₃ w średnim LRFS odłożony); kompensacja zmiennego P₃ jako nowa funkcja, wszystko w jednym
pliku → `Plot.sliding_lrfs` przeniesione do `P3Track.plot_track`, usunięte z `plot.jl`.

Nowe w `modules/p3track.jl`: `harmonic_groups`, `fundamental_track`, `select_groups`, `demodulate`,
`align_phases`, `phase_fold`, `constant_fold`, `analyse` (cały łańcuch), `plot_folds`.
Skrypt `~/claude/work/scripts/p3track_fold.jl`, log `~/claude/work/logs/p3track_fold.log`,
wykresy `~/claude/work/figures/p3track/<PSR>_p3fold_groups.png`, `<PSR>_sliding_lrfs_L<L>_groups.png`.

**Fold:** faza modulacji mierzona w każdym impulsie, nie całkowana z P₃: demodulacja zespolona przy lokalnym
f₃(n) (interpolacja śladu), jądro odniesione do impulsu n → arg Z(n,φ) = θ(n) − ψ(φ); θ(n) względem wspólnego
szablonu grupy T(φ) (iteracja). Szablon niesie kształt pasma dryfu, więc odcinki rozdzielone nullami zgrywają
się automatycznie. **Leave-one-out** (impuls n z wagą 0 we własnym oknie) jest konieczny: bez niego kontrola
z tasowaniem dawała głębokość 0.15–0.18 (J0034) przy 0.06 dla stałego P₃ — szum impulsu ustawiał jego fazę.

| pulsar | L | impulsy | depth zmienne P₃ | kontrola (tasowanie) | depth stałe P₃ | koherencja |
|---|---|---|---|---|---|---|
| J0034-0721 | 26 | 501 (10 serii) | **0.206** | 0.069–0.090 | 0.062 | 0.76 |
| J0151-0635 | 58 | 1039 (z harmon. 400–439) | **0.172** | 0.044–0.052 | 0.065 | 0.90 |
| J0820-1350 | 19 | 1034 | **0.142** | 0.018–0.021 | 0.071 | 0.93 |
| J1825+0004 | 57 | 681 | **0.086** | 0.037–0.039 | 0.043 | 0.79 |
| J1750-3503 | 196 | 897 | **0.205** | 0.150–0.154 | 0.160 | 0.67 |

Foldy pokazują pasma dryfu tam, gdzie stały P₃ daje rozmazany profil (J0034, J0151, J1825). Faza względem stałego
P₃ wędruje o 1–2 cykle w obrębie obserwacji (J0151, J1825, J1750).

Ograniczenie: jeden szablon na grupę → u reversera (J1750-3503) epizody o przeciwnym kierunku dryfu wchodzą do
jednego folda z kształtem dominującego kierunku. Głębokość modulacji to miara porównawcza, nie istotność.

### 2026-10-01 (cd. 7) — P3Track: kontrola na 5 P3-only bez oczekiwanego dryfu; faza szablonu jako dyskryminator

Kontrola: J1401-6357, J1603-2531, J1001-5939 (obejrzane wcześniej: czysta AM), J1146-6030, J2307+2225 (najniższe
z_cv przy k_snr ≥ 1, wpis cd. 2). Skrypt `~/claude/work/scripts/p3track_control.jl`, log
`~/claude/work/logs/p3track_control.log`, wykresy `~/claude/work/figures/p3track/<PSR>_{p3fold_groups,sliding_lrfs_L<L>_groups}.png`.

Zmiany w kodzie: (1) usunięta górna osłona krawędzi w `feature_peak` (P₃ ≈ 2.1 → f ≈ 0.48 była odrzucana w całości;
wystarcza wymóg lokalnego maksimum); (2) `template_phase`: ψ(φ) = −arg T szablonu grupy, nachylenie liczone
**osobno w każdej składowej** (offset fazy między rozdzielonymi składowymi jest określony tylko mod 1 cykl —
jedno dopasowanie przez przerwę dawało dla J0151 0.06 zamiast 0.86); Δψ = Σ|Δψ_składowej|; nowy wiersz w `plot_folds`.

| etykieta | PSR | L | w grupach | grupa P₃ (impulsy) | depth zm. / tasow. / stałe | **Δψ [cyk]** (składowe) |
|---|---|---|---|---|---|---|
| P3-only | J1401-6357 | 16 | 6% | 2.28 (58) | 0.155 / 0.04–0.08 / 0.116 | **0.00** |
| P3-only | J1603-2531 | 194 | 95% | 34.1 (1003), 51.7 (832), 13.0 (181) | 0.08–0.10 / ~0.04 / ~0.05 | **0.01, 0.01, 0.01** |
| P3-only | J1001-5939 | 16 | 0% | – | – | – |
| P3-only | J1146-6030 | 44 | 0% | – | – | – |
| P3-only | J2307+2225 | 16 | 0% (grupa 20 P < 5·P₃) | – | – | – |
| drift | J0034-0721 | 26 | 48% | 6.64 (501) | 0.206 / 0.07–0.09 / 0.062 | **1.53** |
| drift | J0151-0635 | 58 | 100% | 14.27 (1039) | 0.172 / 0.05 / 0.065 | **0.86** (−0.69, −0.17) |
| drift | J0820-1350 | 19 | 98% | 4.77 (1034) | 0.142 / 0.02 / 0.071 | **1.94** |
| drift | J1825+0004 | 57 | 65% | 14.67 (681) | 0.086 / 0.04 / 0.043 | **0.23** (+0.03, −0.20) |
| drift | J1750-3503 | 196 | 87% | 44.9 (897) | 0.205 / 0.15 / 0.160 | **2.03** |

- Głębokość modulacji rośnie po kompensacji także dla P3-only (J1603: P₃ wędruje 13–52 w trzech grupach) — sama w sobie
  nie odróżnia dryfu od AM. Odróżnia faza szablonu: P3-only ≤ 0.01 cyklu, dryfery 0.23–2.03.
- J1001, J1146, J2307: brak okien ze stabilną cechą P₃ (good 3, 10, 12 z ~1000) — modulacja niekoherentna na skali 4·P₃,
  w średnim LRFS cecha szeroka. Metoda zwraca „brak stabilnego P₃”, nie werdykt.

Otwarte: skala szumu Δψ (kontrola tasowaniem lub z błędu ψ na bin) — próg dryf/AM jest na razie tylko empiryczny
(przerwa 0.01 → 0.23); J1401 tylko 6% impulsów w grupach.

### 2026-10-01 (cd. 8) — dokument metody P3Track

`docs/p3track_method.md`: opis metody (kroki 1–5), wyniki na 5 dryferach + 5 P3-only, sprawy do rozstrzygnięcia
(§8: kalibracja progu Δψ, kategoria „brak stabilnego P₃”, L przy wielu trybach, reverserzy, J1401, maska składowych,
harmoniczne, Δψ a P₂, znak dryfu, batch), historia poprawek, użycie. Wykresy w `docs/figures/p3track_*.png`.

### 2026-10-01 (cd. 9) — J1825+0004: drugi reżim P₃ po zmianie modu (uwaga użytkownika)

Użytkownik widzi w sLRFS J1825+0004 dwie grupy P₃, w tym stabilne P₃ po zmianie modu (~715). Przy L = 57 (z P₃ = 14.2
z params) cecha drugiego reżimu jest poza zasięgiem (mierzalne P₃ ≲ L/3 ≈ 19). Skrypt
`~/claude/work/scripts/j1825_mode2.jl`, log `~/claude/work/logs/j1825_mode2.log`,
wykresy `~/claude/work/figures/p3track/J1825+0004_sliding_lrfs_L{96,128,160}_groups.png`, `J1825+0004_L{96,128,160}_p3fold_groups.png`.

- Średnie widmo (okna 128): 1–700 jedna cecha P₃ = 14.6 (kontrast 3.0); 716–1040 kilka: 21.8 (2.3), 10.1, 12.6, 6.7.
- **L = 160: grupa 2, P₃ ≈ 38 (55 → 34 monotonicznie), impulsy 744–972 (229 P, ~6 cykli), Δψ = 0.06** — złożenie
  pokazuje poziome pasy na całym profilu: modulacja amplitudowa, nie dryf. Grupa dryfu P₃ ≈ 14.3 ma Δψ = 0.22.
- L = 128: ta sama grupa odrzucona (182 P < 5·38.7); L = 96: P₃ ≈ 25 (203 P), zmieszana z odcinkiem 220–286.
- Dłuższe L zaostrza tolerancję grup (1/L) i dzieli reżim 1 (L = 160: 16.5 w 1–223 i 14.25 w 224–737).

Wniosek: J1825+0004 = mod dryfu (P₃ ≈ 14.5, impulsy 1–700) + mod jasny z wolną AM (P₃ ≈ 35–55). Pojedyncze L dobrane
do P₃ z params nie widzi drugiego reżimu — potwierdza odłożony problem długości okna przy wielu trybach.

### 2026-10-01 (cd. 10) — P3Track: drugie przejście dla reżimów o długim P₃

Zgoda użytkownika na drugie przejście po impulsach spoza grup (`long_p3_pass`, domyślnie w `analyse`): sonda
(średnie widmo wolnych fragmentów, maksimum o największej wybitności w 3/Lp ≤ f < 3/L₁), drabinka okien
min(4·P₃′, wolne) + L₁·{2,3,4,6,8}, poszukiwanie tylko f < 3/L₁, odcinki przycięte do wolnych impulsów; wybór
okna z największą liczbą impulsów w grupach. Kontrast zawsze względem mediany z f ≥ fmin (w wąskim zakresie
niskich f mediana siedziała na czerwonym kontinuum i nic nie przechodziło). Log `~/claude/work/logs/p3track_control_pass2.log`,
wykresy `*_pass2_p3fold_groups.png`, `*_groups_pass2.png`.

- J1825+0004: **P₃ ≈ 36.7, impulsy 729–987, Δψ = 0.08** (L = 171) — drugi reżim (AM) odzyskany.
- J1146-6030 (P3-only, wcześniej brak grup): 20.8 (159–425, Δψ 0.02) i 16.2 (802–1007, Δψ 0.44 przy L = 132, ale
  0.04 przy L = 84; złożenie prawie poziome, modulacja słaba) — Δψ niestabilne przy słabej modulacji.
- J0034-0721 tryb A nie odzyskany (impulsy zajęte w pierwszym przejściu); J1401, J1001, J2307 bez fałszywych grup.
- Pierwsze przejście bez zmian. Dokument `docs/p3track_method.md` uzupełniony (§4a, §7, §8, §9).

### 2026-10-01 (cd. 11) — wykres zbiorczy obu przejść (`plot_summary`)

Użytkownik na wykresie drugiego przejścia J1825+0004 widział tylko mod końcowy — mod dryfu (P₃ ≈ 14.7, 1–696) jest
w pierwszym przejściu (L = 57), a drugie z założenia szuka tylko w wolnych impulsach i tylko P₃ > L₁/3.
`P3Track.plot_summary(data, res, outdir)`: stos impulsów z paskami grup obu przejść + ślady P₃ wszystkich grup
z etykietą (przejście, L, P₃, impulsy, Δψ). `~/claude/work/figures/p3track/<PSR>_p3track_summary.png`,
skrypt `~/claude/work/scripts/p3track_summary_one.jl <katalog> <plik>`; dodany też do `p3track_control.jl`.

### 2026-10-01 (cd. 12) — test harmonicznej przed scaleniem (`harmonic_test`)

Pytanie użytkownika o krótsze P₃: pierwsze przejście je obejmuje (zakres do f = 0.5; J1603 przy L = 194 dało grupę
P₃ ≈ 13). Ryzyko: osobny mod o P₃ ≈ ½ scalany jako 2. harmoniczna. Zgoda na test przed scaleniem.

Test: impulsy kandydata składane przy f₃/2, składowa Fouriera h = 1 wzdłuż fazy P₃ vs 20 tasowań. Pierwsza wersja
(całkowita głębokość złożenia) błędna: sygnał przy f złożony przy f/2 daje dwa cykle na fold, więc „modulacja”
jest zawsze. Werdykty: harmonic / separate (≥ 10 cykli fundamentalnej) / inconclusive (odcinki usuwane).
Skrypt `~/claude/work/scripts/p3track_harmonic_test.jl`, log `~/claude/work/logs/p3track_harmonic_test.log`.

| przypadek | werdykt | h1(f/2) | max tasowań |
|---|---|---|---|
| syntetyk osobny mod P₃ = 4, szum 0.6 / 1.2 | separate / separate | 0.036 / 0.067 | 0.042 / 0.079 |
| syntetyk harmoniczna, szum 0.6 / 1.2 | harmonic / harmonic | 0.065 / 0.084 | 0.047 / 0.080 |
| syntetyk harmoniczna, szum 2.0 (42 P) | inconclusive | 0.220 | 0.250 |
| J0151-0635, 400–439 (40 P, 2.7 cyklu) | inconclusive | 0.046 | 0.065 |

Próg 10 cykli z dwóch syntetyków — do weryfikacji. Dokument metody uzupełniony (§4, §8.7, §9).

### 2026-10-01 (cd. 13) — przeliczenie 10 pulsarów pełną metodą (oba przejścia + test harmonicznej)

`~/claude/work/scripts/p3track_control.jl`, log `~/claude/work/logs/p3track_control_all.log`, wykresy zbiorcze
`~/claude/work/figures/p3track/<PSR>_p3track_summary.png`. Wyniki identyczne z cd. 7/10/12 (jedyna zmiana: J0151
kandydat 2:1 → inconclusive, grupa 1037 P). Dryfery Δψ 0.23–2.03 (+ mod AM J1825 po ~715: 0.08); P3-only Δψ
0.00–0.02, wyjątek J1146-6030 802–1007: 0.44 przy słabej modulacji (niestabilne względem L); J1001, J2307 bez grup,
J1401 tylko 6% impulsów.

### 2026-10-01 (cd. 14) — kalibracja Δψ (`template_significance`)

Bootstrap blokowy (bloki L/2) po impulsach grupy → σ Δψ każdego fragmentu, χ² → z; werdykt drift (z ≥ 5 i Δψ ≥ 0.1),
am (Σ(|Δψ|+2σ) < 0.1), inconclusive (reszta lub < 5 bloków). Kalibracja na syntetykach (`p3track_dpsi_calib.jl`,
logi `p3track_dpsi_calib*.log`) wymusiła zmianę miary: pierwsza wersja dawała 13/20 fałszywych dryfów w AM
z przeciwfazą (krótkie fragmenty w strefie znoszenia, szum na krawędziach) → gradient z Σ T*_j T_{j+1} bez rozwijania,
fragmenty ≥ 5 binów i ≥ 5% mocy, maska 0.15 (0.2 odcinało ogon J1825). Wynik: 0/80 fałszywych dryfów, dryf 10/10.

10 pulsarów (`p3track_control_dpsi.log`): J0034, J0151, J0820 (z = 37), J1750 (15.4) → drift; J1603 g1/g2, J1401 → am;
**J1825+0004 mod dryfu → inconclusive** (−0.10 ± 0.04, z = 2.5: faza płaska w głównej części składowej, zmienia się
tylko w słabym ogonie); J1825 mod 2, J1603 g3, J1146 obie → inconclusive (za mało bloków / za mały Δψ).
Dokument metody: §6.1, §7, §8.1, §9, §10.

### 2026-10-01 (cd. 15) — J1825+0004: dlaczego mod dryfu wychodzi „inconclusive”

Skrypt `~/claude/work/scripts/j1825_inspect.jl`, log `~/claude/work/logs/j1825_inspect.log`,
wykres `~/claude/work/figures/p3track/J1825+0004_inspect.png`. (Uwaga techniczna: `using PyPlot` przed
`include("modules/data.jl")` wywala ładowanie Glib_jll — systemowa libmount bez MOUNT_2_40; PyPlot ładować po include.)

Szablon grupy (681 P, P₃ = 14.67), faza ψ = −arg T z σ z bootstrapu:
- **172.3–175.1° (szczyt składowej, |T| 0.23–1.0): ψ płaskie, −108° … −85° … −121°, σ 1–4°** — modulacja w fazie (AM).
- **175.1–177.5° (opadające zbocze, |T| 0.46 → 0.11): ψ spada o ~250° (≈ 0.7 cyklu)**: −121 → −173 (bin 25, |T| = 0.18,
  lokalne minimum) → +93 (= −267) → 72 → 46 → 18 → 17 → 7 (= −353); σ 6–13°.
- Kształt **ten sam we wszystkich czterech ćwiartkach czasu** (1–170, 171–355, 356–526, 527–696) — trwały, nie szum.
- Fold − profil: główna łata 172–176° pozioma, łata na zboczu 176–178° opóźniona o ~0.3–0.4 cyklu → wrażenie
  nachylonego pasma w złożeniu.
- Klasyczny LRFS (jeden bin FFT, 1–700): faza zaszumiona, P₃ wędruje, mało informatywny.

Dlaczego Δψ = −0.10 ± 0.04: gradient z Σ T*_j T_{j+1} jest ważony amplitudą w całym fragmencie 17:30, a płaska,
jasna część dominuje wagę; stromy spadek fazy na słabszym zboczu się rozcieńcza. Miara zakłada gradient jednorodny
w składowej. Pytanie definicyjne do użytkownika: czy opóźnienie fazy na zboczu (skok przy minimum |T| + łagodny spadek
w ogonie) to dryf, czy AM z opóźnioną częścią zbocza (degeneracja „dryf” vs „kontinuum składowych opóźnionych w czasie”).

### 2026-10-01 (cd. 16) — kategoria „dryf częściowy” (`partial`)

Decyzja użytkownika: definicja b (dryf = wzór przesuwa się przez dominującą część emisji) + kategoria dryfu częściowego.
`template_significance`: gdy nie `drift`, szukane okno 5 binów w masce z |Δψ| ≥ 0.1, z ≥ 5 (bootstrap), zmianą rozłożoną
(maks. przyrost ≤ 50%), bez głębokiego minimum |T| (`deep_dip`) i z σψ ≤ 20° na bin. Kalibracja (`p3track_dpsi_calib.jl`,
log `p3track_dpsi_calib_partial.log`): bez dwóch ostatnich warunków 2/20 fałszywych `partial` w AM z przeciwfazą przy
szumie 1.5 (najpierw rozmyty skok przy minimum, potem biny szumu na krawędziach); po nich 0/80.

10 pulsarów (`p3track_control_partial.log`): **J1825+0004 mod dryfu → partial** (okno zbocza 175.8–177.2° bez skoku:
Δψ = −0.29, z = 7.9; odporne na okno 4–5 i udział ≤ 0.4–0.5, `j1825_partial_check.jl`); reszta bez zmian — dryfery
drift (lokalne okna też znalezione), P3-only am/inconclusive, żadnego `partial`.

### 2026-10-01 (cd. 17) — batch P3Track: pilot 20 pulsarów i porównanie zap / bez zap

Skrypt `~/claude/work/scripts/p3track_batch.jl` (wznawialny, `--part k/n`, `--limit N`, `--psrs`, `--tag`, `--nozap`),
uruchomienie `p3track_pilot_run.sh`. Wyniki: `~/output/claude/p3track_batch/p3track_{pilot,zaps_on,zaps_off}_part1of1.csv`,
wykresy `~/output/claude/p3track_batch/figures_<tag>/` i `~/claude/work/figures/p3track_batch_<tag>/`, logi
`~/claude/work/logs/p3track_batch_*.log`.

Uwaga: pierwsze, przerwane przez użytkownika uruchomienie pilota zdążyło zapisać 8 pulsarów (stary format kolumn) —
usunięte razem z katalogami wykresów, wszystko przeliczone od nowa.

**Zapy**: w całej próbce (521 z danymi) zapy w params.json ma tylko 4 pulsary, wszystkie drift; w pilocie (20) żaden.
J1524-5706 i J1843-0211 mają te impulsy wyzerowane już w archiwum (264 i 46 wierszy zer) — `Data.zap!` nic nie zmienia;
J1915+0752 też 0 nowych; J2139+2242: 91 impulsów. Wyniki z zapami i bez: **identyczne** (J2139: 906 vs 907 impulsów
w grupie, Δψ 1.74, drift w obu). Błąd przy tym: okno z samych zer → kontrast 0/0 = NaN → `quantile` w `contrast_null`
się wywalał (J1524, J1843 w obu wariantach, J2139 z zapami) — poprawione w `feature_peak` (okno puste: brak cechy).

**Pilot** (10 drift + 10 P3-only, pierwsze z listy z danymi), 10.2 min (średnio 31 s, max 157 s J0659+1414):
- drift: J0034, J0108, J0151, J0255, J0421 → drift; J0134 → partial; J0211, J0304, J0401 → inconclusive; J0152 → nogroup.
- P3-only: 0 × drift, 0 × am; 6 × inconclusive (J0601, J0629, J0659, J0709, J0737, J0831, J0849), 4 × nogroup (J0836, J0837,
  J0855 — P₃ ≈ 2.0–2.3, blisko Nyquista).
- Najczęstsza przyczyna inconclusive: < 5 bloków (krótkie grupy, zwłaszcza drugie przejście z długim L). Kilka takich
  grup P3-only ma formalnie duże z przy Δψ 0.3–0.45 (J0601, J0659, J0709) — reguła minblocks słusznie blokuje werdykt.
- `am` nie wychodzi ani razu: górna granica Σ(|Δψ|+2σ) < 0.1 jest za ostra dla realnych danych (np. J0709 g1:
  Δψ = 0.04, granica 0.17; J0849: 0.04, 0.13).

### 2026-10-01 (cd. 18) — nowe kryterium `am`, bez zapowania

Decyzje użytkownika: zapowanie pominięte (usunięte z `p3track_batch.jl`); poprawić `am` przed pełnym batchem.
`am` = z < 3, brak `partial`, Σ|Δψ_frag| + 2·√(Σσ²) < 0.25 cyklu (wcześniej Σ(|Δψ|+2σ) < 0.1 — ani razu w pilocie).
Kalibracja (`p3track_dpsi_calib_am.log`, dodany syntetyk wolnego dryfu P₂ = 120 binów): AM wspólna faza am 20/20 (szum 0.6),
13/20 (1.5); AM przeciwfaza 19/20 (0.6), 0/20 (1.5, inconclusive); **fałszywe am w klasach z dryfem: 0/50**.
Pilot ponownie (`p3track_pilot_am_part1of1.csv`, 10.1 min): nowe `am` — J0659+1414 (grupa drugiego przejścia P₃ ≈ 48),
J0709-5923 (P₃ 25.8), J0849-6322 (P₃ 8.0), J0304+1932 (drift wg Song, grupa drugiego przejścia P₃ ≈ 28 — mod AM obok
nierozstrzygniętej grupy P₃ ≈ 6.3, jak J1825). Reszta bez zmian; dominująca przyczyna inconclusive: < 5 bloków.

### 2026-10-01 (cd. 19) — pełny batch P3Track v1 (533 pulsary)

`~/claude/work/scripts/p3track_batch_run.sh` (8 procesów, `p3track_batch.jl --part k/8 --tag v1`, bez zapowania),
32–50 min na część (CPU łącznie 5.2 h, mediana 26 s/pulsar, max 1177 s). Wyniki: `~/output/claude/p3track_batch/p3track_v1.csv`
(907 grup), `p3track_v1_pulsars.csv` (werdykt pulsara = najsilniejszy z grup: drift > partial > am > inconclusive),
wykresy `~/output/claude/p3track_batch/figures_v1/` i `~/claude/work/figures/p3track_batch_v1/` (2316 plików).
Logi `~/claude/work/logs/p3track_batch_v1_part*.log`. Błędy: 12 × „brak danych” (6 drift, 6 P3-only).

| etykieta Song+23 | n (analiza) | drift | partial | am | inconclusive | brak grup |
|---|---|---|---|---|---|---|
| drift | 412 | 144 (35%) | 18 (4%) | 39 (9%) | 173 (42%) | 38 (9%) |
| P3-only | 109 | 3 (3%) | 2 (2%) | 28 (26%) | 55 (50%) | 21 (19%) |

- P3-only z drift/partial: J1810-5338, J1543+0929, J1057-5226 (drift), J1016-5345, J1825+0004 (partial).
- 8 pulsarów ma grupę drift/partial i grupę am (mody): J1651-5222, J0905-4536, J1750-3157, J1933+1304, J1735-0724,
  J1705-3423, J1922+1733, J1648-6044.
- **Inconclusive: 504 z 575 takich grup to < 5 bloków** (grupa krótka względem L/2) — główne ograniczenie czułości.
- Pokrycie (ułamek impulsów w grupach): drift mediana 0.41 (kw. 0.14–0.68), P3-only 0.30 (0.05–0.70).
- vs T_cv (z_cv(B=32) ≥ 5), etykieta drift: T_cv+ → drift 122, partial 12, am 18, inconclusive 98, brak 13;
  T_cv− → drift 14, partial 6, am 21, inconclusive 62, brak 18. P3-only: T_cv+ (12) → drift 2, am 1, inconcl. 8, brak 1.
- 49 grup `am` u pulsarów z etykietą drift ma Δψ ≤ 0.08 — do obejrzenia: mod AM, czy dryf z P₂ ≫ W (Δψ ≈ W/P₂ małe,
  kryterium am tego nie odróżnia).

### 2026-10-01 (cd. 20) — przegląd wątpliwych werdyktów batcha v1

Skrypt `~/claude/work/scripts/p3track_analysis/am_drift.py` (tabela grup am u dryferów z P₂ i W3s z travel v4); wykresy z `~/claude/work/figures/p3track_batch_v1/`.

Trzy typy problemów:
1. **Werdykt pulsara z „cudzej” grupy.** Wiele grup `am` u dryferów to drugie przejście z P₃ ≫ P₃ katalogowego
   (J1742-4616: 7 → 48, cała obserwacja; J1913+0936: 4 → 56; J1055-6905: 2.4 → 27; J1848+0604: 2.1 → 33) — wolna
   quasi-okresowa modulacja (burst/null), a cecha dryfu z katalogu nie dostała grupy. Pulsar dostaje `am` z tej grupy.
2. **Dryf w słabej składowej przy płaskiej fazie dominującej** — wg definicji b to raczej `partial`, a dostaje `drift`,
   bo Δψ sumuje składowe bez wag mocy: J1057-5226 (P3-only; 4924 P; składowa 101° płaska, słabsza 86–91°: −0.19 ± 0.03),
   J1543+0929 (P3-only; dominująca −0.09 ± 0.02, słaba +0.16 ± 0.04).
3. **`am` przy granicznym gradiencie.** J1511-5414 (drift wg Song): ψ rośnie wyraźnie przez składową,
   +0.09 ± 0.03 (z = 2.7 < 3) → `am`. Wolny dryf / P₂ ≫ W wygląda tak samo.
J1946-2913 (drift wg Song, T_cv z = 10.9): główna składowa płaska (±10°), słaba składowa ~185° poniżej maski —
uporządkowanie widziane przez T_cv może pochodzić z przesunięcia fazy między składowymi (degeneracja, jak J1825).

### 2026-10-01 (cd. 21) — batch v2: reguła mocy (B), am przy z < 2 (C), werdykt z P₃ katalogowego (A)

Batch v2 (8 procesów, 31–49 min, 12 × brak danych): `~/output/claude/p3track_batch/p3track_v2.csv`, werdykty pulsarów
`p3track_v2_pulsars.csv` (`batch_summary_v2.py`: werdykt z grup o P₃ ±30% katalogowego; „nocat” = grupy tylko przy innym
P₃; inne mody osobno). Kalibracja przed batchem: `p3track_dpsi_calib_v2.log` (fałszywe drift/am 0).

| etykieta | n | drift | partial | am | inconclusive | nocat | brak grup |
|---|---|---|---|---|---|---|---|
| drift | 412 | 127 (31%) | 26 (6%) | 16 (4%) | 111 (27%) | 94 (23%) | 38 (9%) |
| P3-only | 109 | 1 (1%) | 4 (4%) | 24 (22%) | 34 (31%) | 25 (23%) | 21 (19%) |

(v1 z tą samą regułą A: drift 135/18/21/106/94/38, P3-only 3/2/26/32/25/21.)
- Przypadki z przeglądu: J1057-5226, J1543+0929 → partial; J1511-5414 → inconclusive; J1742-4616, J1055-6905 → nocat
  (am tylko jako inny mod); J1946-2913 → am (płaska główna składowa). P3-only z drift: tylko J1810-5338 (73 P, z = 5.4).
- **nocat = 23% w obu etykietach**: grupy istnieją, ale przy P₃ innym niż katalogowe. W pierwszym przejściu mediana
  P₃grp/P₃kat = 0.50 → zwykle dominuje 2. harmoniczna, a fundamentalna nie ma własnej grupy, więc test harmonicznej
  (wymaga obu grup) nie ma czego łączyć. W drugim przejściu: wolna modulacja (×1.4 … ×16).
- Inconclusive (grupy): 504 × < 5 bloków, 46 × 2 ≤ z < 5, 31 × granica ≥ 0.25, 3 × z ≥ 5 przy za małym Δψ.

### 2026-10-02 — przegląd losowych przykładów v2; bi-drift; werdykt krótkich grup z folda

Użytkownik: J1537-4912 dryf (bi-drift — sprawdzić), J1907+0740 AM, J1528-4109 wygląda na dryf (inconclusive przez < 5 bloków).
- **J1537-4912 bi-drift potwierdzony** (`j1537_bidrift.jl`): składowa 163–177° (89% mocy) −0.19 ± 0.02, 183–193° (11%)
  +0.17 ± 0.03; w ćwiartkach czasu −0.24/+0.22, −0.14/+0.14, −0.24/+0.33, −0.10/+0.11. Flaga `bidrift` w kodzie.
- **Krótkie grupy**: nowy estymator z folda (bootstrap impuls po impulsie), używany gdy < 5 bloków. Kalibracja
  `p3track_short_calib.jl` (log `p3track_short_calib.log`): AM bez fałszywego dryfu, krótki dryf czulej. Dane
  (`p3track_short_real.jl`, log `p3track_short_real.log`): J1528 → drift; krótkie grupy P3-only (J0601, J0659, J0709, J0849),
  które miały z ≈ 5–13 z bootstrapu blokowego, → inconclusive (z_fold ≤ 1.9) — tamte z były artefaktem 2–3 bloków.
- Batch: kolumny verdict_block, verdict_fold, z_fold, dpsi_fold, verdict_src, bidrift.

### 2026-10-02 (cd.) — co kryje się pod „nocat” (23% pulsarów, reguła A: P₃ grupy ±30% P₃ z params)

Rozkład 119 pulsarów nocat (v2): 64 × tylko grupy z drugiego przejścia (cecha z params bez grupy, jest wolna
modulacja); 20 × grupa przy ≈ ½ P₃ (kandydat na harmoniczną); 19 × przy ≈ 2× P₃; 16 × inne stosunki.
Przykłady (`~/claude/work/scripts/p3track_example_one.jl`, wykresy `~/claude/work/figures/p3track/`):
- **J1807+0756** (P₃ params/Song 19.0): dominująca cecha w danych przy P₃ ≈ 6 (f ≈ 0.17), grupa 567 P z werdyktem
  drift; w średnim widmie brak piku przy f = 1/19. Nie harmoniczna — P₃ z katalogu nie zgadza się z cechą w danych.
- **J1915+0738** (P₃ 37): P₃ wędruje ciągle 17 → 27 → 35 w obserwacji; grupy 16.8 i 24.6 to odcinki jednej wędrówki,
  nie harmoniczna.
- **J1159-6409** (P₃ 13.5): grupy 5.9 i 9.3.
Wniosek: reguła A (zgodność z P₃ katalogowym) myli „inną cechę niż w katalogu” i „wędrujące P₃” z „brakiem grupy”.
Harmoniczna bez fundamentalnej to tylko część (≤ 20 z 119). Poprawka `plot_folds` dla nshuffle = 0.

### 2026-10-02 (cd.) — P–Ṗ: ułamek obserwacji z dryfem (P3Track v2)

Zapisane do rozważenia: werdykt pulsara z pierwszego przejścia (`docs/p3track_method.md` §8.11).
Nowa funkcja `Plot.ppdot_p3track(outdir; results, verdicts=("drift","partial"))`: kolor = f_drift = Σ impulsów w grupach
drift/partial ÷ N (część obserwacji z ciągłym P₃, foldem z kompensacją i istotnym gradientem fazy), puste szare = brak grupy
dryfu. Skrypt `~/claude/work/scripts/p3track_ppdot.jl`, wykres `~/claude/work/figures/ppdot_p3track_fdrift.png`
(+ QNAP `p3track_batch/`). Wyniki v2 (bez werdyktu krótkich grup z folda — v3 niepoliczone).
- drift (Song): 162/412 z grupą dryfu, mediana f = 0.42; P3-only: 5/109, mediana 0.18.
- Wg Ė (drift): 1e29–31: 59% z dryfem (med. f 0.45); 1e31–32: 47% (0.40); 1e32–33: 32% (0.34); 1e33–34: 15% (0.64);
  > 1e34: 22% (0.56, n = 18). Spearman log Ė vs f: −0.29 (wszyscy), −0.05 (tylko z dryfem) — maleje udział pulsarów
  z wykrytym dryfem, nie ułamek czasu dryfu u tych, które go mają.

### 2026-10-02 (cd.) — batch v3: miary stabilności, werdykt krótkich grup z folda, bi-drift

Batch v3 (8 procesów, 32–51 min, 12 × brak danych): `~/output/claude/p3track_batch/p3track_v3.csv`, `p3track_v3_pulsars.csv`,
`p3track_v3_drift_metrics.csv` (dominująca grupa drift/partial na pulsar). Analiza `~/claude/work/scripts/p3track_analysis/v3_metrics.py`.
- Werdykt pulsara (reguła A): drift 130 (v2: 127), partial 27 (26), am 20 (16), inconcl. 103 (111) wśród 412 dryferów;
  P3-only: drift 1, partial 4, am 28 (24), inconcl. 30 (34).
- Krótkie grupy (werdykt z folda): 454 inconclusive, 40 am, 8 drift, 2 partial; zmiany względem blokowego: 40 × inconcl.→am,
  8 × →drift, 2 × →partial. (40 nowych am — do przejrzenia: w kalibracji fold rzadko potwierdzał AM.)
- Bi-drift: J1418-3921, J1537-4912, J1239+2453, J1921+1948, J1843-0211.
- Miary (mediana [kwartyle]): drift — p3_wander 0.089 [0.05–0.12], phase_wander 6.3 [2.9–11.7] cykli/1000 P, koherencja
  0.76 [0.67–0.84]; am — 0.094, 2.9, 0.85; inconclusive — 0.023, 1.8, 0.79.
- Pulsary z dryfem (n = 173), Spearman z log Ė: p3_wander +0.20, phase_wander +0.21, koherencja −0.11; z P₃: +0.39, **−0.55**,
  +0.09 — phase_wander w cyklach/1000 P silnie zależy od P₃ (więcej cykli przy krótkim P₃); do rozważenia normalizacja
  na cykl P₃ (np. cykle fazy na 100 cykli P₃).
- `Plot.ppdot_p3track(...; quantity=:fdrift|:p3_wander|:phase_wander|:coherence)`, wykresy
  `~/claude/work/figures/ppdot_p3track_{fdrift,p3_wander,phase_wander,coherence}.png` (+ QNAP).

### 2026-10-02 (cd.) — zgłoszenie z sesji claude-ac: P3Track gubi dryf przy P₃ ≈ 2 (Nyquist) — zweryfikowane

Druga sesja (tylko odczyt; skrypty `~/claude/work/scripts/play/`, logi `~/claude/work/logs/play/`) wskazała przyczynę:
`feature_peak` szuka maksimum tylko w `srange[2:end-1]`, więc bin f = 0.5 nigdy nim nie jest; cecha przy f₃ i jej lustro
1 − f₃ zlewają się w listku Hanna (±2/L) z maksimum dokładnie w 0.5 → brak maksimum wewnętrznego → `edge` → okno odrzucone.
Sprawdzone tu na v3: P₃ z params 1.8–2.13 (17 pulsarów): 0 × drift (nogroup 7, nocat 8, partial 1, inconclusive 1); 2.13–2.25:
drift 4/16; 2.25–2.5: drift 8/19; najmniejsze P₃ grupy w v3 = 2.15 (J1534-4428). Odcięcie ostre → część „P3-only / brak
dryfu” przy P₃ ≈ 2 to artefakt. (Moja wcześniejsza poprawka — usunięcie górnej osłony — tego nie naprawiła: maksimum nadal
musiało być wewnętrzne.) Testy claude-ac: (A) dopuszczenie maksimum w f = 0.5 z progiem z tasowania; (B) fold przy P₃ = 2
w blokach 32 P (A(φ) = Σ(−1)ⁿxₙ/B vs tasowanie) — czulszy; B + fold w obu aliasach 0.5 ∓ δ: J0846-3533 → drift
(P₃ 2.025, Δψ 0.61, z 7.9; w v2 nogroup). Kierunek dryfu przy Nyquiście nieokreślony (alias). Decyzja o wdrożeniu: użytkownik.

### 2026-10-02 (cd.) — ścieżka Nyquista wdrożona (`nyquist_pass`, `plot_nyquist`)

Decyzja użytkownika: wdrożyć wg rekomendacji claude-ac. Dla P₃ z params ≤ 2.2: test B (bloki 32 P przy P₃ = 2 vs tasowanie)
→ odcinki → f₃ z periodogramu → pełny werdykt z folda tylko gdy M·2δ ≥ 2, inaczej `nyquist`; kierunek dryfu nieokreślony.
Walidacja (`p3track_nyq_check.jl`): J0846-3533 drift (P₃ 2.025, |Δψ| 0.61, z 7.9), J0943+2253 nyquist, kontrola J0924-5302
0 bloków, syntetyki P₃ 2.05: AM 12/12 am, dryf 6/6 drift. Batch: wiersze pass = 3; test na 2 pulsarach OK (usunięty).
Dokument metody §4b. Do zrobienia: batch v4.

### 2026-10-02 (cd.) — porządki na dysku, batch tylko PNG

Usunięte wykresy batchy v1, v2, pilotaży i porównania zapów (QNAP `p3track_batch/figures_*` i `~/claude/work/figures/p3track_batch_*`,
po ~0.7 GB w każdym miejscu); zostały wykresy v3 i wszystkie CSV (QNAP `p3track_batch/` = 343 MB). Wolne: /home 499 GB,
QNAP 1.9 TB (97% zajęte). Funkcje wykresów P3Track mają `pdf=true|false`; batch zapisuje tylko PNG (~½ miejsca),
PDF dowolnego pulsara z `p3track_summary_one.jl`. Skrypt `p3track_batch_run.sh` ustawiony na tag v4.

### 2026-10-02 (cd.) — batch v4: ścieżka Nyquista, sLRFS w wykresach, tylko PNG

Batch v4 (8 procesów, 34–54 min): `~/output/claude/p3track_batch/p3track_v4.csv` (915 grup), `p3track_v4_pulsars.csv`,
wykresy PNG (zbiorczy, sLRFS pierwszego i drugiego przejścia, foldy, Nyquist) w `figures_v4/` (585 MB) i
`~/claude/work/figures/p3track_batch_v4/` (584 MB). Dwa błędy (J1524-5706, J1843-0211: wyzerowane impulsy → 0/0 w normalizacji
panelu widm `plot_track`) poprawione i przeliczone osobno (`p3track_v4_fix_J1524_J1843.csv`, wiersze podmienione w v4.csv).
- Werdykt pulsara (reguła A), drift (412): drift 134 (v3 130), partial 27, am 21, inconcl. 104, nyquist 5, nocat 89, brak 32;
  P3-only (109): drift 1, partial 4, am 28, inconcl. 30, nyquist 6, nocat 22, brak 18.
- Ścieżka Nyquista: wynik u 23 pulsarów, aliasy rozdzielone u 7 → drift 4 (J1502-6128 |Δψ| 1.23 z 12.2, J0846-3533 0.61 / 7.9,
  J1425-5723 1.17 / 12.9, J1848+0604 1.75 / 9.4 — trzy z nich w v3 bez grup), partial 1 (J1517-4356), am 1 (J0624-0424),
  inconclusive 1; `nyquist` 16.
- **Słabość**: 5 werdyktów `nyquist` opiera się na jednym bloku (32 P), 3 na dwóch — przy progu 99% i ~63 blokach na pulsar
  oczekiwane ~0.6 fałszywie istotnego bloku; potrzebny warunek łącznej istotności (np. ≥ 2–3 bloków albo test liczby bloków).
- P–Ṗ odświeżone z v4 (`ppdot_p3track` domyślnie v4; grupa do miar stabilności = największa z miarą, wiersze Nyquista bez miar).

### 2026-10-02 (cd.) — warunek łącznej istotności ścieżki Nyquista; v4b

`nyquist_pass`: liczba istotnych bloków niezachodzących vs Binomial(n, 0.01), p < 10⁻³ (pola `n_ind`, `k_ind`, `p_global`;
kolumny CSV `nyq_blocks`, `nyq_p`). Przeliczone tylko pulsary z P₃ ≤ 2.2 (29, tag v4nyq, ~3 min w 4 procesach) i wstawione:
`~/output/claude/p3track_batch/p3track_v4b.csv`, `p3track_v4b_pulsars.csv`; wykresy tych 29 podmienione w `figures_v4/`
(stare usunięte, w tym nieaktualne `_nyquist.png`).
- Ścieżka zgłasza wynik u 11 pulsarów (było 23), wszystkie p ≤ 7e-7 (≥ 4 bloki): drift 4 (J0846-3533 30/32 bloków,
  J1425-5723 6/32, J1502-6128 4/7, J1848+0604 8/29), partial 1 (J1517-4356), am 1 (J0624-0424), nyquist 5 (J1539-4828,
  J1716-4111, J0855-3331, J0943+2253, J1826-1131). Odpadły werdykty z 1–2 bloków.
- Werdykt pulsara (reguła A), drift: 134/27/21/103 + nyquist 0, nocat 92, brak 35; P3-only: 1/4/28/30 + nyquist 3, nocat 23,
  brak 20. P–Ṗ odświeżone z v4b.

### 2026-10-02 (cd.) — zależności od Ė w v4b

Skrypt `~/claude/work/scripts/p3track_analysis/v4b_edot.py` (S/N: k_snr z travel v4; Ė = 4π²IṖ/P³, I = 10⁴⁵).
- **Wykrycie dryfu (etykieta drift Song+23)**, odsetek z grupą drift/partial: Ė 1e29–31: 61% [56–67]; 1e31–32: 51% [46–55];
  1e32–33: 33% [28–37]; 1e33–34: 17% [13–22]; > 1e34: 22% [14–33] (n = 18). Spearman −0.33; cząstkowa (P₃ kat., k_snr) −0.33;
  logit det ~ log Ė + log k_snr + log P₃: log Ė −0.65 (z −5.9), log k_snr +0.24 (z +2.9), log P₃ −0.40 (z −1.4).
  k_snr i P₃ prawie nie korelują z Ė (+0.11, +0.17) → spadek nie wynika z S/N ani P₃.
- Ułamek czasu z dryfem u tych z dryfem: 0.45 / 0.40 / 0.30 / 0.54 / 0.56 — bez trendu.
- **Stabilność** (dominująca grupa, n ≈ 170), cząstkowa z log Ė przy kontroli P₃ grupy i k_snr [68% bootstrap]:
  p3_wander +0.25 [+0.18, +0.32]; phase_wander +0.21 [+0.13, +0.28]; wędrówka fazy na 100 cykli P₃ +0.28 [+0.19, +0.36];
  koherencja −0.21 [−0.28, −0.13]. Przy wyższym Ė dryf — tam, gdzie jest — mniej stabilny i mniej koherentny (słabo).
- P3-only: drift/partial 0–10% w każdym przedziale (5 pulsarów); werdykt am rośnie z Ė: 1/9, 4/20, 10/39, 9/31, 4/10.
Zastrzeżenie: etykieta Song+23 sama zależy od Ė (selekcja).

### 2026-10-02 (cd.) — aktualizacja `docs/p3track_method.md`

Nagłówek i §0 na stan batcha v4b (tabela werdyktów pełnej próbki, zależność od Ė, Nyquist, bi-drift, przegląd Ė > 10³⁴);
nowe §5.1 (miary stabilności) i §7b (pełna próbka: przebieg, pliki, werdykty, Nyquist, Ė, P–Ṗ, dryfery przy wysokim Ė);
§8 przepisane (12 spraw, w tym reguła werdyktu pulsara, wąskie składowe, werdykty tuż nad progiem, 40 am z folda);
§10 uzupełnione (batch, Nyquist, pola, P–Ṗ). Nowe wykresy w `docs/figures/p3track_{ppdot_fdrift,ppdot_p3_wander,
J0846-3533_nyquist,J1453-6413_fold,J1537-4912_fold}.png`.

### 2026-10-02 (cd.) — wniosek: wobec Song+23 metoda zmienia bardzo niewiele

Przegląd P3-only z drift/partial: J1810-5338 (jedyny drift) wątpliwy — 73 P w 3 odcinkach o różnym P₃, z 5.4 tuż nad progiem;
J1825+0004 przekonujący (dryf na zboczu, mod 1–696); J1057-5226 i J1543+0929 umiarkowane; J1016-5345 wątpliwy. Żaden P3-only
nie przechodzi pewnie do dryfu; 26% P3-only potwierdzone jako am; wśród dryferów 33% drift, 7% partial, ~55% bez
rozstrzygnięcia. Dopisane do §0 i §7b `docs/p3track_method.md` (wartość metody: gdzie/kiedy/jak stabilny dryf, Ė, P₃ ≈ 2).

### 2026-10-02 (cd.) — usunięte wykresy batcha v3

Usunięte `~/output/claude/p3track_batch/figures_v3` i `~/claude/work/figures/p3track_batch_v3` (po 340 MB; zastąpione przez v4).
Zostaje: QNAP `p3track_batch/` 587 MB (figures_v4 + wszystkie CSV + P–Ṗ), `~/claude/work/figures/` 667 MB (p3track_batch_v4,
p3track/ testy, P–Ṗ). Razem ~1.25 GB.

### 2026-10-02 (cd.) — zamknięcie sesji P3Track

Skrypty analizy CSV przeniesione z katalogu tymczasowego sesji do `~/claude/work/scripts/p3track_analysis/` (odnośniki
w dokumentach poprawione). Stan prac i wskazówki dla nowej sesji (szukanie lepszej metody): `docs/p3track_method.md` §0, §8, §9;
`docs/travel_test_method.md`; ten dziennik.

### 2026-10-02 (cd.) — nowa sesja: czy ciągłe P₃(n) z `p3fold_coherent` odróżnia dryf od P3-only?

Pomysł użytkownika: `coherent_fold` daje ciągłe P₃(n), może to ułatwi wykrycie dryferów. Test na zestawie kontrolnym P3Track
(5 dryferów + 5 P3-only), te same parametry co `SpaTs.p3fold_coherent` (low-pass 1/300, rząd 6, jackknife 4 grupy).
Dodatkowo Δψ z jednego złożenia koherentnego na pulsar (C = (I − ⟨I⟩)·e^{−iθ(n)}, `template_phase` + `template_significance`
z P3Track, bootstrap blokami L/2 i impuls po impulsie). Skrypt `~/claude/work/scripts/coh_control.jl`, log
`~/claude/work/logs/coh_control.log`, wykresy `~/claude/work/figures/coh_control/<PSR>_coh.png`.

| etykieta | PSR | SNR | wędr. P₃ std/med | wędr. fazy / 100 cykli P₃ | Δψ (składowe) | z | werdykt | P3Track |
|---|---|---|---|---|---|---|---|---|
| P3-only | J1401-6357 | 27.0 | 0.009 | 0.27 | 0.17 (−0.17 ± 0.18) | 0.4 | inconclusive | am |
| P3-only | J1603-2531 | 23.9 | 0.113 | 3.61 | 0.08 (+0.08 ± 0.05) | 1.0 | am | am, am, inconcl. |
| P3-only | J1001-5939 | 20.5 | 0.010 | 0.26 | — (składowa 5 binów, brak fragmentu) | — | inconclusive | brak grup |
| P3-only | J1146-6030 | 7.2 | 0.027 | 1.44 | 0.17 (3 skł., σ 0.14–0.20) | −1.5 | inconclusive | inconclusive |
| P3-only | J2307+2225 | 3.1 | 0.010 | 0.41 | 0.17 (−0.17 ± 0.12) | 1.1 | inconclusive | brak grup |
| drift | J0034-0721 | 20.4 | 0.019 | 0.57 | 1.67 (−1.67 ± 0.12) | 14.1 | drift | drift |
| drift | J0151-0635 | 5.0 | 0.038 | 1.81 | 1.12 (−0.87, −0.24) | 31.1 | drift | drift |
| drift | J0820-1350 | 89.5 | 0.007 | 0.33 | 1.55 (−1.55 ± 0.03) | 37.0 | drift | drift |
| drift | J1825+0004 | 3.1 | 0.063 | 0.63 | 0.37 (−0.34 ± 0.15, −0.04) | 1.4 | inconclusive | partial |
| drift | J1750-3503 | 2.1 | 0.061 | 2.44 | 1.57 (+1.57 ± 0.21) | 7.5 | drift | drift |

- **Ciągłe P₃(n) nie rozdziela klas**: wędrówka P₃ 0.007–0.063 (drift) vs 0.009–0.113 (P3-only), wędrówka fazy 0.33–2.44 vs
  0.26–3.61 — zakresy się pokrywają (jak p3_wander w P3Track v4b: 0.089 vs 0.094). Powód z konstrukcji: filtr dopasowany
  Σ_φ x·conj(L) zjada ψ(φ), więc θ(n) i P₃(n) są takie same dla dryfu i AM; dryf siedzi tylko w arg L(φ) / w złożeniu.
- P₃(n) przy dużej wędrówce jest zawodne: J1603 (P3Track: P₃ 13–52) daje 45–53 z ostrymi pikami w miejscach poślizgu fazy
  (spadek amplitudy po low-passie → skok rozwiniętej fazy). Low-pass 1/300 nie nadąża za szybkimi zmianami.
- **Δψ z jednego złożenia na pulsar** daje te same werdykty co P3Track dla 4 dryferów i J1603 (am), 0 fałszywych dryfów;
  gorzej tam, gdzie dryf jest tylko w części obserwacji (J1825: rozmyty przez impulsy po 715) i przy P₃ ≈ 2 (J1401: am → inconcl.).
- Wniosek: P₃(n) samo nie jest dyskryminatorem; złożenie koherentne + faza szablonu ≈ uproszczony P3Track bez grup.

### 2026-10-02 (cd.) — Δψ(t): gradient fazy w poprzek składowej rozdzielczy w czasie (test)

Pomysł: zamiast jednego szablonu na obserwację — lokalne szablony w niezachodzących oknach W = max(16, 4·P₃),
T_w(φ) = ⟨(I − ⟨I⟩)·e^{−iθ(n)}⟩, θ(n) z `coherent_fold` (jak `p3fold_coherent`); gradient Δψ_w w każdej składowej
(Σ T*_j T_{j+1}, maska z mocy niekoherentnej minus podłoga szumu z off-pulse'u). Trwałość bez znaku: okno dzielone na połowy,
Z_w = Δψ_a·Δψ_b/(σ_a σ_b) (σ z bootstrapu impulsów w połowie), z = ΣZ_w/√(n_okien·n_skł) — dryf daje Z > 0 także przy zmianie
kierunku między oknami. Ze znakiem: średnia ważona 1/σ², χ² po składowych. Skrypt `~/claude/work/scripts/dpsi_time.jl synth|real`,
logi `~/claude/work/logs/dpsi_time_{synth,real}.log`, wykresy `~/claude/work/figures/dpsi_time/` (+ `synth/`).
(Pierwsza wersja z nullem z losowania znaków okien nasycała się: z ≤ 3.9 przy ~31 oknach.)

Syntetyki (1000 P, P₃ = 8): AM wspólna faza / przeciwfaza / dudnienie 8.0 + 8.5 → trwałość z ≤ 1.1 (100 przypadków; statystyka
konserwatywna, rozrzut pod AM ≈ 0.2 zamiast 1); dryf z = 95 (szum 0.6) / 3.4 (1.5, ze znakiem 9.4); reverser (epizody ~100 P)
75 / 2.8, ze znakiem 10.8 / 1.9; epizodyczny (dryf w 1/3 obserwacji) 30 / 1.3; wolny dryf 4.9 / 0.2.

Dane (zestaw kontrolny):

| etykieta | PSR | W | okna | Δψ glob. | trwałość z | znak z |
|---|---|---|---|---|---|---|
| P3-only | J1401-6357 | 16 | 64 | −0.04 | −2.4 | 1.3 |
| P3-only | J1603-2531 | 194 | 10 | +0.01 | −0.9 | −0.8 |
| P3-only | J1001-5939 | 16 | 32 | — (składowa < 5 binów) | — | — |
| P3-only | J1146-6030 | 44 | 24 | +0.12 | −0.1 | 0.9 |
| P3-only | J2307+2225 | 16 | 64 | −0.04 | −1.0 | 0.8 |
| drift | J0034-0721 | 26 | 40 | −0.73 | 10.0 | 10.3 |
| drift | J0151-0635 | 58 | 17 | −0.76, −0.18 | 53.1 | 27.2 |
| drift | J0820-1350 | 18 | 58 | −0.70 | 64.9 | 23.2 |
| drift | J1825+0004 | 56 | 18 | −0.19, +0.06 | 0.0 | 1.4 |
| | J1825 1–696 / 716–1040 | 56 | 12 / 5 | −0.12 / +0.27 | 0.6 / −0.1 | 1.9 / 1.0 |
| drift | J1750-3503 | 196 | 5 | −0.09, +0.69 | 0.7 | 3.0 |

- Silne dryfery wykryte, P3-only bez fałszywych detekcji — ale to samo daje P3Track i złożenie koherentne.
- **Przypadki, dla których metoda była pomyślana, nie wychodzą**: J1825 (dryf epizodyczny, słaby, tylko na zboczu składowej —
  gradient całej składowej go rozcieńcza; na wykresie okna 100–300 mają Δψ ≈ −1, ale pojedyncze okna za szumne);
  J1750 (reverser): P₃ ≈ 49 → W = 196, 5 okien, a epizody 28 ± 4 P ≪ P₃ — fazy modulacji nie da się zmierzyć w czasie krótszym
  niż cykl P₃, więc żadna metoda fazowa tego nie rozdzieli.
- Δψ glob. zaniżone wobec P3Track (J0034: −0.73 vs −1.61) — okna z nullami i słabą modulacją wchodzą z wagą 1/σ², ale ich Δψ
  jest bliskie przypadkowemu. J1401 (P₃ 2.2): z = −2.4 — przy P₃ ≈ 2 człon 2f po demodulacji aliasuje blisko 0 i nie uśrednia się
  w połowie okna (8 P).
- Wniosek: Δψ(t) daje czytelny obraz „gdzie w czasie” dla jasnych pulsarów, ale nie zwiększa czułości tam, gdzie P3Track i
  Song+23 się rozjeżdżają (słabe, epizodyczne, długie P₃).

### 2026-10-02 (cd.) — śledzenie podpulsów impuls po impulsie (`subtrack`), automatyczna wersja Szary+2022 §3.2

Metoda z pracy o J1750-3503 (ApJ 934, 23): splot impulsu z Gaussem o szerokości podpulsu, maksima z S/N ≥ 5, łączenie w pasma
(tam wzrokowo), tempo dryfu z dopasowania ścieżek. Tu łączenie automatyczne. Skrypt `~/claude/work/scripts/subtrack.jl`
(`[katalog plik]...`, domyślnie zestaw kontrolny), siatka parametrów `grid_subtrack.jl`, syntetyki `subtrack_synth.jl`;
logi `~/claude/work/logs/subtrack_*.log`, wykresy `~/claude/work/figures/subtrack/<PSR>_subtrack[_zoom].png` (+ `synth/`).
- Detekcja: σ podpulsu z autokorelacji fluktuacji (HWHM/(√2·√(2 ln 2))), S/N względem off-pulse'u po tym samym splocie,
  maksima bliższe niż FWHM tłumione, pozycja z paraboli.
- Łączenie: zachłanne jeden-do-jednego do pozycji przewidywanej (ostatnia + nachylenie z ostatnich ≤ 5 punktów), skok ≤ FWHM
  podpulsu, luka ≤ 2 impulsy; dopasowanie Theil–Sen. Null: tasowanie kolejności impulsów (te same podpulsy; AM daje te same
  pionowe ścieżki, dryf się rozpada).
- **Kalibracja na J1750-3503** (siatka S/N 4/5 × skok P₂/2, 1.5, 1.0, 0.6 FWHM × luka 1/2 × LS/TS): D₊ = +0.395, D₋ = −0.338 °/P
  (praca: 0.388, −0.314), D > 0 przez 74% czasu ze ścieżkami (praca ~78%), P₂ = 18.1° (18.6°). Zawyżone D (0.55–0.60) z pierwszej
  wersji brało się z selekcji ścieżek po |D|/σ ≥ 3 (krótkie strome ścieżki); D raportowane jako średnia ważona długością ścieżek
  ≥ 8 punktów bez selekcji.
- **Statystyka dryfu**: f_ls = udział podpulsów w ścieżkach ≥ 8 punktów z |D|/σ ≥ 3; z(f_ls) względem 20 tasowań.
  (z(f_long) — same długie ścieżki — myli pionowe ścieżki AM, J1825: 7.2; z(Q) = ⟨|D|/σ⟩ słabe dla J1750: 1.7.)

| etykieta | PSR | P₂ [°] | f_ls (tas.) | z(f_ls) | D₊ / D₋ [°/P] (ścieżki) |
|---|---|---|---|---|---|
| P3-only | J1401-6357 | 4.7 | 0.24 (0.24) | −0.1 | +0.13 (26) / −0.11 (35) |
| P3-only | J1603-2531 | 11.4 | 0.22 (0.18) | 1.2 | +0.14 (65) / −0.12 (32) |
| P3-only | J1001-5939 | 4.2 | 0.00 (0.00) | −0.2 | ≈ 0 |
| P3-only | J1146-6030 | 5.5 | 0.30 (0.28) | 1.0 | +0.11 (60) / −0.14 (55) |
| P3-only | J2307+2225 | 11.9 | 0.00 (0.02) | −0.3 | ≈ 0 |
| drift | J0034-0721 | 18.1 | 0.67 (0.13) | 14.8 | −2.36 (66) |
| drift | J0151-0635 | 14.1 | 0.56 (0.19) | 17.0 | −0.40 (109) |
| drift | J0820-1350 | 3.6 | 0.85 (0.35) | 8.9 | −0.73 (184) |
| drift | J1825+0004 | 9.0 | 0.03 (0.03) | 0.2 | ≈ 0 (wykrywa jasny szczyt AM, dryf słaby na zboczu) |
| drift | J1750-3503 | 18.1 | 0.43 (0.09) | 13.5 | +0.40 (27) / −0.34 (13) — reverser |

Syntetyki (1000 P, podpulsy σ = 3 biny, P₂ = 24, 20 ziaren AM / 5 dryf; szum 1.0 / 2.0 ≈ S/N podpulsu 9 / 4.6): AM stałe pozycje,
AM przeciwfaza, losowe podpulsy → max z(f_ls) = 2.2 w 120 przypadkach, 0 × z ≥ 5. Dryf z = 30 / 1.5 (D = 2.96 przy 3.0), reverser
33 / 1.6 (D ±2.9), wolny dryf (P₃ = 40) 28 / 33 (D 0.59–0.60 przy 0.6), alias P₃ = 2.1: z = 22, ale D błędne (0.27 / −0.57 zamiast
11.4 lub −12.6) — przy P₃ ≈ 2 łączenie łapie alias. Przy S/N podpulsu < 5 metoda nic nie widzi (to ograniczenie z założenia).

P3-only z oznakami dryfu w travel/P3Track (`subtrack_p3only_cand.log`): J1810-5338 z 0.9, J1543+0929 0.9, J1016-5345 −1.5,
J0837+0610 −0.1, J1633-4453 0.3, J1701-3130 1.9; J1057-5226 z = 6.8 i J1048-5832 4.3, ale efekt znikomy (f_ls 0.07 vs 0.05,
0.29 vs 0.25) przy 27 401 i 7 270 impulsach (mały rozrzut tasowań), D > 0 w 49% / 54% czasu — bez spójnego kierunku.
Potrzebny warunek na wielkość efektu (np. f_ls − f_ls,tas ≥ 0.1; dryfery kontrolne 0.34–0.54, P3-only ≤ 0.04).
Wniosek wstępny: żaden z tych P3-only nie ma pasm podpulsów — zgodnie z Song+23.

### 2026-10-02 (cd.) — pilot: modulacja P₃ w polaryzacji (I, Q′, U′, V) na zestawie kontrolnym

Pomysł: metody fazowe (travel, P3Track, Δψ(t)) mierzą ψ(φ), czyli to samo co 2DFS Song+23. Polaryzacja to nowa informacja
(Song+23 jej nie używają). Archiwa `_16/pulsar.spCf16` mają pełne Stokesy, `rmc=1`, `polc=1`, `scale=FluxDensity`.
Skrypt `~/claude/work/scripts/pol_pilot.jl`, log `~/claude/work/logs/pol_pilot.log`, wykresy `~/claude/work/figures/pol/<PSR>_pol.png`,
cache Stokesów (pdv -t -F) `~/claude/work/pol/cache/*.jld2` (~17 MB/pulsar).

Konstrukcja: wielkości liniowe I, V, Q′, U′ (Q, U obrócone do χ_ref(φ) z ⟨L e^{4iχ}⟩ — niewrażliwe na OPM; Q′ > 0 mod główny,
Q′ < 0 ortogonalny). Reszty X_res = X − (X̄/Ī)(φ)·I oraz X_res2 = reszta z regresji X ~ c + aI + bI² w każdym binie (kontrola
zależności ułamka polaryzacji od jasności). Cecha P₃ w resztach = modulowany jest **ułamek** polaryzacji w danej długości, nie
tylko jasność. Statystyka **nie używa zależności fazy od długości** (każdy bin osobno) — ortogonalna do Δψ / 2DFS.
z: nadwyżka mocy LRFS w paśmie f₃ (FWHM piku I) względem 100 tasowań kolejności impulsów; też w paśmie 2f₃.

| etykieta | PSR | P₃ | z_I | z Q′res2 | z U′res2 | z Vres2 | E_Q′res2/E_I | 2f₃: z_I | 2f₃: max z res2 |
|---|---|---|---|---|---|---|---|---|---|
| P3-only | J1401-6357 | 2.21 | 4.6 | 2.3 | −1.4 | 0.4 | 0.054 | −0.5 | 1.1 |
| P3-only | J1603-2531 | 48.6 | 18.2 | 1.2 | 3.2 | −1.1 | 0.013 | 9.2 | 1.3 |
| P3-only | J1001-5939 | 2.09 | 9.7 | **8.2** | 1.3 | 1.5 | 0.087 | 0.9 | 1.7 |
| P3-only | J1146-6030 | 10.9 | 3.8 | 2.2 | 2.4 | 4.6 | 0.151 | −0.1 | 0.5 |
| P3-only | J2307+2225 | 3.48 | 4.7 | 1.4 | 0.1 | −0.3 | 0.122 | −0.5 | 1.1 |
| drift | J0034-0721 | 6.57 | 9.7 | **22.3** | 5.0 | 2.9 | 0.179 | −1.9 | 4.6 (U′) |
| drift | J0151-0635 | 14.3 | 72.2 | **65.5** | 5.0 | 8.3 | 0.104 | 5.2 | 10.3 (Q′) |
| drift | J0820-1350 | 4.77 | 89.0 | **105.5** | 36.6 | 75.1 | 0.162 | 1.3 | **23.7 (V)** |
| drift | J1825+0004 | 14.2 | 7.1 | 2.1 | 3.1 | 4.2 | 0.069 | 0.6 | 1.1 |
| drift | J1750-3503 | 49.0 | 32.7 | **12.6** | 3.6 | 2.6 | 0.265 | 21.7 | 6.1 (Q′) |

- **Najczystszy kontrast: J1603 (pewne AM) vs dryfery.** J1603: silna modulacja I (z 18), a ułamek polaryzacji prawie niemodulowany
  (Q′res2 z 1.2, E/E_I 0.013). Dryfery: Q′res2 z 13–105, E/E_I 0.10–0.27. Interpretacja: przy AM podpuls jaśnieje w miejscu i
  ułamek polaryzacji w danej długości się nie zmienia; przy dryfie podpuls ze swoją strukturą polaryzacji (OPM na brzegach,
  zmiana znaku V) przechodzi przez daną długość.
- **2f₃ w polaryzacji bez 2f₃ w I**: J0820 (V z 24, U′ 23 przy I 1.3), J0034 (U′ 4.6, V 3.6 przy I −1.9) — struktura polaryzacji
  w obrębie cyklu podpulsu (np. OPM dwa razy na cykl). Wśród P3-only brak (max 1.7).
- J0034: modulacja Q′ skupiona w długości mieszania OPM (spadek L przy bin ~490) — ułamek modów zmienia się z P₃.
- **J1001-5939 (P3-only, P₃ 2.09)**: Q′res2 z 8.2 — kandydat (alias dryfu przy Nyquiście albo okresowe OPM); P3Track: brak grup.
- **Zastrzeżenia**: z skaluje się z S/N (dryfery w zestawie jaśniejsze); E/E_I zależy od średniego L/I; J1825 (słaby, epizodyczny)
  i słabe P3-only (J1146, J2307) poniżej czułości. Statystyka nie odróżni dryfu od okresowego przełączania OPM bez dryfu —
  to raczej pytanie fizyczne niż wada. Wynik na 10 obiektach, progi niewyznaczone.

### 2026-10-02 (cd.) — subtrack: długie obserwacje przycięte do 1000 impulsów

Decyzja użytkownika: przy długich danych brać tylko część. `subtrack.jl`: `MAXPULSES = 1000`, analizowany ciągły fragment
1–1000 (z porównywalne między pulsarami; zestaw kontrolny ma ~1000–2100 P). Log `~/claude/work/logs/subtrack_maxp.log`.
- J1057-5226 (27 401 P): z(f_ls) 6.8 → **−0.1**; J1048-5832 (7270 P): 4.3 → **1.3**; J1603-2531 0.7, J1810-5338 0.8, J1701-3130 1.0.
- Zestaw kontrolny bez zmian jakościowych: dryfery z = 14.3 (J0034), 13.0 (J0151), 11.1 (J0820), 11.7 (J1750); J1825 0.5;
  P3-only −0.3…0.7.
Batch na pełnej próbce jeszcze nie puszczony (decyzja użytkownika).

### 2026-10-02 (cd.) — subtrack na losowej próbce 10 drift + 10 P3-only: dryferów nie wykrywa

Losowanie (Python `random.Random(20261002)`, bez 18 pulsarów już oglądanych; pula 408 drift / 95 P3-only z danymi), lista
`~/claude/work/subtrack_random20.txt`, log `~/claude/work/logs/subtrack_random20.log`, test krótkich ścieżek `subtrack_short.jl` /
`subtrack_short.log`. Każdy pulsar: pierwsze 1000 P. `subtrack.jl` przyjmuje teraz listę `.txt` (etykieta katalog plik) i łapie błędy.

| etykieta | PSR | σ_sub [°] | imp. z podp. | P₂ [°] | z(f_ls) | z(f_sig, ścieżki ≥ 4) | P3Track v4b |
|---|---|---|---|---|---|---|---|
| drift | J1119-7936 | 0.60 | 875 | 4.5 | 0.1 | 0.0 | inconclusive |
| drift | J1910+0714 | 0.54 | 789 | 2.4 | −2.9 | 1.4 | **drift** |
| drift | J1055-6905 | 1.88 | 447 | 9.7 | 0.0 | −1.5 | nocat (am) |
| drift | J1823-0154 | 0.76 | 879 | 2.9 | −0.7 | −0.2 | brak grup |
| drift | J1614+0737 | 0.93 | 815 | 3.7 | −0.5 | −0.5 | inconclusive |
| drift | J1645-0317 | 0.85 | 1000 | 7.2 | 0.0 | 0.1 | **drift** |
| drift | J1916+1030 | — | — | — | ACF płaska (szum) | — | brak grup |
| drift | J1549+2113 | 2.36 | 497 | — | −0.2 | 0.2 | inconclusive |
| drift | J0932-3217 | 0.47 | 155/156 | — | 0.0 | 0.4 | **drift** |
| drift | J1903+2225 | 1.57 | 376 | 6.9 | 1.3 | 2.5 | nocat |
| P3-only (10) | J1801-0357 … J1717-3425 | | | | −0.9…2.6 | −0.4…2.5 | |

- **0/10 dryferów wykrytych** (P3Track: 3/10 drift); P3-only bez fałszywych (max 2.6, J1801-0357 f_ls 0.14 vs 0.06).
  Krótkie ścieżki (≥ 4 punkty) nie pomagają.
- Przyczyna (wykresy J1910+0714, J1645-0317): typowy dryfer próbki ma wąską składową (σ_sub 0.5–0.9°, P₂ 2–7°), 1–2 podpulsy
  na impuls i krótkie P₃ — pasmo żyje ~W/|D| ≈ 3–5 impulsów (J1910: nachylone kreski widoczne gołym okiem, ale ścieżki ≤ 4 P),
  a przesunięcie podpulsu w jego czasie życia jest porównywalne z jego szerokością. Do tego słabe S/N (J1916: ACF fluktuacji
  płaska). Dryf jest tam mierzalny tylko statystycznie (gradient fazy / 2DFS), nie jako ruch pojedynczych podpulsów.
- Wniosek: subtrack działa dla „podręcznikowych” jasnych dryferów z szerokim profilem i wieloma podpulsami (zestaw kontrolny:
  J0034, J0151, J0820, J1750), a nie dla typowego pulsara z listy Song+23. Jako metoda klasyfikacji całej próbki — nie;
  jako narzędzie do szczegółowego opisu (D(t), reverserzy, P₂) jasnych dryferów — tak.

### 2026-10-02 (cd.) — pairshift: przesunięcia podpulsów między kolejnymi impulsami (bez ścieżek)

Użytkownik: w J1910+0714 dryf dobrze widoczny, a subtrack go nie widzi (pasma żyją 3–5 P, ścieżki ≥ 8 P niemożliwe).
Nowa statystyka: detekcja jak w subtrack; w każdej parze (n, n+1) dopasowanie najbliższego podpulsu w obu kierunkach (n → n+1
i n+1 → n, przesunięcia zawsze w kierunku czasu), promień R; s_n = Σ sign(Δφ). Odwrócenie czasu pary zmienia znak s_n
dokładnie, więc null = losowanie znaków par, z = Σs/√Σs² (bez modelu zmienności, ~1000 par, bez nasycenia).
Skrypty `~/claude/work/scripts/pairshift.jl` (kontrola + losowe 20), `pairshift_grid.jl` (wagi), `pairshift_synth.jl`;
logi `~/claude/work/logs/pairshift*.log`, wykresy `~/claude/work/figures/pairshift/<PSR>_pairs.png`.
- Pierwsza wersja (tylko n → n+1, waga Δ/R): fałszywe dryfy w P3-only (J0629+2415 −11.0, J1825 −7.2, J1401 5.3) — dopasowanie
  jednokierunkowe nie jest antysymetryczne przy różnej liczbie podpulsów w impulsach. Po dopasowaniu dwukierunkowym null czysty,
  ale waga liniowa daje J1910 tylko −2.2 (dalekie dopasowania szumu); siatka wag (`pairshift_grid.log`): **znak, R = 1.5·FWHM**.
- Syntetyki (jak subtrack): AM stałe / przeciwfaza / losowe → max |z| 2.4 w 120, 0 × |z| ≥ 3; dryf 25.8 / **10.6** (szum 1 / 2;
  subtrack przy szumie 2: 1.5), reverser 8.0 / 2.2, wolny dryf 8.6 / 1.7, alias P₃ 2.1: 5.5 / 2.0. Wersja blokowa (bez znaku, dla
  reverserów) nasyca się przy √31 ≈ 5.6 i niewiele daje.

| zestaw | etykieta | PSR | mediana Δφ [°/P] | z | P3Track v4b |
|---|---|---|---|---|---|
| kontrola | P3-only (6, w tym J1825) | | | −1.1 … 0.9 | |
| kontrola | drift | J0034-0721 | −1.96 | **−13.6** | drift |
| kontrola | drift | J0151-0635 | −0.32 | **−9.7** | drift |
| kontrola | drift | J0820-1350 | −0.71 | **−27.9** | drift |
| kontrola | drift | J1750-3503 (reverser) | +0.25 | 1.7 | drift |
| losowe | drift | J1910+0714 | −0.41 | **−10.9** | drift |
| losowe | drift | J0932-3217 | +0.22 | **3.3** | drift |
| losowe | drift | J1614+0737 / J1903+2225 / J1645-0317 | | −2.6 / −2.2 / −1.9 | inconcl. / nocat / drift |
| losowe | drift | pozostałe 4 (+ J1916 bez detekcji) | | −1.4 … 0.3 | |
| losowe | P3-only (10) | | | −2.4 … 1.1 | |

- Przy |z| ≥ 3: dryfery kontrolne 3/4, losowe dryfery 2/9 (P3Track: 3/10 — J1910, J1645, J0932); P3-only 0/16.
- J1910+0714 rozwiązany (−10.9, D ≈ −0.41 °/P), J0820 D zgodne z subtrack (−0.71 vs −0.73).
- Wniosek: pairshift wykrywa krótkie pasma, których subtrack nie widzi, i ma czysty null; na losowej próbce mniej więcej tyle co
  P3Track. Ograniczenie zostaje to samo: S/N pojedynczych podpulsów.

### 2026-10-02 (cd.) — batch pairshift v1 na pełnej próbce (533) + subtrack

`~/claude/work/scripts/pairshift_batch.jl` + `pairshift_batch_run.sh` (8 procesów, wznawialny, pierwsze 1000 P), podsumowanie
`pairshift_summary.py`. Wyniki `~/output/claude/pairshift_batch/pairshift_v1.csv` (wiersz na pulsar; części `_partKof8`;
`pairshift_test_part1of1.csv` = test na 3 pulsarach, do usunięcia), wykresy (tylko |z| ≥ 3) `~/claude/work/figures/pairshift_batch_v1/`.
Błędy: 12 × brak danych, 29 × ACF płaska (za słabe na detekcję podpulsów).

| Song+23 | policzone | pairshift \|z\| ≥ 3 | subtrack z(f_ls) ≥ 5 |
|---|---|---|---|
| drift | 388 | **165 (43%)** | 58 |
| P3-only | 103 | **0** | 6 |

- **Null na P3-only**: z od −2.6 do 2.5, mediana −0.1, |z| ≥ 2 u 6/103 (oczekiwane ~5 dla N(0,1)).
- Dryfery krzyżowo z P3Track v4b: pairshift+ i P3Track drift/partial 100, tylko pairshift 65 (P3Track inconclusive 33, nocat 26,
  am 4, nogroup 2), tylko P3Track 56; **suma 221/388 (57%)** wobec 156 (40%) z samego P3Track. Nowe detekcje: mediana |z| 4.1
  (30 z nich w 3–4), najsilniejsze J0533+0402 −16.2 (P3Track inconclusive; krótkie opadające kreski widoczne na stosie),
  J0924-5814 +11.0, J1627-5936 −10.0, J1807+0756 −9.8 (nocat), J1850+0026 −9.1.
- **Kierunek dryfu**: znak z pairshift vs znak Δψ (ważony mocą) dominującej grupy drift/partial P3Track — zgodny 97, przeciwny 2
  (J1741-0840, J1819+1305). Niezależne potwierdzenie konwencji znaku P3Track (sprawa otwarta §8.11 w p3track_method.md).
- P₃ ≤ 2.2: pairshift 2/20 dryferów — przy Nyquiście dopasowanie najbliższego łapie alias.
- Sprzeczne: pairshift+ przy P3Track `am`: J1740+1311, J1946-2913, J1839-1238, J1709-4429 (|z| 3.4–3.7, Δφ ~0.1 °/P) — do obejrzenia.
- **subtrack: fałszywe z(f_ls) w P3-only** (6, np. J1531-5610 z = 64, f_ls 0.06 vs 0.00): długie pionowe ścieżki mają tak małe σ_D,
  że |D|/σ ≥ 3 przy D ≈ 0. Poprawka do zrobienia: wymagać przesunięcia |D|·długość ≥ FWHM podpulsu.
  **J1651-1709** (P3-only, P3Track `am`): wykres `~/claude/work/figures/subtrack/J1651-1709_subtrack_zoom.png` — wiodąca składowa
  (177°) stała, ale w końcowej (183–186°) widać rosnące pasma co ~20–25 P (D ≈ +0.1 °/P); pairshift z = 0.35 (pary w jasnej
  składowej bez przesunięcia rozcieńczają statystykę). Kandydat na dryf w jednej składowej — do oceny wzrokowej.
