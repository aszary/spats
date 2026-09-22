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

# START TUTAJ — stan na 2026-09-17

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
