# Analiza separations — dziennik i przewodnik

Wyznaczanie Δseparation (offset składowych profilu między 1023 a 1523 MHz) dla listy
`input/separations_todo.csv`. Wyniki lądują w `input/separations.csv`.
Powiązany dokument: [`component_position_methods.md`](component_position_methods.md).

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
