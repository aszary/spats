# Separacje Macieja — porównanie z pomiarami w `input/`

**Stan na 2026-09-18.** Konfrontacja niezależnego zestawu pomiarów separacji składowych
(`input/separations_maciej.csv`, 50 pulsarów) z pomiarami już obecnymi w repozytorium.

Repozytorium: `github.com/aszary/spats`, gałąź `claude`.
Dziennik roboczy: `docs/separations_analysis_log.md`.

Skrypty użyte do analizy (poza repo, w `~/claude/work/scripts/`):
`cmp_maciej.py` (tabele porównawcze), `profile_check.py` (profile średnie).
Log: `~/claude/work/logs/cmp_maciej.log`. Wykres: `~/output/claude/cmp_maciej_profile.png`.

---

## 1. Plik wejściowy

`input/separations_maciej.csv` — 50 pulsarów, format jak `input/separations.csv`, ale
**bez kolumn `dsep`, `dsep err`** (12 kolumn zamiast 14). Zawiera tylko pozycje składowych
i separację, nie zawiera offsetu częstotliwościowego.

Nazwę pliku poprawiono przy wczytaniu z `separartions_maciej.csv` (literówka w oryginale).

## 2. Pokrycie

| plik | wierszy | wspólnych z Maciejem | brakujących |
|---|---|---|---|
| `separations.csv` | 91 | 47 | J1627-5936, J1819+1305, J1932+1059 |
| `separations_merged.csv` | 195 | 48 | J1627-5936, J1932+1059 |
| `separations_todo.csv` | 106 | **50** | — |

Wszystkie 50 pulsarów Macieja jest już na liście `separations_todo.csv`, i **wszystkie 50 mają
tam puste pole `zrobione`**. Żaden pulsar nie jest nowy dla próbki.

`ncomp` zgadza się we wszystkich 47 wspólnych przypadkach (wartości 2, 3, 4).

## 3. Zgodność separacji

Względem `separations.csv`: **41/47 zgodnych w granicach 3σ**, z czego **32 identyczne co do
czwartego miejsca po przecinku** (diff = 0.0000). Pozostałe zgodne mieszczą się poniżej 0.6σ.
To znaczy, że w większości przypadków chodzi o ten sam pomiar, a nie o niezależne powtórzenie.

### 3.1 Przypadki rozbieżne (>3σ)

| pulsar | Maciej | `separations.csv` | diff | σ | źródło różnicy |
|---|---|---|---|---|---|
| J1733-3716 | 34.3649 | 39.5506 | −5.19 | 12.1 | lon 2 (−5.18°) |
| J1901+0716 | 7.5265 | 5.6621 | +1.86 | 10.1 | lon 1 (−1.91°) |
| J1757-2421 | 20.6813 | 18.7415 | +1.94 | 7.8 | lon 1 (−1.96°) |
| J1714-1054 | 11.2534 | 11.8583 | −0.60 | 5.9 | lon 2 (−0.61°) |
| J1803-3329 | 4.7241 | 4.1425 | +0.58 | 4.2 | lon 1 (−0.72°) |
| J1808-3249 | 11.1026 | 9.6328 | +1.47 | 3.9 | lon 1 (−1.37°), lon 3 (+0.10°) |

We wszystkich przypadkach rozjeżdża się **jeden skrajny komponent**, nie oba — to inne
dopasowanie zewnętrznej składowej, a nie przesunięcie fazy całego profilu.

Względem `separations_merged.csv` rozbieżnych >3σ jest siedem: powyższa szóstka (dla której
`merged` bierze wartości ze źródła `stare`) plus **J1819+1305** (18.5010 vs 13.3531, 11.9σ),
którego nie ma w `separations.csv`.

---

## 4. Rozstrzygnięcie — rekonstrukcja z zapisanych fitów per-pulse

Cztery z siedmiu spornych pulsarów mają zachowany plik `~/output/claude/<PSR>_16/component_offsets.txt`
z surowymi `mu` dopasowań gaussowskich dla każdego pulsu P3-fold, osobno dla low i high.
Odtworzono na nich procedurę `Plot._offset_summary`: `lon = (mu_high + mu_low)/2`, waga `1/err²`,
błąd zewnętrzny `lon_ext = (1/√Σw)/2 · √max(1, χ²/dof)`.

| pulsar | rekonstrukcja z fitów | `separations.csv` | Maciej | `merged` |
|---|---|---|---|---|
| J1757-2421 | **18.652 ± 0.132** | 18.7415 | 20.6813 | 18.7415 |
| J1803-3329 | **4.075 ± 0.082** | 4.1425 | 4.7241 | 4.1425 |
| J1808-3249 | **9.591 ± 0.112** | 9.6328 | 11.1026 | 9.6328 |
| J1819+1305 | **17.403 ± 0.134** | — | 18.5010 | 13.3531 |

Dla pierwszych trzech pulsarów `separations.csv` odtwarza się z surowych danych; różnice
0.05–0.09° odpowiadają maskowaniu części pulsów (`keep`) przy oryginalnym przebiegu.
Wartości Macieja odbiegają od rekonstrukcji o 2–13σ. To nie jest kwestia wyboru metody —
liczby w `separations.csv` są spójne z zapisanymi fitami, liczby Macieja nie.

### 4.1 J1819+1305 — osobny problem, niezależny od danych Macieja

Rekonstrukcja daje separację **17.403 ± 0.134** przy pozycjach składowych
**169.74 / 175.30 / 187.14**. Wartość **13.3531** wpisana w `separations_merged.csv`
(źródło `przeglad`, uwaga: „wsad 2026-09-17, przegląd śladów mu: czyste po odrzuceniu pików")
**nie wychodzi z żadnej pary tych składowych**. Separacja rzędu 13.4 jest osiągalna wyłącznie
przy skrajnym podzbiorze pulsów — odpowiada minimum zakresu, w jakim wędruje G3
(`G3_high` schodzi do ~521.7 binu w pojedynczych, obarczonych dużym błędem pulsach).
Podejrzenie: cięcie „po odrzuceniu pików" usunęło zbyt wiele pulsów.

Wartość Macieja (18.50) leży 2.7σ od rekonstrukcji, czyli bliżej niż `merged`.

**Wymaga oddzielnej weryfikacji.**

---

## 5. Profile średnie — pozostałe trzy przypadki

Profile uzyskane przez `pdv -FTt` z `pulsar.low` / `pulsar.high` (znormalizowane do maksimum,
mediana odjęta). Wykres zbiorczy: `~/output/claude/cmp_maciej_profile.png` — czerwone ciągłe
to pozycje Macieja, zielone przerywane z `separations.csv`, fioletowe kropkowane to
`bin_st` / `bin_end` z `params.json`.

- **J1733-3716 — rozstrzygnięte na korzyść `separations.csv`.** Drugi (szeroki) komponent ma
  szczyt w 201.5° (low) i 200.4° (high). Wartość 201.87 trafia w maksimum, wartość Macieja
  196.69 leży na zboczu narastającym. Różnica −5.19° jest realnym błędem dopasowania.
- **J1714-1054 — nierozstrzygalne.** Drugi komponent jest słaby (amplituda względna 0.31 w low,
  0.17 w high), szczyty w 185.6° i 184.6°. Obie wartości (184.86 vs 185.47) mieszczą się
  w komponencie. Deklarowane błędy 0.07–0.12° wyglądają na zaniżone.
- **J1901+0716 — nierozstrzygalne z profilu średniego.** Sporny komponent to szerokie ramię
  bez własnego maksimum; brak `component_offsets.txt` dla tego pulsara.

---

## 6. Wnioski i rekomendacje

1. **Nie scalać danych Macieja do `separations_merged.csv` jako pomiarów zastępujących.**
   W każdym rozstrzygalnym przypadku wypadają gorzej od wartości już zapisanych,
   a w 32/47 przypadków są po prostu tą samą liczbą.
2. **Użyć ich jako flagi jakości.** Rozbieżne pulsary pokrywają się z przypadkami, w których
   dopasowanie jednego skrajnego komponentu jest niestabilne — to użyteczna, niezależna
   wskazówka, gdzie błędy w `separations.csv` są zaniżone.
3. **Zweryfikować J1819+1305 w `separations_merged.csv`** — wartość 13.3531 jest
   najprawdopodobniej błędna, niezależnie od danych Macieja.
4. **Brak kolumn `dsep` / `dsep err`** w pliku Macieja oznacza, że nie da się go użyć
   bezpośrednio w analizie zwężenia profili (`Plot.ppdot_separations`, `_read_separations`)
   bez uzupełnienia.

### Stan

`input/separations_maciej.csv` leży w repozytorium jako materiał referencyjny.
**Nie jest scalony** z żadnym z pozostałych plików; żaden istniejący plik nie był modyfikowany.
