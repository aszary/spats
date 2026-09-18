# Zwężenie profili pulsarów z częstotliwością — analiza MeerKAT

**Stan na 2026-09-17.** Pomiar zmiany separacji skrajnych składowych profilu między dolnym
a górnym podpasmem obserwacji MeerKAT, dla próbki pulsarów wykazujących dryf podpulsów.

Repozytorium: `github.com/aszary/spats`, gałąź `claude`.
Dziennik roboczy: `docs/separations_analysis_log.md` (= `~/claude/work/NOTES.md`, symlink).

---

## 1. Dane

### 1.1 Próbka macierzysta

533 pulsary z listy `input/pulsars.txt` — pulsary z wykrytym dryfem podpulsów
(Song et al. 2023). Dla każdego istnieje katalog `~/output/claude/<PSR>_16/`
(w kontenerze `/home/psr/output/<PSR>_16/`).

Obserwacje: **MeerKAT**, odbiornik L-band, PSRFITS z rozdzielczością pojedynczych pulsów.
Typowy plik: 1024 biny fazy na obrót, ~1000–1100 sub-integracji (pojedynczych pulsów),
długość obserwacji rzędu 150 s.

### 1.2 Podział na częstotliwości — skąd biorą się „1023" i „1523 MHz"

Surowe dane są kanalizowane na **16 kanałów częstotliwościowych** (plik `pulsar.spCf16`;
stąd sufiks `_16` w nazwie katalogu). Funkcja `Data.multifrequency_split` tnie je na trzy
podpasma przez zerowanie kanałów (`paz -Z`), po czym scala każde podpasmo w jeden kanał
(`pam -F`):

| podpasmo | zachowane kanały | plik |
|---|---|---|
| **low** | 0–2 (trzy najniższe) | `pulsar.low` |
| **mid** | 7–8 | `pulsar.mid` |
| **high** | 13–15 (trzy najwyższe) | `pulsar.high` |

Dla typowej konfiguracji (środek 1283.6 MHz, pasmo 642 MHz, kanał 40.1 MHz) daje to
podpasmo dolne o środku ~1023 MHz i górne ~1541 MHz. Stąd etykiety w kodzie.

**Uwaga ważna dla interpretacji.** Środki podpasm **nie są stałe w całej próbce** —
zależą od tego, ile pasma padło ofiarą RFI i jak wyglądała konfiguracja obserwacji.
Przeskanowane nagłówki 528 katalogów:

- `ν_low`: mediana 1022 MHz, zakres 603–1286 MHz
- `ν_high`: mediana 1541 MHz, zakres 1037–1658 MHz
- **stosunek `ν_high/ν_low`: mediana 1.507, zakres 1.00–1.75**

W próbce, która ostatecznie dała pomiary (195 pulsarów), stosunek mieści się w
1.363–1.745 przy medianie 1.507 — czyli żaden przypadek patologiczny tam nie trafił.
Etykieta „1523/1023 = 1.489", używana w kodzie na sztywno, zaniża prawdziwy stosunek
o ~1%, co przekłada się na ~5% błędu w wyprowadzonym wykładniku widmowym.
W liczbach podanych niżej wykładnik liczony jest **z częstotliwości każdego pulsara osobno**.

### 1.3 Przetwarzanie do p3-foldu

Dla każdego podpasma osobno:

1. **Odjęcie linii bazowej** (`pmod -debase`), przy okazji wyznaczane jest okno on-pulse
   `bin_st`–`bin_end` na podstawie pasma dolnego.
2. **LRFS / 2DFS** — wyznaczenie okresu modulacji podpulsów **P3** (w jednostkach okresu
   obrotu) wraz z błędem.
3. **Składanie z okresem P3** (p3-fold): pojedyncze pulsy są wrzucane do `p3_ybins`
   przedziałów fazy dryfu i uśredniane. Wynik to macierz `p3_ybins × 1024`, czyli
   „średni cykl dryfu" — kilka do kilkudziesięciu profili o wysokim S/N zamiast tysiąca
   zaszumionych pojedynczych pulsów.

Wersje: `p3fold_refine` (z douszczelnianiem P3) i `p3fold_norefine`. **W całej analizie
używany jest `norefine`** — dla spójności z wcześniejszą kampanią.

Parametry lądują w `params.json` (`p3`, `p3_error`, `p3_ybins`, `bin_st`, `bin_end`,
`nbin`, `nsubint`); p3-foldy zapisywane są jako ASCII w formacie PSRSalsa
(nagłówek + wiersze `sub chan bin wartość`).

**Stan danych:** 510 z 533 pulsarów ma komplet p3-foldów low + high. Etap kosztowny jest
więc policzony dla praktycznie całej próbki — ponowny pomiar offsetów nie wymaga dotykania
danych surowych.

---

## 2. Co dokładnie jest mierzone

### 2.1 Dopasowanie gaussów

Dla **każdego wiersza p3-foldu osobno** (czyli dla każdej fazy dryfu) i **osobno w każdym
podpaśmie** dopasowywanych jest `n_comp` gaussów plus stała, w oknie `bin_st`–`bin_end`
(`GaussianFit.fit_gaussians`). Składowe są następnie **sortowane po pozycji μ**
(`GaussianFit.component_offsets`), więc G1 to zawsze składowa najbardziej z lewej.

Dla i-tej składowej liczony jest offset między pasmami:

```
offset_i = μ_i(high) − μ_i(low)          [stopnie długości]
```

a jej pozycja to średnia z obu pasm, `μ_i = (μ_i(high) + μ_i(low)) / 2` — czyli wartość
środkowopasmowa.

### 2.2 Separacja skrajnych składowych i wielkość docelowa

Bierzemy **wyłącznie składową pierwszą i ostatnią**:

```
W      = μ_last − μ_first                       separacja mid-band [stopnie]
Δsep   = offset_last − offset_first  =  W(high) − W(low)
ΔW/W   = Δsep / W                               wielkość bezwymiarowa
```

**Dlaczego różnica, a nie same offsety.** Błąd DM albo złe wyrównanie pasm przesuwa
*cały* profil o tę samą wartość; w różnicy skrajnych składowych to się kasuje
(Hassall i in. 2012). Δsep jest więc odporne na resztkowy błąd dedyspersji, czego
pojedynczy offset nie jest.

**Dlaczego dzielimy przez W.** Δsep w stopniach silnie koreluje z szerokością profilu
(ρ_S = +0.62, p = 9×10⁻¹¹) — wykres kolorowany Δsep pokazuje w dużej mierze to, jak szeroki
jest profil, a nie fizykę. Po normalizacji ta zależność znika (ρ_S = −0.14, p = 0.21).
**ΔW/W jest wielkością porównywalną między pulsarami o różnej szerokości profilu**, Δsep nie.

Ilustracja: J1834-0426 miał czwarte co do wielkości poszerzenie w stopniach (+1.313°),
ale przy W = 95.9° to zaledwie +1.4% szerokości. J1224-6407, siódmy w stopniach (+0.922°),
przy W = 5.1° daje +18.0% — frakcyjnie jest pierwszy.

**Znak.** ΔW/W < 0 oznacza profil **węższy** przy wyższej częstotliwości, czyli zwykłe
radius-to-frequency mapping. ΔW/W > 0 to poszerzenie.

### 2.3 Uśrednianie po pulsach

Dla każdej składowej liczona jest średnia ważona offsetu po wszystkich zachowanych
wierszach p3-foldu, wagi `1/σ²`. Błąd kwotowany dla separacji to **σ_ext** (rozdmuchany
rozrzutem między pulsami), nie σ_int — bo to rozrzut profil-do-profilu jest realnym
ograniczeniem. Wyniki trafiają do `input/separations.csv`
(kolumny: `lon i`, `lon i err`, `sep`, `sep err`, `dsep`, `dsep err`).

Wielkość `n_comp` **nie jest wyznaczana automatycznie** — patrz §3.3.

---

## 3. Przebieg prac i kryteria odrzucania

### 3.1 Kampania pierwsza (lista `separations_todo.csv`)

106 pulsarów, ręczny przegląd wspierany kryterium liczbowym. Zapisanych **91**,
odrzuconych 15. Główne tryby awarii: systematyczny rozjazd low↔high we wszystkich pulsach,
niestabilne etykietowanie składowych, martwe wiersze p3-foldu, za niskie S/N.

Przy okazji poprawiono trzy rzeczy:

- **Luka w `_offset_summary`**: pojedynczy pomiar z `err = NaN` lub `0` zamieniał w NaN
  całą sumę ważoną, unieważniając cały pulsar. Naprawione filtrem per składowa.
- **Δsep było liczone i wyrzucane** — plik zapisywał tylko `sep`. Dodane kolumny
  `dsep`, `dsep err`, a 91 istniejących wierszy uzupełnione z logów (walidacja: `sep`
  z logu zgadza się z wierszem CSV co do 1×10⁻⁴ stopnia dla 91/91).
- **21 wartości w `offsets.csv` okazało się nieaktualnych**, w tym artefakty:
  J1733−3716 −12.208° → −1.827°, J1714−1054 −2.241° → −0.944°, oraz zmiany znaku
  J1824−0127 i J1825−0935.

### 3.2 Kampania druga (wsad na 380 pulsarów)

Reszta próbki macierzystej. Skrypt `~/claude/work/scripts/batch_run.jl` woła
`Data.Plot.analyse_p3folds4_agent` (wariant **nieinteraktywny**; zwykły `analyse_p3folds4`
czeka na klawisz i zawiesiłby sesję wsadową). Czas: 380 pulsarów w ~20 minut.

Bramki automatyczne:

| bramka | próg | uzasadnienie |
|---|---|---|
| odrzucenie pulsów odstających | `median + 4·MAD` na max\|Δμ\| | realny offset bywa duży (J1843−0459 ma 11 binów i to sygnał), więc próg musi być adaptacyjny |
| rozjazd systematyczny | mediana \|Δμ\| > 10 binów | kryterium adaptacyjne nie łapie awarii dotyczącej *wszystkich* pulsów; zdrowe pulsary mają 0.2–4.9 bina, zepsute 16–23 |
| za mało pulsów | < 5 zachowanych | |
| składowe zlane | W < 1.0° lub W/σ_W < 5 | w 91 zweryfikowanych minimum to W = 2.61° przy 8.8σ |
| ekstremum | \|ΔW/W\| > 0.4 | w 91 zweryfikowanych maksimum to 0.32 |

Wynik: **149 pomiarów ze 380 (39%)**.

| status | liczba |
|---|---|
| rozjazd systematyczny low↔high | 184 |
| ok, skierowane do przeglądu | 102 |
| ok, czyste automatycznie | 47 |
| za mało pulsów | 26 |
| składowe zlane | 21 |

### 3.3 Dlaczego `n_comp` jest ustawione na 2 dla wszystkich

Liczba składowych nie jest zapisana w `params.json`, a funkcja dopasowująca wymaga jej
jako argumentu. Sprawdzono dwie metody automatycznego wyznaczania na 105 pulsarach
o znanym `n_comp`:

| metoda | trafność |
|---|---|
| **stała n = 2 (linia bazowa)** | **79%** |
| licznik pików na profilu scalonym | 61% |
| kryterium BIC na dopasowaniu 1–4 gaussów | 58% |
| konsensus obu, uczciwa walidacja 2-fold | 79% |

Żadna nie bije linii bazowej (McNemar p = 0.42). Przyczyna: `n_comp` to liczba składowych
**dających się śledzić w p3-foldzie**, a profil scalony wyrzuca dokładnie tę informację
(dryf), która ją definiuje. Rozkład prawdy: 83× n=2, 21× n=3, 1× n=4.

Wniosek operacyjny: selektor jest niepotrzebny, bo błędy są wykrywalne dalej —
niedopasowane `n_comp` objawia się jako zlane składowe albo wysokie χ².
Cena: pulsary faktycznie trójskładnikowe są mierzone jako dwuskładnikowe, co dla Δsep
(liczonego i tak ze składowych skrajnych) jest mniej szkodliwe, niż się wydaje,
ale nie jest neutralne.

### 3.4 Przegląd ręczny 123 pulsarów

Zamiast oglądać ~1500 obrazków per-puls, z logu odtworzono ślady μ każdej składowej
pulse po pulsie (składowe posortowane po pozycji, tak jak robi to kod) i narysowano jeden
panel na pulsara, w arkuszach po 12: `~/output/claude/review_sheets/sheet_01..11.png`.
Skrypt: `scripts/review_sheets.py`, werdykty z uzasadnieniem: `~/claude/work/review_verdicts.csv`.

**Werdykt: 57 zachowanych, 66 odrzuconych.** Dominujące powody odrzucenia:

- bistabilne etykietowanie — ślad μ przeskakuje między dwiema pozycjami,
- pasmo górne rozrzucone przy stabilnym dolnym,
- ślady zbiegające się lub pokrywające (dwa gaussy na jednej składowej),
- trwały rozjazd low↔high w jednej tylko składowej.

**Razem: 104 nowe pomiary + 91 starych = 195.**

---

## 4. Wyniki

Tabela finalna: `~/claude/work/separations_final_merged.csv` (195 wierszy, kolumna `zrodlo`).
**Nie jest jeszcze scalona** do `input/separations.csv`.

### 4.1 Zwężenie jest efektem populacyjnym

| wielkość | wartość |
|---|---|
| N | 195 |
| mediana ΔW/W | **−0.024**, 68% CI [−0.028, −0.016] |
| pulsary zwężające profil | **125 z 195 (64%)**, test dwumianowy **p = 1.0×10⁻⁴** |
| powyżej 3σ | 60 zwężeń, 22 poszerzenia |
| rozrzut próbki σ | 0.107 |
| mediana błędu pomiaru | 0.021 (5× mniej niż rozrzut) |

To jest wynik odporny: przetrwał korektę 21 wartości, powiększenie próbki z 91 do 195
i ręczny przegląd.

**Wykładnik widmowy** `W ∝ ν^-a`, liczony poprawnie dla wielkości środkowopasmowej
(`W_high/W_low = (1 + frac/2)/(1 − frac/2)`) i z częstotliwościami każdego pulsara osobno:

> **a = +0.058**, 68% CI [+0.036, +0.063]; a > 0 dla 125 z 195 pulsarów.

To wyraźnie mniej niż wartości 0.2–0.3 zwykle cytowane w literaturze. Zastrzeżenie:
tutejsze W to odległość skrajnych dopasowanych gaussów, a nie szerokość na 10% maksimum,
i pasmo obejmuje zaledwie czynnik ~1.5 w częstotliwości — dźwignia jest mała, więc
porównanie z literaturą wymaga ostrożności.

### 4.2 Podpróbki różnią się — i to jest istotne

| podpróbka | N | mediana ΔW/W | zwężeń |
|---|---|---|---|
| stare 91 (z `todo.csv`) | 91 | −0.045 | 73% |
| nowe, auto-czyste | 47 | −0.011 | 60% |
| nowe, zachowane w przeglądzie | 57 | −0.006 | 54% |

Nowe pulsary zwężają się **słabiej** (KS względem starych: p < 0.001). Najbardziej
prawdopodobna przyczyna: stare 91 pochodzą z `separations_todo.csv`, listy
wyselekcjonowanej z przypadków, w których offset dawał się zmierzyć i był istotny —
czyli z obciążeniem w stronę dużych |ΔW/W|. **Nowa, nieselekcjonowana próbka jest
prawdopodobnie bliższa prawdziwemu rozkładowi populacji**, a to znaczy, że efekt jest
słabszy, niż sugerowała pierwsza kampania.

### 4.3 Ruch obu krawędzi profilu — to nie jest błąd dedyspersji

Rozbicie Δsep na ruch każdej krawędzi osobno (na 90 pulsarach pierwszej kampanii):

| kanał | mediana | znaki | p |
|---|---|---|---|
| składowa wiodąca | +0.184° | 63/90 dodatnich | 1.9×10⁻⁴ |
| składowa tylna | −0.065° | 59/90 ujemnych | 4.2×10⁻³ |
| centroid (wspólne przesunięcie) | +0.069° | 54/90 | 0.073 — nieistotne |

Obie krawędzie idą do środka **niezależnie**, wspólne przesunięcie jest zgodne z zerem,
a zwężenie jest symetryczne (Wilcoxon na |lead| − |trail|, p = 0.57). Gdyby przyczyną był
błąd DM albo złe wyrównanie pasm, obie składowe przesunęłyby się w tę samą stronę.

### 4.4 Zależność od parametrów pulsara — wynik niepewny

Testowano ΔW/W przeciw log P, log Ṗ, log B_d, log τ_c, log Ė, log B_lc oraz log W.
Wagi `1/σ²` są tu bezużyteczne (χ²/dof wokół średniej ważonej ≈ 64), więc użyto korelacji
rangowej Spearmana i bootstrapu.

Jedyny sygnał pojawia się w Ė / τ_c / Ṗ — to jedna i ta sama zależność, bo te wielkości
są współliniowe. W P i B_d nie ma nic.

**Historia tego wyniku, w kolejności:**

| etap | próbka | N | ρ_S | p |
|---|---|---|---|---|
| pierwsza kampania | stare, zweryfikowane | 91 | +0.257 | 0.014 |
| po wsadzie, przed przeglądem | nowe, tylko auto-czyste | 46 | +0.350 | 0.017 |
| **po pełnym przeglądzie** | **nowe, rozłączne ze starymi** | **103** | **+0.158** | **0.110** |
| łącznie po przeglądzie | wszystko | 194 | +0.207 | 0.0037 |

**Interpretacja.** Rozłączna, niezależna próbka **nie potwierdza** trendu na poziomie
p < 0.05 — ale go też nie obala. 95% przedział ufności dla nowej próbki to
**ρ ∈ [−0.04, +0.34]**, czyli obejmuje zarówno zero, jak i wartość +0.25 ze starej próbki.
Przy N = 103 i prawdziwym ρ = 0.25 szansa na wykrycie wynosiła 72%, więc wynik p = 0.11
zdarzyłby się w 28% przypadków nawet przy w pełni realnym trendzie.

Łączne p = 0.0037 przy N = 194 jest najlepszym dostępnym oszacowaniem, ale **nie jest
niezależnym potwierdzeniem** — jest zdominowane przez tę samą próbkę, z której hipoteza
się wzięła.

Niewiadoma, której nie udało się zamknąć: podpróbka 46 auto-czystych dała ρ = 0.350,
a po dołożeniu 57 przejrzanych spadło do 0.158. Albo kryterium auto-czystości selekcjonuje
coś, co korelację zawyża, albo to fluktuacja (przy N = 46 przedział ufności był szeroki:
[+0.07, +0.58]).

**Wielkość efektu, gdyby był realny:** mediana ΔW/W −0.062 w dolnym kwartylu Ė wobec
−0.035 w górnym — 2.7 punktu procentowego szerokości profilu, czyli ~28% rozrzutu próbki
(Cohen d = 0.39, rozkłady nakładają się w 85%). Z pojedynczego pulsara nie da się
powiedzieć, z którego kwartyla pochodzi.

### 4.5 Test hipotezy geometrycznej

Gdyby zwężenie było uniwersalne (to samo ΔW/W dla wszystkich), a szerokość profilu W
zależała od Ė, to trend w Δsep *w stopniach* pojawiłby się jako artefakt. Sprawdzone:
**W rzeczywiście koreluje z Ė przy ustalonym P** (ρ_S = −0.355, p = 0.001; wyższe Ė →
węższy profil). Ale ponieważ trend mierzymy w wielkości **już podzielonej przez W**,
nie może być artefaktem tej korelacji. To zastrzeżenie zostaje zamknięte —
w odróżnieniu od §4.4, które zostaje otwarte.

---

## 5. Czego ta analiza nie rozstrzyga

1. **Rozrzut jest 5× większy od błędów pomiarowych** (σ = 0.107 wobec mediany błędu 0.021).
   ΔW/W rządzi więc coś spoza płaszczyzny P–Ṗ. Najbardziej naturalny kandydat to geometria:
   kąt nachylenia α i parametr zderzenia β decydują, jak szerokość stożka przekłada się
   na szerokość profilu. W `~/output/claude/` są katalogi `*_rvm` — to najbardziej
   obiecujący następny krok.
2. **184 pulsary odrzucone jako „rozjazd systematyczny"** to 48% wsadu. Jeśli choć część
   dałoby się uratować lepszą metodą dopasowania (np. parowaniem składowych po pozycji,
   a nie przez sortowanie), próbka urosłaby do ~300 i pytanie z §4.4 byłoby rozstrzygalne
   (przy N = 300 i ρ = 0.25 moc wynosi 99%).
3. **Próbka macierzysta jest wyselekcjonowana po dryfie podpulsów.** Wnioski dotyczą
   pulsarów dryfujących, nie pulsarów w ogóle. Powiększenie statystyki tego nie zmieni.
4. **Δsep bierze tylko składowe skrajne**, więc pulsary o trzech i czterech składowych
   oddają część informacji o wnętrzu profilu.
5. **J1717−3425 (ΔW/W = −0.742)** przeszedł jako auto-czysty, bo bramka `MAX_FRAC`
   powstała już po przebiegu. Nie był oglądany. Bez niego: nowe p = 0.064, łącznie p = 0.0020.
6. **Oceny jakości w `offsets.csv` opisują stare dopasowania.** J1714−1054 wypada
   z wykresów przez ocenę 3, choć jego poprawiony wynik ma χ² = 0.91.

---

## 6. Gdzie co leży

| co | gdzie |
|---|---|
| wyniki pierwszej kampanii | `input/separations.csv` (91 wierszy, w repo) |
| wyniki wsadu | `~/claude/work/separations_batch_full.csv` (149) |
| tabela finalna | `~/claude/work/separations_final_merged.csv` (195) |
| werdykty przeglądu | `~/claude/work/review_verdicts.csv` (123) |
| arkusze przeglądowe | `~/output/claude/review_sheets/sheet_01..11.png` |
| diagram P–Ṗ (ΔW/W) | `~/output/claude/ppdot_separations.png/.pdf` |
| ΔW/W vs parametry | `~/output/claude/frac_vs_params.png/.pdf` |
| logi | `~/claude/work/logs/` |
| skrypty | `~/claude/work/scripts/` |
| dziennik chronologiczny | `docs/separations_analysis_log.md` |

Kod w repo (gałąź `claude`): `Plot.ppdot_separations`, `Plot._read_separations`,
kwarg `plots` w `Plot.analyse_p3folds4_agent`, kolumny `dsep` w `Plot._offset_summary`,
liniowa skala koloru w `Plot._ppdot`.

---

## 7. Podsumowanie w trzech zdaniach

Profile pulsarów dryfujących **zwężają się z częstotliwością** — 125 ze 195 pulsarów,
p = 10⁻⁴ — a zwężenie jest symetryczne i nie daje się wytłumaczyć błędem dedyspersji,
bo obie krawędzie profilu idą do środka niezależnie, przy wspólnym przesunięciu zgodnym
z zerem. Efekt jest jednak **słaby**: mediana ΔW/W = −0.024, co odpowiada wykładnikowi
widmowemu a ≈ 0.06, kilkakrotnie poniżej wartości cytowanych w literaturze, a rozrzut
między pulsarami pięciokrotnie przewyższa błędy pomiarowe. **Sugerowana zależność siły
zwężenia od spin-down luminosity nie została potwierdzona na niezależnej próbce**
(ρ_S = +0.158, p = 0.11 przy N = 103) i pozostaje otwarta — do rozstrzygnięcia brakuje
około stu pomiarów, które najprawdopodobniej dałoby się odzyskać z 184 pulsarów
odrzuconych dziś jako „rozjazd systematyczny".
