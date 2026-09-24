# Test „travel": dryf podpulsów czy modulacja amplitudowa

**Stan na 2026-09-24.** Metoda rozstrzygania, czy wzór podpulsów **przemieszcza się** w długości
(dryf), czy tylko **jaśnieje i gaśnie w miejscu** (P3-only) — niezależna od kryterium
Song et al. (2023) i wolna od jego głównego obciążenia.

Repozytorium: `github.com/aszary/spats`, gałąź `claude`.
Kod: `modules/travel.jl` (moduł `Travel`), `Plot.travel`, `SpaTs.travel_test`.
Skrypty: `~/claude/work/scripts/travel_*.jl`, `check_onpulse.jl`.
Dziennik: `docs/separations_analysis_log.md` (= `~/claude/work/NOTES.md`, symlink).

---

> **UWAGA (2026-09-24, późno): model zerowy (§7.1) nie jest skalibrowany dla nieseparowalnej
> zmienności impuls-do-impulsu** (jitter, podpulsy w losowych pozycjach). Pola bez żadnego ruchu
> dają T/T_inc do tysięcy σ przy sile modulacji spotykanej w danych, a kontrola off-pulse tego nie
> łapie. Frakcje detekcji z §5.3 nie są więc dowodem uporządkowania czasowego. Skan spójności (§7.4)
> jest na ten efekt odporny. Szczegóły: dziennik, wpis 2026-09-24 (cd.).
>
> **Naprawa (2026-09-24, cd. 2): statystyka krzyżowa między blokami `T_cv` z nullem z randomizacji
> znaków bloków** — bez modelu zmienności, skalibrowana na wszystkich syntetykach bez ruchu, czułość na
> sztywny dryf porównywalna z T. Przebieg v3: trwałe uporządkowanie u **271/406 dryferów (67%)** i
> **12/107 P3-only (11%)**. Szczegóły, walidacja i tabele: dziennik, wpis 2026-09-24 (cd. 2).
> **To jest teraz główna statystyka detekcyjna; T i T_inc zostają jako diagnostyka.**

## 1. Problem

Song et al. (2023) klasyfikują cechę w 2DFS jako dryf, jeśli jej **centroida mocy** jest istotnie
przesunięta względem osi 1/P₂ = 0. Słabość jest w estymatorze, nie w idei: stochastyczna zmienność
kształtu pulsu wrzuca moc wzdłuż osi 1/P₂ = 0, a gdy ta moc jest asymetryczna — zwykle jest —
centroida przesuwa się od zera i przy niedoszacowanym błędzie przekracza próg. Pozorny offset bierze
się z **biasu centroidy**, nie z ruchu podpulsów. Dodatkowo test istotności przez shuffle bada
hipotezę „to tylko szum", a nie właściwą „to mogłaby być czysta modulacja amplitudowa".

---

## 2. Tożsamość zerowa

Dla fluktuacji `δI(n, φ)` po odjęciu profilu statycznego:

```
K(Δ, τ) = Σ_{n,φ}  δI(n,φ) · δI(n+τ, φ+Δ)
```

**Każda** modulacja amplitudowa jest *separowalna*: `δI(n,φ) = a(φ)·w(n)` z rzeczywistym (także
ujemnym) a(φ). Wtedy `K = [Σ_φ a(φ)a(φ+Δ)]·[Σ_n w(n)w(n+τ)]`, a drugi czynnik jest **dokładnie
parzysty w τ**:

```
Ĉ_w(−τ) = Σ_n w(n)·w(n−τ)  = [m = n−τ] =  Σ_m w(m+τ)·w(m) = Ĉ_w(τ)
```

To ta sama suma po tych samych parach — **żadnego założenia o w(n)**. Stąd

```
A(Δ, τ) = K(Δ, τ) − K(Δ, −τ)  ≡  0     przy modulacji amplitudowej
```

**tożsamościowo**, dla zrealizowanych danych. Bez zakładania separowalności
`A = Σ_φ [C_{φ,φ+Δ}(τ) − C_{φ+Δ,φ}(τ)]`, czyli netto „czy długość φ wyprzedza φ+Δ?". Zeruje się
dla składowych w antyfazie, niezależnych modulacji o różnych P3, impulsów zapowanych i po filtrze
wysokopasmowym.

---

## 3. Związek z 2DFS

To nie jest nowy kanał informacji — transformata K to widmo mocy 2D, a A to jego asymetria względem
k_φ → −k_φ, formalnie ten sam kanał co kryterium 2DFS. Różnica jest w estymatorze: część symetryczna
jest **rzutowana do zera algebraicznie** zamiast uśredniana w centroidzie, używana jest cała
płaszczyzna (Δ, τ) zamiast ręcznego prostokąta, a model zerowy to modulacja amplitudowa, nie szum.

---

## 4. Dwa różne pytania

Rozdzielenie ich jest kluczowe i długo je myliłem.

| | pytanie | statystyka | założenia |
|---|---|---|---|
| **A** | czy jest **jakikolwiek** ruch? (czyli: czy to na pewno nie jest czysta modulacja amplitudowa) | `T`, `T_inc` | brak — tożsamość z §2 |
| **B** | czy modulacja to **zasadniczo** ruch, o charakterze sztywnej translacji? | `ρ` | wymaga mierzalnej geometrii (P₂ wewnątrz profilu) **oraz dostatecznie stabilnego P₃** (§10.3) |

Klasyfikacja Song et al. jest jakościowa i binarna, więc odpowiada jej **pytanie A**. `ρ` to
charakterystyka dodatkowa, cenna tam, gdzie da się ją policzyć, ale **nie jest produktem głównym** —
przez pewien czas błędnie ją za taki uważałem.

Różnica ma znaczenie praktyczne: `T` jest statystyką detekcyjną, więc jej istotność rośnie z S/N.
Jasny pulsar, w którym wędruje 1% modulacji, da ogromne `T`. „T = 300σ" znaczy „ruch jest
wykrywalny", a nie „ten pulsar jest dryferem".

---

## 5. Produkt główny: T, T_inc i kontrole

### 5.1 Statystyki

**`T = Σ A(Δ,τ)²`** — omnibus, bez założeń o geometrii.

**`T_inc = Σ_b Σ A_b²`** — suma niekoherentna po blokach impulsów. `T` liczy się na jednej globalnej
mapie, więc dryfer spędzający tyle samo czasu w każdym kierunku kasuje się; kalibracja: zbalansowany
reverser daje **T/T_inc = 3.4·10⁻⁶**. `T_inc` płaci wyższym progiem szumu, więc dla dryfu stałego
jest mniej czuła. Raportować obie.

### 5.2 Kontrole, bez których wynik nie znaczy nic

**Off-pulse.** Ta sama statystyka na pasku bez sygnału — musi wyjść zgodna z zerem. Jedyna realna
podatność metody to szum niesymetryczny w czasie (dryf wzmocnienia, RFI, zła linia bazowa).
W pełnym przebiegu odrzuciła 3 z 515.

**Spójność blokowa, skanowana po długości bloku** (`block_consistency_scan`, §7.4). Leave-one-out;
projekcja na mapę globalną zawierałaby człon własny i dawała ~1/√n dla szumu. Pojedyncza wartość
**nie wystarcza** — dudnienie dwóch bliskich okresów udaje spójność 0.95, jeśli bloki są dłuższe od
okresu dudnienia. Raportowaną wielkością jest **`min(cons)`** po całym skanie, czytana **wyłącznie
łącznie z siłą detekcji** (§7.4).

### 5.3 Wynik na pełnej próbce

Przebieg v2 (`travel_batch_v2.csv`, 2026-09-24): `max_dphi = W₃σ ÷ 2` (§8.3) i skan spójności
(§7.4). W nawiasie pierwszy przebieg z `max_dphi = M ÷ 2`.

| etykieta Song+23 | w analizie | z detekcją ≥5σ |
|---|---|---|
| drift | 406 | **368 (91%)** (330/405, 82%) |
| p3only | 107 | **79 (74%)** (57/107, 53%) |

Żaden pulsar nie stracił detekcji. Istotność wzrosła medianowo 1.8×, a nowe detekcje pochodzą
z pulsarów o najbardziej przewymiarowanym oknie. Stary zakres Δ dokładał tam biny bez emisji.
Model zerowy pozostał skalibrowany: σ kontroli off-pulse 1.03 (T) i 0.91 (T_inc).

**Tej frakcji nie wolno podawać bez spójności.** `min(cons)` wśród detekcji:

| siła detekcji | drift: n, mediana | p3only: n, mediana |
|---|---|---|
| 5–20σ | 54, 0.086 | 31, 0.018 |
| 20–100σ | 101, 0.231 | 23, 0.035 |
| >100σ | 213, **0.419** | 25, **0.088** |

Przy tej samej sile detekcji P3-only mają spójność kilkakrotnie niższą, a około ⅓ ma ujemną.
Poniżej 0.2 jest 71/79 detekcji p3only i 129/368 detekcji drift. Odczyt: 74% P3-only nie jest
czystą modulacją amplitudową, ale u zdecydowanej większości uporządkowanie czasowe **nie jest
trwałe**. Nie ma tu stałego wyprzedzania, jakie widać u dryferów. Przejście przez zero przy
detekcji ≥20σ: 0/314 drift, 1/48 p3only (J1001-5939). Próg na `min(cons)` zależny od S/N nie
jest jeszcze wyznaczony (§11).

---

## 6. Charakterystyka dodatkowa: ρ

### 6.1 Konstrukcja

Dla wzoru `cos(2π(φ/P₂ − n/P₃))` rozwinięcie cosinusa różnicy daje dwie połowy **o równych
współczynnikach**:

```
cos(2π(Δ/P₂ − τ/P₃)) = cos(2πΔ/P₂)cos(2πτ/P₃)  +  sin(2πΔ/P₂)sin(2πτ/P₃)
                       └── parzysta w Δ ──┘        └── nieparzysta w Δ ──┘
```

Nieparzystą mierzy `A = antisym_map(K)`, parzystą `E = sym_map(K)`. Modulacja amplitudowa, będąc
separowalną, wkłada wszystko w parzystą. Rzutując na oba szablony (z trójkątnym taperem
`(1−τ/N)(1−Δ/M)` korelacji liniowej — jego pominięcie zaniża projekcję o 23% przy N=600, M=40):

```
ρ = √2 · frac_odd / √(frac_odd² + frac_even²)
```

ρ = 1 dla sztywnego dryfu, 0 dla modulacji amplitudowej, ograniczone przez √2, zawsze zdefiniowane.
Klasyfikatorem jest |ρ|; znak niesie kierunek dryfu. (Wcześniejszy `R = odd/even` to tangens tego
samego kąta i ma biegun — stąd wartości 4.03, NaN-y i potrzeba dwóch progów jakości. Zostaje w
wyniku, bo relacje wyżej są w nim sformułowane.)

### 6.2 Po co iloraz — cyrkularność

Ocena samego `frac_odd` wymagałaby wiedzy, ile zwraca prawdziwy dryfer, a oczywista próbka
referencyjna (pulsary oznaczone `drift`) to zbiór podejrzany o skażenie. ρ bierze miarę **z tego
samego pulsara**: siła modulacji, harmoniczne, zanik koherencji i taper mnożą obie połowy tak samo
i kasują się. Wśród dryferów `frac_odd` rozciąga się na 25×, a ρ na 2×.

### 6.3 Geometria

**|P₂| dopasowywane** skanem maksymalizującym projekcję **parzystą** — parzysta mierzy koherentną
modulację niezależnie od ruchu, więc wybór geometrii nie może wyprodukować ruchu. Znak P₂ wychodzi
ze znaku projekcji nieparzystej.

**P₃ NIE jest dopasowywane** — brane z pomiaru LRFS. Swobodny fit ucieka w róg dużych P₂/P₃
(J2053-7200: 63 zamiast 3.06).

**`max_dphi` powinno pochodzić z zasięgu emisji (W₃σ), nie z zadeklarowanego okna** — patrz §8.3.

Poza geometrią ρ zakłada też **dostatecznie stabilne P₃** — tolerancja ~±20%, powyżej tego wynik
przestaje być interpretowalny jako miara sztywności (§10.3).

### 6.4 Wynik tam, gdzie geometria jest mierzalna

Reżim P₂fit ≤ M/2, 179 pulsarów:

| | n | mediana ρ | kwartyle |
|---|---|---|---|
| **drift** | 172 | **1.011** | 0.830–1.087 |
| p3only | 7 | 0.324 | 0.092–0.926 |

**Mediana 1.011 na 172 niezależnych pulsarach przy przewidywaniu teorii dokładnie 1 i zerowej
liczbie parametrów swobodnych.** Sama geometria też rozróżnia: mediana P₂fit/M = **0.39** dla drift
(63% poniżej 1) wobec **1.39** dla p3only (21% poniżej 1).

ρ **nie zależy od S/N** — pozorna zależność w próbce zbiorczej była paradoksem Simpsona: w obu
grupach geometrii ρ jest płaskie (0.912/1.006/1.032 oraz 0.342/0.326/0.318 dla rosnącego S/N),
a S/N steruje tylko tym, do której grupy pulsar trafia (40%/54%/61% w grupie P₂/M < 0.5).

---

## 7. Model zerowy

### 7.1 Konstrukcja

Pod H₀ sygnał wnosi do A dokładnie zero, więc cała wariancja pochodzi z członów szumowych. Surogat =
wiodący mod SVD (czyli dokładnie separowalny) przeskalowany do odszumionej mocy, plus szum
**bootstrapowany z własnego off-pulse'u** pulsara (ciągły pasek tej samej szerokości, losowe
przesunięcie cykliczne w czasie; przy profilu szerszym niż najdłuższy ciągły fragment — cyklicznie
po liście binów, flaga `offpulse_wrapped`).

Surogat jest rank-1 **celowo i poprawnie**: „separowalny" znaczy dokładnie „rank 1", a pole rank-1
z definicji nie może wędrować. Reprezentuje więc ściśle tę hipotezę, którą ma reprezentować.
Sprawdzone: dla pola rank-1 wychodzi T = −0.3σ, T_inc = −1.3σ.

### 7.2 Uporządkowanie to nie zawsze dryf: dudnienie

Rozważ pulsar z dwiema nakładającymi się składowymi, z których **każda pulsuje własnym okresem**.
Nic się nie przemieszcza — to nadal modulacja amplitudowa. Ale

```
K = Ĉ_{a1}(Δ)Ĉ_{w1}(τ) + Ĉ_{a2}(Δ)Ĉ_{w2}(τ)  +  [Σ_φ a₁a₂]·[Σ_n w₁(n)w₂(n+τ)] + (sym.)
```

Dwa pierwsze człony są parzyste w τ i znikają w A. Człon skrośny zawiera **korelację wzajemną**,
która parzysta nie jest. Dwa bliskie okresy dudnią jak dwie rozstrojone struny: przez część cyklu
dudnienia jedna składowa błyska wcześniej, przez resztę druga. Zmierzone na syntetyku, **przy
całkowitym braku ruchu**:

| P₃ składowych | okres dudnienia | T | T_inc | spójność blok. (4 bloki) |
|---|---|---|---|---|
| 7.0 i 7.4 | 129 P | **37.9σ** | **73.3σ** | **0.95** |
| 7.0 i 8.0 | 56 P | 1.7σ | 3.1σ | −0.36 |
| 7.0 i 9.0 | 32 P | −0.1σ | 9.7σ | 0.76 |
| 7.0 i 13.0 | 15 P | 0.9σ | 1.0σ | 0.20 |

Groźne są **bliskie** okresy: przy dudnieniu 129 P obserwacja mieści ich tylko ~8, więc efekt się nie
uśrednia. Przy dudnieniu 15 P przechodzi 67 razy i znika.

**To nie jest fałszywy alarm statystyki, tylko prawdziwy alarm na coś innego.** W tych danych
uporządkowanie czasowe **naprawdę jest** — na odcinku krótszym niż dudnienie jedna długość dosłownie
wyprzedza drugą. T = 74σ jest liczbą poprawną, a surogat nie skłamał: pulsar z jednym zegarem
faktycznie nigdy by tyle nie wyprodukował. Błędny był krok rozumowania **„jest uporządkowanie ⇒ jest
dryf"** — uporządkowanie jest dla dryfu konieczne, ale niewystarczające.

Rozróżnienie ma konsekwencje praktyczne, bo wskazuje, gdzie naprawiać: **nie w modelu zerowym, tylko
po detekcji**.

### 7.3 Odrzucona naprawa: surogat rank-r z randomizacją faz

Naturalny odruch to poszerzyć null: zachować r modów SVD stojących ponad progiem szumu
(σ(√N+√M) dla czystego szumu) i każdemu zrandomizować fazy Fouriera, co zachowuje widmo mocy
(a więc rytm i autokorelację), a niszczy wzajemne ustawienie modów. Zaimplementowane
(`rank_r_modes`, `phase_randomize!`, przełącznik `surrogate_rank`). Wynik:

| | surogat rank-1 | surogat rank-r |
|---|---|---|
| dudnienie 7.0/7.4 (brak ruchu) | T = 38σ, T_inc = 75σ | **T = −1.2σ** |
| dudnienie 7.0/9.0 | T_inc = 8.9σ | −1.9σ |
| **kontrola: prawdziwy dryf** | **T = 19112σ** | **T = 1.9σ** |

Fałszywki znikają — ale razem z sygnałem. To nie jest niedoróbka implementacji, tylko rzecz
nieusuwalna: **dryf *jest* określoną relacją fazową między dwoma modami** (dwa mody w kwadraturze
przy tej samej częstości dają falę bieżącą, przy przesunięciu zerowym — stojącą). Randomizacja faz
losuje dokładnie tę relację, więc produkuje null „wzór może biec albo stać z równym
prawdopodobieństwem", a nie „wzór nie biegnie". Poszerzając null tak, by objął dudnienie, obejmuje
się nim **również dryf**.

**Wniosek: null zostaje rank-1.** `surrogate_rank` domyślnie 1; opcja zachowana w kodzie jako zapis
sprawdzonego i odrzuconego wariantu.

### 7.4 Właściwy dyskryminator: skan po długości bloku

Różnica nie tkwi w tym, *czy* jest uporządkowanie, tylko czy jest **trwałe**: dryf ma A wyprzedzające
B przez całą obserwację, dudnienie odwraca znak co pół okresu dudnienia. Stąd spójność blokowa
liczona przy **kilku długościach bloku**:

| przypadek | 4 bloki (250 P) | 10 bloków (100 P) |
|---|---|---|
| dudnienie 7.0/7.4 (129 P) | +0.98 +0.97 +0.97 +0.95 | **−0.12 −0.68 −0.99 −0.97 −0.05 +0.93 −0.99 −0.99 +0.90 +0.61** |
| prawdziwy dryf | +1.00 ×4 | **+1.00 ×10** |

Przy blokach dłuższych od dudnienia uśrednienie **udaje idealną spójność** — dlatego pojedyncza
wartość `block_consistency` nie wystarcza, a reguła „T_inc bez zgodności blokowej nie jest
kandydatem" jest **niewystarczająca**: najgorszy przypadek ma zgodność 0.95. Prawdziwy dryf jest
niewzruszony przy każdej długości.

**Blokada, która to uniemożliwiała, i jej usunięcie.** Liczba bloków była przycinana do
`N ÷ (4·max_lag)`, bo blok musi pomieścić opóźnienia do `max_lag`. Przy realnych danych dawało to
najwyżej 4 bloki — dokładnie reżim, w którym dudnienie udaje spójność, więc w pełnym przebiegu test
nie miał szans zadziałać. Rozwiązanie: **mapa bloku nie potrzebuje tego samego zasięgu τ co mapa
globalna**. `block_consistency_scan` bierze dla każdego podziału `lag_b = min(max_lag, L÷4)`, więc
drobniejsze podziały automatycznie używają krótszych opóźnień. Zwracane są `block_scan_nb`,
`block_scan_len`, `block_scan_lag`, `block_scan_cons` oraz **`block_scan_min`**.

Wynik na syntetyku:

| przypadek | 2 bl. | 4 bl. | 8 bl. | 16 bl. | 32 bl. |
|---|---|---|---|---|---|
| dudnienie 7.0/7.4 (129 P) | +0.95 | +0.97 | +0.97 | **−0.28** | **−0.35** |
| prawdziwy dryf | +1.00 | +1.00 | +1.00 | +1.00 | +1.00 |
| rank-1 (ścisłe H₀) | +0.22 | −0.03 | +0.15 | +0.11 | +0.04 |

i na pulsarach:

| pulsar | rola | T | **min(cons)** | przebieg |
|---|---|---|---|---|
| J0151-0635 | dryfer wzorcowy | >999σ | **0.85** | 0.91 → 0.85 |
| J0034-0721 | dryfer wzorcowy | >999σ | **0.53** | 0.71 → 0.53 |
| J0837+0610 | kand. promocji | 98σ | 0.46 | 0.86 → 0.46 |
| J0304+1932 | kand. degradacji | >999σ | **0.22** | 0.43 → 0.22 |
| J0629+2415 | P3-only, ρ = 1.36 | 123σ | **0.15** | 0.54 → 0.15 |
| J0601-0527 | P3-only, ρ = 1.37 | 100σ | **0.07** | 0.20 → 0.07 |
| J1907+0731 | P3-only, T_inc = 4.9σ | 2σ | **−0.01** | ≈0 wszędzie |

J0601-0527 i J1907+0731 mają spójność bliską zeru **przy każdej długości bloku**, mimo detekcji 100σ
w przypadku pierwszego — uporządkowanie w nich nie odtwarza się w ogóle, więc nie jest ani dryfem,
ani dudnieniem o długim okresie. Obie miały przy tym wysokie ρ (1.37 i 1.36), czyli skan mówi tu coś,
czego ρ nie mówiło.

**Dwa zastrzeżenia.**

*Spadek sam w sobie nie jest dowodem.* Prawdziwe dryfery też schodzą z długością bloku (J0034-0721:
0.71 → 0.53), bo krótsze bloki mają szumniejsze mapy. Sygnaturą dudnienia jest **przejście przez
zero** albo utrzymywanie się blisko zera, nie samo opadanie. Stąd reguła: **`min(cons)` wolno czytać
tylko łącznie z siłą detekcji** — przy T > 999σ wartość 0.22 jest znacząca, przy detekcji 10σ ta
sama liczba nie znaczy nic. (Warto odnotować, że syntetyczny dryf trzyma +1.00, a realne dryfery
0.53–0.85: realna modulacja nigdy nie jest tak koherentna jak model.)

*Krótkie dudnienia uciekają.* Dudnienie 7.0/9.0 (okres 32 P) **nie zostało złapane** — zostaje
+0.94, bo najkrótsze bloki mają 31 P, czyli tyle co samo dudnienie, a ich τ ≤ 7 już ledwie pokrywa
P3. Jego T = 1.6σ, ale T_inc = 10.5σ, więc przeszłoby jako detekcja. To pozostaje luką.

---

## 8. Dane i okna

### 8.1 Pełne pasmo

Katalogi `<PSR>_16/` zawierały tylko podpasma z `paz -Z` (zapuje wymienione kanały, zostaje reszta):
`low` = kanały 0–2, `mid` = 7–8, `high` = 13–15. Podpasmo kosztuje czynnik ~2.25 w RMS (J0601-0527:
0.0275 wobec 0.062, zgodnie z √(16/3)). Pierwszy przebieg mieszał 85 pulsarów pełnopasmowych z 430
na 3/16 pasma — niedopuszczalne dla frakcji detekcji. Pełne pasmo odtwarzane raz i **zostaje w
katalogu `_16`**: `pulsar.full` / `pulsar_full.debase.gg` / `pulsar_full_debase.txt`, ~7 s i ~72 MB
na pulsara.

Wstępne przetwarzanie: filtr wysokopasmowy biegnącą średnią, `hp_halfwin = 50` — ten sam dla każdej
długości, więc separowalność (a z nią tożsamość zerowa) zostaje nienaruszona.

### 8.2 Okna on-pulse są hojne, ale poprawne

Na 515 pulsarach: M medianowo **116 binów** wobec zasięgu emisji W₃σ = **56**, stosunek 0.49, szczyt
wewnątrz okna u **512/515**. Okna są dwukrotnie szersze niż kontur 3σ, ale to normalna praktyka —
próg 3σ obcina skrzydła — a nie błąd automatu `pmod`.

### 8.3 Mechanizm zależności od okna i jego naprawa

Okno **nie psuje ρ bezpośrednio**: przy P₂ ustalonym na sztywno ρ wynosi 0.871/0.906/0.977/0.894 dla
M = 100…250, czyli bez trendu. Rozcieńczenie szumem kasuje się w ilorazie.

Łańcuch przyczynowy jest inny:

1. `max_dphi = M/2`, więc szersze okno = dalszy przeszukiwany zakres Δ, wchodzący w obszar bez
   emisji.
2. Fit P₂ maksymalizuje projekcję parzystą, a w rogu dużych P₂ szablon parzysty jest prawie stały
   i dopasowuje się do gładkiego tła. P₂fit idzie 150 → 186 → 232 → **309** przy prawdziwym 150.
3. Zawyżone P₂ kaleczy **asymetrycznie**: przy Δ → 0 `cos(2πΔ/P₂) → 1` (pełna waga tam, gdzie sygnał
   jest), a `sin(2πΔ/P₂) → 0` (zerowa waga tam, gdzie sygnał jest). Kanał nieparzysty traci pokrycie
   z sygnałem, parzysty nie. ρ leci w dół.

**Naprawa: `max_dphi` z zasięgu emisji, nie z okna.** Zweryfikowane — P₂fit staje się idealnie
stabilne (150.2 przy M = 100, 150, 200, 250), a zjazd ρ spada z 0.87→0.65 do 0.87→0.78. Bez nowych
parametrów: W₃σ jest zmierzone dla wszystkich 515 pulsarów. **Wdrożone w przebiegu v2** — mediana
ρ drift (P₂fit ≤ M/2) 0.990 przy n = 200 (wcześniej 1.011 przy n = 172).

Poszerzanie samego okna szkodzi niezależnie (rozcieńcza sygnał, zawyża P₂fit do 309 przy M = 250,
zjada obszar off-pulse — przy M = 300 metoda przestaje działać). Zwiększanie samego zasięgu Δ też
nie pomaga. Powód zasadniczy: **informacja o okresie w długości pochodzi wyłącznie z obszaru, który
świeci.**

---

## 9. Reżimy: kiedy P₂ jest w ogóle mierzalne

Liczba jednocześnie widocznych podpulsów to ≈ **W/P₂**. Jeśli **P₂ > W**, w danym impulsie widać
najwyżej jeden podpuls, wędrujący przez profil. To jest fizycznie dopuszczalne: P₂ jest duże, gdy
iskier w karuzeli jest mało albo gdy linia widzenia przecina pierścień emisji blisko prostopadle
(wąski zakres azymutu karuzeli). Przy przejściu stycznym jest odwrotnie.

**W tym reżimie P₂ nie jest mierzalne, a jedynie ograniczone od dołu.** Nie widać dwóch sąsiednich
podpulsów naraz, więc w danych nie ma periodyczności w długości. Stąd:

- ucieczka fitu na kraniec siatki **nie jest usterką**, tylko poprawną odpowiedzią „P₂ ≥ tyle";
  wartości P₂fit powyżej M należy raportować jako **limity, nie pomiary**;
- ρ degraduje się w sposób zmierzony: 0.85 przy P₂/M = 2, 0.64 przy 4 — częściowo z powodu
  nieusuwalnego ubytku informacji (to samo ogranicza 2DFS), częściowo z powodu wyboru szablonu
  wymagającego pełnego cyklu;
- `T` i `T_inc` działają **bez zmian**, bo nie zakładają geometrii. W tym reżimie mamy więc
  detekcję ruchu, ale klasyfikacja przez ρ jest niepewna i **nie wolno z niej robić degradacji**.

W granicy dużego P₂ `sin(2πΔ/P₂) → 2πΔ/P₂`, czyli szablon nieparzysty staje się liniowy w Δ — to
ta sama informacja, którą mierzy gradient fazy w `PhaseDrift.drift_test`. Dodanie tego wariantu
granicznego jako osobnej statystyki dla reżimu P₂ > W jest naturalnym uzupełnieniem.

---

## 10. Walidacja

### 10.1 Syntetyczna (`Travel.selftest()`)

| test | wynik |
|---|---|
| FFT z paddingiem zerowym vs suma wprost | 4.3·10⁻¹⁶ |
| pole separowalne: wędrujące P3 + nulling + skok znaku a(φ) | max\|A\| = 1.1·10⁻¹⁶ × K(0,0) |
| syntetyczna fala bieżąca | P₂ = 18.0 (prawda 18), P₃ = 12.0 (prawda 12), rank1 = 1.00 |
| projekcja matched w prawdziwej geometrii | frac = 1.000 |
| ρ dla czystego dryfu | 1.000 |
| ρ dla fali bieżącej + stojącej | 0.607 przy przewidzianym analitycznie 0.600 |
| zbalansowany reverser | T/T_inc = 3.4·10⁻⁶ |

### 10.2 Krzywa odniesienia ρ(P₂/M)

| P₂/M | 0.25 | 0.5 | 1.0 | 1.5 | 2.0 | 3.0 | 4.0 |
|---|---|---|---|---|---|---|---|
| ρ | 1.00 | 1.00 | 0.98 | 0.87 | 0.85 | 0.69 | 0.64 |

Materiał odniesienia, nie kryterium stosowane przez kod.

### 10.3 ρ wobec P₃: wartość nie szkodzi, wędrówka szkodzi

**Wartość P₃ jest bez znaczenia.** Sztywny dryf przy P₂/M = 0.4 i P₃ = 2.05 / 2.5 / 3.5 / 5 / 9 daje
ρ = **0.988 / 0.988 / 0.989 / 0.989 / 0.989**. Wersja szablonowa jest na to zupełnie niewrażliwa.

**Wędrówka P₃ już nie**, i to niemonotonicznie. Sztywny dryf, P₃ średnie 5, wahające się
sinusoidalnie o okresie 300 impulsów:

| wędrówka P₃ | zakres | ρ | T |
|---|---|---|---|
| ±0 (stałe) | 5.0 | **0.989** | 8916σ |
| ±10% | 4.5–5.5 | 1.012 | 16265σ |
| ±20% | 4.0–6.0 | 1.075 | 9841σ |
| ±40% | 3.0–7.0 | **1.385** | 7922σ |
| ±60% | 2.0–8.0 | **0.510** | 6404σ |
| monotonicznie 4→6 | — | 1.038 | 10670σ |

Mechanizm widać z rozwinięcia: przy wędrującym okresie faza nagromadzona do opóźnienia τ ma rozkład,
więc `K ∝ cos(2πΔ/P₂)·⟨cos Φ(τ)⟩ + sin(2πΔ/P₂)·⟨sin Φ(τ)⟩`. Równość współczynników, na której stoi
ρ = 1, wymaga Φ(τ) = 2πτ/P₃ z **jednym** P₃; rozmycie fazy tłumi oba człony niejednakowo i miesza
kanały.

Trzy wnioski:

1. **Detekcja jest odporna** — T trzyma 6400–16000σ w każdym przypadku, zgodnie z §2. Pytanie A nie
   zakłada nic i nadal nie zakłada.
2. **ρ toleruje umiarkowaną wędrówkę**: do ±20% błąd nie przekracza 8%, a powolna monotoniczna
   zmiana P₃ jest praktycznie niewidoczna (1.038). To istotne, bo wolny dryf P₃ jest pospolity.
3. **Przy silnej wędrówce ρ psuje się w obie strony**: ±40% → 1.385, ±60% → **0.510**. Ten drugi
   przypadek jest groźny, bo silnie wędrujący dryfer zostałby odczytany jako „nie wędruje", czyli
   fałszywie zdegradowany. Zapaść przy ±60% wiąże się z tym, że P₃ schodzi tam do 2.0, gdzie czynnik
   τ szablonu nieparzystego degeneruje się (‖sin‖²/‖cos‖² = 0.160 przy P₃ = 2.05).

**Konsekwencja interpretacyjna: ρ ≠ 1 nie znaczy „to nie jest sztywny dryf", dopóki nie wiadomo, jak
stabilne jest P₃.** Obserwowany ogon 1.1–1.35 u 25 dryferów odpowiada ilościowo wędrówce rzędu
±25–40%, co dla realnych pulsarów jest typowe — hipoteza z §11 pkt 2 jest więc zgodna co do rzędu
wielkości, nie tylko co do kierunku.

Brakuje **niezależnej miary stabilności P₃**. Naturalną jest szerokość cechy f₃ w LRFS (wędrujący
okres ją poszerza; pliki `pulsar_*.debase.lrfs` są w katalogach). Wąska cecha → ρ interpretowalne
jako miara sztywności; szeroka → ρ raportować, ale nie wyciągać z niego wniosków o rygidności.
To zarazem test hipotezy o ogonie: jeśli te 25 obiektów ma systematycznie szersze cechy f₃, sprawa
jest zamknięta. `p3_error` z `params.json` się nie nadaje — wartości rzędu 0.0013 przy P₃ ~2–5 to
formalne błędy dopasowania piku, nie szerokości. **Niezmierzone.**

### 10.4 Test stabilności okna

ρ przy oknie 1.0 / 1.5 / 2.0 × W₃σ: prawdziwe dryfery są odporne (J0034-0721: 1.042/1.067/1.048;
J0151-0635: 0.970/1.011/1.006), ale obiekty o niskim ρ chwieją się nawet 40-krotnie (J1239+2453:
0.011/0.374/0.461). Walidacja populacyjna stoi, ale **pojedynczy obiekt o niskim ρ nie jest
zmierzony bez tego testu**. Uwaga: kryterium „rozrzut/σ_ρ" **nie działa** — formalne σ_ρ łapie tylko
szum surogatów i jest o rząd wielkości za małe dla jasnych (J0034-0721: rozrzut 0.025 przy
σ_ρ = 0.001) i za duże dla słabych. Liczy się rozrzut bezwzględny.

### 10.5 Wariant bez szablonu — sprawdzony i odrzucony

Skoro sztywny dryf daje obu połowom równe współczynniki, ich **normy** też powinny być równe, więc
`√(‖A‖²/‖E‖²)` byłoby klasyfikatorem bez żadnej geometrii. Na syntetyku działa i jest całkowicie
odporne na szerokość okna (0.819/0.829/0.829/0.829 tam, gdzie wersja szablonowa spada 0.872 → 0.653).
Ma też policzalną korektę na próbkowanie τ: bez niej przy P₃ = 2.05 czyta 0.407 zamiast 1, bo
‖sin‖²/‖cos‖² = 0.160; po podzieleniu przez √(‖sin‖²/‖cos‖²) wraca do 1.017–1.021 w całym zakresie
P₃ = 2–9.

**Ale na danych rzeczywistych zawodzi:**

| | mediana drift | mediana p3only | rozdzielczość |
|---|---|---|---|
| z szablonem | **1.048** | **0.361** | **2.9×** |
| bez szablonu + korekta | 0.746 | 0.599 | 1.25× |

J0034-0721 (B0031−08, podręcznikowy dryfer, T > 999σ) czyta **0.203**. Powód: `‖E‖²` zbiera **całe
tło mapy** — autokorelację profilu przy małych Δ, harmoniczne, składowe niezwiązane z dryfem —
a `‖A‖²` nie, bo A znika przy Δ = 0 z konstrukcji. Mianownik jest zawyżony o rzeczy, które z
geometrią dryfu nie mają nic wspólnego.

**To jest właściwa odpowiedź na pytanie „po co w ogóle fit P₂": szablon jest jedyną rzeczą, która
wycina z mapy część związaną z dryfem i odrzuca resztę.** Syntetyk tego nie pokazał, bo pojedyncza
sinusoida daje A i E identyczną strukturę. `rho_free` zostaje w module jako diagnostyka (ile w mapie
jest struktury poza modelem dryfu), nie jako klasyfikator.

---

## 11. Ograniczenia, w kolejności ważności

0. **Null rank-1 nie niesie nieseparowalnej zmienności on-pulse** — T rośnie ~liniowo z siłą modulacji
   dla pól bez ruchu (σ ≈ 60·mod). **Zaadresowane przez `T_cv`** (`crossblock_test`). Jego własne
   ograniczenia: istotność ograniczona liczbą bloków; długie P₃ słabo pokryte przy lag ≤ L÷4;
   reverser o losowych epizodach niewidoczny (`T_adj` go łapie, ale myli z dudnieniem).
1. **Dudnienie dwóch bliskich okresów udaje dryf** (§7.2). Null jest poprawny i T jest poprawne;
   rozdzielenie następuje po detekcji, skanem spójności po długości bloku (§7.4) — **wdrożone**,
   `block_scan_min`. Pozostałe luki: (a) dudnienia o krótkim okresie (≲ najkrótszy blok) uciekają,
   (b) `min(cons)` spada też z powodu szumu, więc jest czytelne tylko przy silnej detekcji;
   **brak progu/odniesienia zależnego od S/N** (kandydat: skan na surogatach rank-1 i na
   syntetycznym dryfie o zadanym S/N) — w v2 raportowany jest tylko rozkład (§5.3).
2. **ρ zakłada dostatecznie stabilne P₃, a stabilności nie mierzę** (§10.3). Tolerancja sięga
   ~±20% wędrówki, ale przy ±40% ρ rośnie do 1.385, a przy ±60% zapada się do 0.510 — czyli silnie
   wędrujący dryfer może zostać **fałszywie zdegradowany**. To warunek stosowalności ρ, nie tylko
   źródło rozrzutu. Potrzebna niezależna miara: szerokość cechy f₃ w LRFS.
3. **Niewyjaśnione ρ > 1 u 25 dobrych dryferów** (J1519-6106: 1.315 ± 0.008 przy rank1 = 0.97).
   Model przewiduje ρ ≤ 1. Wiodąca hipoteza — wędrówka P₃ rzędu ±25–40% — jest ilościowo zgodna
   (§10.3). Inne kandydatki: składowa stojąca odejmująca od kanału parzystego (na syntetyku
   R = 1.32 dla dryfu z rampą), złe uwarunkowanie przy P₃ ≈ 2, ruch nie-sztywny, bi-drifting.
4. **Reżim P₂fit ≤ M/2 liczony względem zadeklarowanego M**, choć `max_dphi` pochodzi już z W₃σ
   (§8.3, wdrożone w v2); spójniej byłoby P₂fit ≤ W₃σ/2.
5. **Listy kandydatów wymagają testu stabilności okna** (§10.4); dotąd żadna go nie ma.
6. **Reżim P₂ > W**: detekcja działa, klasyfikacja przez ρ nie (§9).
7. **Geometria dopasowywana na tych samych danych**, co zawyża `frac_even` i zaniża ρ.
8. **Degeneracja nieusuwalna**: „dryf" i „kontinuum składowych opóźnionych w czasie" to ten sam
   obserwabl.
9. **P₂ raportowane w binach**, w literaturze w stopniach (`360·P₂/nbin`).
10. **Kolejność kanałów niezweryfikowana** (zakładam kanał 0 = dół pasma za nazewnictwem w kodzie).

---

## 12. Historia poprawek

| błąd | objaw | przyczyna |
|---|---|---|
| brak tapera w szablonie | `frac` = 0.767 zamiast 1 | korelacja liniowa sumuje po (N−τ)(M−\|Δ\|) parach wobec NM w zerze |
| projekcja blokowa z członem własnym | 0.48 dla czystego szumu | projekcja na mapę globalną zawiera wkład bloku |
| znak P₂ ze zgadywanki `ridge` | 31 fałszywych kandydatów do degradacji | kryterium musi być na module; znak jest fizyczny |
| swobodny fit P₃ | J2053-7200: 63 zamiast 3.06 | maksimum ucieka w róg dużych P₂/P₃ |
| iloraz zamiast kąta | R = 4.03, NaN-y, potrzeba dwóch progów | tangens ma biegun, sinus nie |
| batch mieszał dane | 85 pełnopasmowych, 430 na 3/16 | różnica czułości 2.25× |
| zmiana sygnatury `_travel_maps` | test reversera czytał 2000 zamiast 3.4·10⁻⁶ | rozpakowanie 4-krotki do 2 zmiennych |
| sprawdzenie okna po fakcie | 7 pulsarów z `ArgumentError` z `Cmd` | `nothing` interpolowane do polecenia |

**Błędne wnioski, wycofane:**

1. Ujemne R przypisane niezgodności katalogowego P₃ — `p3_meas` liczone jest z mapy A, więc przy
   braku ruchu korelacja **wynika** z niskiego R, zamiast je powodować. Przyczyną był znak.
2. Ogłoszenie, że korelacja ρ z P₂fit/M unieważnia klasyfikację — przeceniałem; łańcuch jest
   etykieta → czy istnieje znajdywalna geometria → ρ, czyli metoda działająca.
3. Ograniczenie list kandydatów do reżimu P₂ ≤ M/2 — za ostre.
4. ρ > 1 jako wyjątkowa kategoria na podstawie dwóch obiektów z pilota — w pełnej próbce jest ich
   25, głównie dryferów.
5. **Ogłoszenie, że fit P₂ jest zbędny** na podstawie samego syntetyku — syntetyk był za czysty,
   żeby to rozstrzygnąć (§10.5). Dane realne rozstrzygnęły w drugą stronę.
6. **Postawienie ρ jako produktu głównego** — właściwym produktem dla pytania Song et al. jest
   detekcja `T`/`T_inc`, a ρ jest charakterystyką dodatkową (§4).
7. **Zdiagnozowanie dudnienia jako błędu kalibracji nullu** („PRIORYTET: surogat jest rank-1") —
   to był problem **klasyfikacji**, nie kalibracji. Null jest poprawny, liczba 74σ prawdziwa,
   a błędny był krok „jest uporządkowanie ⇒ jest dryf". Próba naprawy nullu (§7.3) zniszczyła
   detekcję prawdziwego dryfu, co tę diagnozę rozstrzygnęło.

Wspólny mianownik: **sam pomiar nie zmienił się ani razu**; zmieniała się ocena, kiedy wolno go
interpretować. Każda korekta wyszła z testu, nie z rozumowania.

---

## 13. Użycie

```julia
SpaTs.travel_test(vpmout*"J0820-1350"; max_lag=15, show_=false)
SpaTs.travel_test(vpmout*"J0820-1350"; max_lag=15, p2_template=:auto, show_=false)
SpaTs.travel_test(vpmout*"J0601-0527_16"; datafile="pulsar_full_debase.txt", show_=false)
julia --project=. -e 'include("modules/travel.jl"); Travel.selftest()'
```

Argumenty: `max_lag` (2–3·P₃), `max_dphi` (**docelowo z W₃σ**), `hp_halfwin` (50), `nblocks` (4),
`p2_template=:auto` z `p3_template` z LRFS, `p2_cap_frac`, `orth_even`.

Skrypty (`~/claude/work/scripts/`): `travel_batch_full.jl` (533 pulsary, wznawialny, tryb
`--limit N --out PLIK`, `max_dphi = W₃σ ÷ 2` z `onpulse_check.csv`), `travel_v2_summary.py`
(detekcja, `min(cons)` wg siły detekcji, ρ; porównanie v2 z pierwszym przebiegiem), `travel_rho_summary.jl`, `travel_rho_pilot.jl`, `travel_stability.jl`,
`travel_variants.jl`, `check_onpulse.jl`, `travel_check.jl`, `travel_stress.jl`.

Wyniki: `~/output/claude/travel_batch_v3.csv` (aktualny, z `cv_z`/`cv_p`/`adj_z`/`adj_p` jako
listy po nb = 8;16;32;64), `travel_batch_v2.csv` (bez `T_cv`), `travel_batch_full.csv` (pierwszy
przebieg, `max_dphi = M/2`, bez skanu), `onpulse_check.csv`, `travel_stability.csv`,
`travel_rho_all.png`, `travel_rho_pilot.png`.

### Kolejność czytania wyniku

1. **Kontrola off-pulse** — |σ| > 3 dla T lub T_inc unieważnia wszystko poniżej.
2. **Detekcja** — max(T, T_inc) ≥ 5σ. To jest odpowiedź na pytanie A: czy to nie jest czysta
   modulacja amplitudowa.
3. **`min(cons)` ze skanu po długości bloku** — to jest trzecia standardowa wielkość, obok T i ρ.
   Blisko 1 przy każdym podziale: uporządkowanie trwałe, czyli dryf. Przejście przez zero: dudnienie
   dwóch niezależnych zegarów. Blisko zera wszędzie: uporządkowanie nieodtwarzalne, czyli ani dryf,
   ani dudnienie. Oba znaki w `block_proj` to reverser. **Czytać tylko przy silnej detekcji** —
   przy słabej spadek pochodzi z szumu. Bez tego punktu detekcja z punktu 2 nie jest jeszcze
   kandydatem na dryf.
4. **ρ** — tylko gdy P₂fit mieści się w profilu. Blisko 1: ruch ma charakter sztywnej translacji;
   blisko 0: nie. Odniesienie dla szerokiego P₂ w §10.2. Znak to kierunek dryfu (dodatni = od
   wcześniejszych do późniejszych długości, Szary+2022).
5. **Stabilność okna** — dla pojedynczego obiektu obowiązkowa (§10.4).
