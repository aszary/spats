# Test „travel": dryf podpulsów czy modulacja amplitudowa, bez użycia P3

**Stan na 2026-09-23.** Metoda rozstrzygania, czy wzór podpulsów **przemieszcza się** w długości
(dryf, modulacja fazowa), czy tylko **jaśnieje i gaśnie w miejscu** (P3-only, modulacja
amplitudowa) — niezależna od kryterium Song et al. (2023) i wolna od jego głównego obciążenia.

Repozytorium: `github.com/aszary/spats`, gałąź `claude`.
Kod: `modules/travel.jl` (moduł `Travel`), `Plot.travel`, `SpaTs.travel_test`.
Dziennik roboczy: `docs/separations_analysis_log.md` (= `~/claude/work/NOTES.md`, symlink).

---

## 1. Problem

Song et al. (2023) klasyfikują cechę w 2DFS jako dryf, jeśli jej **centroida mocy** jest istotnie
przesunięta względem osi 1/P₂ = 0; brak przesunięcia przy maksimum na niezerowym 1/P₃ to P3-only.
Warunek dodatkowy (Appendix B): znak wyznaczonego P₂ musi być jednoznaczny w granicach błędu.

Słabość jest w estymatorze, nie w idei. Stochastyczna zmienność kształtu pulsu wrzuca moc
skoncentrowaną **wzdłuż osi 1/P₂ = 0**. Gdy ta moc jest asymetryczna — a zwykle jest — centroida
liczona w ręcznie zakreślonym prostokącie przesuwa się od zera i przy niedoszacowanym błędzie
przekracza próg. Pozorny offset bierze się wtedy z **biasu centroidy**, nie z faktycznego
przemieszczania się podpulsów. Dodatkowo test istotności przez shuffle kolejności impulsów bada
hipotezę „to tylko szum", a nie właściwą hipotezę zerową „to mogłaby być czysta modulacja
amplitudowa".

---

## 2. Idea: separowalność i tożsamość zerowa

Niech `δI(n, φ)` będzie fluktuacją natężenia po odjęciu profilu statycznego (`n` — numer impulsu,
`φ` — bin długości). Dwuwymiarowa autokorelacja:

```
K(Δ, τ) = Σ_{n,φ}  δI(n,φ) · δI(n+τ, φ+Δ)
```

**Każda** modulacja amplitudowa jest *separowalna*: wszystkie długości podążają za jedną i tą samą
falą czasową, przeskalowaną rzeczywistym (także ujemnym) współczynnikiem, `δI(n,φ) = a(φ)·w(n)`.
Wtedy

```
K(Δ, τ) = [Σ_φ a(φ)a(φ+Δ)] · [Σ_n w(n)w(n+τ)]
```

a drugi czynnik jest **dokładnie parzysty w τ**:

```
Ĉ_w(−τ) = Σ_n w(n)·w(n−τ)  = [m = n−τ] =  Σ_m w(m+τ)·w(m) = Ĉ_w(τ)
```

To ta sama suma po tych samych parach, tylko przeindeksowana. **Nie ma tu żadnego założenia o
w(n)** — ani stacjonarności, ani okresowości, ani istnienia P3. Stąd statystyka

```
A(Δ, τ) = K(Δ, τ) − K(Δ, −τ)  ≡  0     przy modulacji amplitudowej
```

zeruje się **tożsamościowo**, dla zrealizowanych danych, a nie tylko w wartości oczekiwanej. Fala
bieżąca `δI = f(φ − v·n)` daje `K(Δ,τ) = C_f(Δ − vτ)`, co parzyste w τ nie jest.

Ponieważ `K(Δ,τ) = K(−Δ,−τ)` (przeindeksowanie sum), A jest antysymetryczne w obu argumentach oraz
`A(Δ,0) = A(0,τ) = 0` — cała informacja siedzi w ćwiartce Δ ≥ 1, τ ≥ 1.

### Co dokładnie mierzy A

Bez zakładania separowalności:

```
A(Δ, τ) = Σ_φ [ C_{φ,φ+Δ}(τ) − C_{φ+Δ,φ}(τ) ]
```

czyli A to netto odpowiedź na pytanie **„czy długość φ wyprzedza φ+Δ, czy odwrotnie?"**. Zeruje się
zawsze, gdy żadna długość systematycznie nie wyprzedza innej. To obejmuje więcej niż podręcznikowy
przypadek P3-only:

- dwie składowe **w antyfazie** (a(φ) zmieniające znak) pozostają separowalne → dokładnie zero.
  To ta klasa, która w `PhaseDrift.drift_test_sliding` czyta 25σ fałszywego dryfu;
- dwie składowe z **niezależnymi** modulacjami o różnych P3 → zero w wartości oczekiwanej;
- **zapowane impulsy** są nieszkodliwe z tego samego powodu — wyzerowanie impulsu mnoży każdą
  długość przez ten sam czynnik, więc pole pozostaje separowalne;
- **wysokopasmowy filtr** (ten sam dla każdej długości) też zachowuje separowalność.

---

## 3. Związek z 2DFS — uczciwie

To **nie jest nowy kanał informacji**. Transformata Fouriera z K to dwuwymiarowe widmo mocy, a to,
co przeżywa w A, to asymetria tego widma względem k_φ → −k_φ — formalnie ten sam kanał, którego
używa klasyczne kryterium 2DFS („modulacja amplitudowa daje 2DFS symetryczne względem osi
pionowej" to standardowa wiedza).

Różnica jest w **estymatorze** i jest zasadnicza:

| | Song et al. | test travel |
|---|---|---|
| część symetryczna | wchodzi do centroidy i ją obciąża | rzutowana do zera **algebraicznie**, przed jakąkolwiek estymacją |
| obszar | ręcznie zakreślony prostokąt | cała płaszczyzna (Δ, τ) |
| założenie o P3 | trafienie w bin f₃ | żadne |
| model zerowy | szum (shuffle) | modulacja amplitudowa |

---

## 4. Statystyki detekcji

### 4.1 T — omnibus, koherentna

```
T = Σ_{Δ,τ} A(Δ,τ)²
```

Bez założeń: żadnego tempa dryfu, żadnego P3, żadnego szablonu. Zero w średniej przy modulacji
amplitudowej, dodatnie gdy jakaś długość systematycznie wyprzedza inną.

### 4.2 T_inc — niekoherentna, po blokach

```
T_inc = Σ_b Σ A_b²
```

gdzie `A_b` to mapa policzona w b-tym bloku impulsów. Powód: T liczy się na jednej globalnej mapie,
więc **dryfer, który spędza tyle samo czasu w każdym kierunku, kasuje się i czyta jako brak ruchu**.
W kalibracji syntetycznej zbalansowany reverser daje **T/T_inc = 3.4·10⁻⁶**. J1750-3503 przeżywa
samo T tylko dlatego, że jego epizody są niesymetryczne (28 P w jedną stronę wobec 88 P w drugą).

T_inc płaci wyższym progiem szumu (nb bloków szumu zamiast jednego uśrednienia), więc dla dryfu
stałego jest mniej czuła. Raportować obie; branie lepszej wymaga korekty na trials.

---

## 5. Geometria i iloraz R

### 5.1 Szablon

Dla sztywno dryfującego wzoru `cos(2π(φ/P₂ − n/P₃))`, tożsamość
`cos(x−y) − cos(x+y) = 2 sin x sin y` daje mapę o postaci iloczynu dwóch sinusoid:

```
A(Δ, τ) / K(0,0) = (1 − τ/N)(1 − Δ/M) · 2 · sin(2πΔ/P₂) · sin(2πτ/P₃)
```

Czynnik `(1 − τ/N)(1 − Δ/M)` to **trójkątny taper korelacji liniowej**: przy opóźnieniu (Δ,τ)
sumuje się po (N−τ)(M−|Δ|) parach wobec NM w zerze. To czysta geometria, znana dokładnie.
Pominięcie jej zaniża projekcję o 23% przy N = 600, M = 40 (patrz §9).

### 5.2 Dwie połowy modelu

```
cos(2π(Δ/P₂ − τ/P₃)) = cos(2πΔ/P₂)·cos(2πτ/P₃)  +  sin(2πΔ/P₂)·sin(2πτ/P₃)
                       └──── parzysta w Δ ────┘     └──── nieparzysta w Δ ────┘
```

Dla fali bieżącej obie połowy mają **równe współczynniki**. Modulacja amplitudowa, będąc
separowalną, wkłada wszystko w parzystą i nic w nieparzystą. Definiujemy więc:

- `A = antisym_map(K) = K(Δ,τ) − K(−Δ,τ)` — nieparzysta, mierzy ruch
- `E = sym_map(K)   = K(Δ,τ) + K(−Δ,τ)` — parzysta, mierzy koherentną modulację niezależnie
  od tego, czy wędruje

i projekcje

```
frac_odd  = ⟨A, T_odd⟩  / (K(0,0)·‖T_odd‖²)
frac_even = ⟨E, T_even⟩ / (K(0,0)·‖T_even‖²)
R         = frac_odd / frac_even
```

**R = 1 dla czystego dryfu, 0 dla modulacji amplitudowej.** Znak R niesie kierunek dryfu, więc
klasyfikatorem jest **|R|**.

### 5.3 Dlaczego iloraz, a nie sama projekcja

Ocena `frac_odd` wymagałaby wiedzy, ile zwraca prawdziwy dryfer — a oczywista próbka referencyjna
(pulsary już oznaczone jako `drift`) to dokładnie zbiór podejrzany o skażenie. Kalibracja na nim
jest **cyrkularna**: ucząc się na mieszance dowiadujemy się, że „dryfer może mieć frac ≈ 0", próg
spada i nic nie da się zdegradować.

R bierze miarę **z tego samego pulsara**. Siła modulacji, zawartość harmonicznych, zanik koherencji
z opóźnieniem i taper mnożą obie połowy identycznie i **kasują się w ilorazie**. Skala problemu:
wśród dryferów `frac_odd` rozciąga się na 25× (0.0047 do 0.118), a R tylko na 2× (0.63 do 1.32).
Na progu bezwzględnym J1110-5637 i J2053-7200 wyglądałyby jak „prawie zero" obok J0820-1350 i
zostałyby błędnie zdegradowane.

### 5.4 Wyznaczanie geometrii

`|P₂|` dopasowywane skanem (siatka logarytmiczna, 48 punktów) **maksymalizującym projekcję
parzystą**. To metodologicznie istotne: parzysta mierzy koherentną modulację niezależnie od ruchu,
więc wybór geometrii **nie może wyprodukować ruchu**. Znak P₂ wychodzi potem ze znaku projekcji
nieparzystej — to odpowiedź fizyczna, nie wybór dopasowania.

`P₃` **nie jest dopasowywane** — brane jest z pomiaru LRFS w `params.json`. Swobodny fit ucieka w
róg dużych P₂/P₃, gdzie szablon jest prawie stały i dopasowuje się do gładkiego tła mapy (patrz §9).

Skan liczony w postaci zamkniętej, dwoma mnożeniami macierzy: z `W[t,d] = E[t,d]·(1−t/N)(1−d/M)`
licznik faktoryzuje się do iloczynów macierzowych, a `‖szablon‖²` rozpada się dokładnie na czynnik
τ razy czynnik Δ.

---

## 6. Model zerowy i kontrole

### 6.1 Surogaty

Pod H₀ sygnał wnosi do A **dokładnie zero** — nie w przybliżeniu — więc cała wariancja pochodzi z
członów szumowych. To radykalnie upraszcza null:

- **sygnał**: wiodący mod SVD danych (czyli dokładnie separowalny, czyli dokładnie H₀),
  przeskalowany tak, by całkowita moc równała się odszumionej `Σ X² − N·M·σ²`. Kształt nie ma
  znaczenia, liczy się amplituda, bo sygnał wchodzi tylko przez człon skrośny sygnał×szum;
- **szum**: bootstrapowany z **własnego obszaru off-pulse pulsara** — ciągły pasek tej samej
  szerokości co okno on-pulse, z losowym przesunięciem cyklicznym w czasie. Zachowuje realny poziom
  szumu, jego korelację czasową *oraz* korelację między sąsiednimi binami długości. Gdy profil jest
  szerszy niż najdłuższy ciągły fragment off-pulse, pasek brany jest cyklicznie po liście binów
  off-pulse (flaga `offpulse_wrapped`).

### 6.2 Kontrola off-pulse

Ta sama statystyka liczona na pasku off-pulse, gdzie sygnału nie ma — **musi wyjść zgodna z zerem**.
Jedyna realna podatność metody to szum niesymetryczny w czasie (dryf wzmocnienia, narastające RFI,
źle odjęta linia bazowa), i właśnie to ta kontrola wyłapuje. W przebiegu po 515 pulsarach odrzuciła
**1**.

### 6.3 Spójność blokowa

Projekcja mapy każdego bloku impulsów na sumę pozostałych (**leave-one-out** — projekcja na mapę
globalną zawierałaby człon własny i dawała ~1/√n_bloków dla samego szumu). Zero dla szumu, jeden dla
ruchu powtarzalnego przez całą obserwację, oba znaki dla dryfera odwracającego kierunek.

---

## 7. Dane

Analiza idzie na **pełnym paśmie**. Katalogi `<PSR>_16/` zawierały tylko podpasma, tworzone z
16-kanałowego `pulsar.spCf16` przez `paz -Z` (zapuje wymienione kanały, zostaje reszta):

| plik | `paz -Z` | zostaje | szerokość |
|---|---|---|---|
| `low` | `3-15` | kanały 0–2 | 3/16 |
| `mid` | `0-6 9-15` | kanały 7–8 | 2/16 |
| `high` | `0-12` | kanały 13–15 | 3/16 |

Użycie podpasma kosztuje czynnik ~2.25 w RMS (J0601-0527: 0.0275 pełne pasmo wobec 0.062 dla `low`,
zgodnie z √(16/3)). Pełne pasmo odtwarzane jest raz i zostaje w katalogu `_16`:

```
pulsar.full             pam -F -u <dir> -e full pulsar.spCf16
pulsar_full.debase.gg   pmod -onpulse "<bin_st> <bin_end>" -device /NULL -debase
pulsar_full_debase.txt  pdv -t -F
```

Katalogi bez `_16` mają `pulsar.debase.txt` już pełnopasmowe (`process_psrdata` nie robi `paz -Z`).

**Wstępne przetwarzanie**: wysokopasmowy filtr biegnącą średnią, `hp_halfwin = 50` impulsów —
świadomie łagodny, usuwa wędrowanie wzmocnienia bez tykania modulacji. Ten sam filtr dla każdej
długości, więc separowalność (a z nią tożsamość zerowa) zostaje nienaruszona.

---

## 8. Walidacja

### 8.1 Syntetyczna (`Travel.selftest()`)

| test | wynik |
|---|---|
| FFT z paddingiem zerowym vs suma wprost | zgodność 4.3·10⁻¹⁶ |
| pole separowalne: wędrujące P3 + nulling + skok znaku a(φ) | max\|A\| = 1.1·10⁻¹⁶ × K(0,0) |
| syntetyczna fala bieżąca | P₂ = 18.0 (prawda 18), P₃ = 12.0 (prawda 12), rank1 = 1.00 |
| projekcja matched w prawdziwej geometrii | frac = 1.000 |
| R dla czystego dryfu | 1.000 |
| R dla fali bieżącej + stojącej o tych samych okresach | 0.607 przy przewidzianym analitycznie 0.600 |
| zbalansowany reverser | T/T_inc = 3.4·10⁻⁶ |

Przewidywanie dla mieszaniny wyprowadza się z tego, że fala stojąca to pół bieżącej w przód plus pół
w tył: `R = (a_f² − a_b²)/(a_f² + a_b²)`.

### 8.2 Na danych rzeczywistych (pełne pasmo)

| pulsar | charakter | T | \|R\| |
|---|---|---|---|
| J0820-1350 (B0818−13) | podręcznikowy dryfer | >999σ | **1.08** (R < 0, dryf ujemny) |
| J2053-7200 | wobble P3 opróżnia bin LRFS | 110σ | **1.32** |
| J1750-3503 | dryf odwracający kierunek | 163σ | **0.63** |
| J1110-5637 | silny koherentny dryf | 81σ | **0.75** |
| J1907+0731 | sklasyfikowany P3-only | 1.7σ | **0.31 ± 0.11** |

J0820-1350 odtwarza P₃ = 4.8 wobec katalogowych 4.78, nie dostawszy tej wartości. J2053-7200 to
przypadek, w którym `PhaseDrift.phase_modulation` daje 0.5σ, bo wobble P3 opróżnia globalny bin FFT.

### 8.3 Populacyjna (przebieg podpasmowy, 515 pulsarów)

Etykiety Song et al. wchodzą wyłącznie jako **zbiór testowy** — nic się na nich nie uczy.

| etykieta | z detekcją ≥5σ | mediana R | mediana R przy zgodnym P₃ |
|---|---|---|---|
| drift | 293 | 0.89 | **1.03** (n=204) |
| p3only | 49 | 0.12 | **0.16** (n=14) |

Te liczby pochodzą z przebiegu na **niejednorodnych** danych (85 pulsarów pełnopasmowych, 430 na
3/16 pasma) i przed poprawką znaku — traktować jako wstępne. Kierunek separacji jest jednak na tyle
wyraźny, że raczej się nie odwróci.

---

## 9. Ograniczenia i problemy otwarte

1. **R nie jest ścisłym ułamkiem i może przekroczyć 1** (zaobserwowane 1.08, 1.32). Mianownik
   zbiera też parzystą projekcję modulacji nie-wędrującej, a ta nie ma ustalonego znaku. Odczyt jest
   **porządkowy, nie dosłowny**. Patologia siedzi przy górnym końcu, a degradacja rozgrywa się przy
   dolnym, gdzie do zera dąży licznik i iloraz zachowuje się dobrze.

2. **Reverser zaniża R** (J1750-3503: 0.63, najniżej z czwórki) — globalna mapa częściowo się kasuje.
   Przed odczytem R sprawdzić `T_inc` i `block_proj`; przy obu znakach ograniczyć się do jednego
   epizodu przez `pulse_st`/`pulse_end`.

3. **Null jest rank-1, a pole może być rank ≥ 2.** J1907+0731 daje `T_inc` = 4.8σ przy `T` = 1.7σ,
   przy czystej kontroli off-pulse (T_off = −0.3σ, T_inc_off = 0.2σ) — więc to nie asymetria szumu.
   Projekcje blokowe ≈ 0.00, czyli mapy bloków **nie zgadzają się ze sobą**: sygnatura losowego
   uporządkowania między niezależnymi modami, a nie dryfu. Jeśli pole jest rank ≥ 2 z niezależnymi
   modami czasowymi, surogat zaniża wariancję i zawyża istotność. Możliwa poprawka: surogat rank-r
   z niezależnie przesuwanymi w czasie modami. **Do czasu rozstrzygnięcia: `T_inc` bez zgodności
   blokowej nie jest kandydatem na dryf.**

4. **Geometria dopasowywana na tych samych danych**, na których mierzona jest projekcja, co zawyża
   `frac_even` i tym samym lekko **zaniża R**. Surogaty używają ustalonego szablonu i tego wyboru
   nie odtwarzają.

5. **Degeneracja nieusuwalna**: „dryf" i „kontinuum składowych opóźnionych w czasie" to *ten sam
   obserwabl*. Test tego nie łamie i nie powinien.

6. **Podłoga czułości**: dryf tak wolny, że wychylenie fazy przez cały profil tonie w szumie, jest
   niewykrywalny. Kwantyfikowalne jako limit.

7. **Kolejność kanałów niezweryfikowana** — zakładam kanał 0 = dół pasma za nazewnictwem w kodzie,
   ale `vap` i `psredit` nie wystawiają częstotliwości per kanał dla tych plików. Bez znaczenia przy
   pełnym paśmie, istotne przy ewentualnym porównaniu międzyczęstotliwościowym.

---

## 10. Historia poprawek

Błędy wyłapane w trakcie — zapisane, bo każdy dawał wynik wyglądający wiarygodnie.

| błąd | objaw | przyczyna |
|---|---|---|
| brak tapera w szablonie | `frac` = 0.767 zamiast 1 dla syntetyku o znanej geometrii | korelacja liniowa sumuje po (N−τ)(M−\|Δ\|) parach wobec NM w zerze |
| projekcja blokowa z członem własnym | 0.48 dla czystego szumu | projekcja na mapę globalną zawiera własny wkład bloku → ~1/√n |
| znak P₂ ze zgadywanki `ridge` | ujemne R czytane jako „brak ruchu", 31 fałszywych kandydatów do degradacji | kryterium musi być na \|R\|; znak jest fizyczny |
| swobodny fit P₃ | J2053-7200 dopasowało 63 zamiast 3.06, R spadło do 0.12 | maksimum ucieka w róg dużych P₂/P₃, gdzie szablon jest prawie stały i łapie gładkie tło; ortogonalizacja względem stałej nie wystarczyła |
| batch mieszał dane | 85 pulsarów pełnopasmowych, 430 na 3/16 pasma | różnica czułości ~2.25× w amplitudzie, niedopuszczalna dla frakcji detekcji |
| zmiana sygnatury `_travel_maps` | test reversera czytał 2000 zamiast 3.4·10⁻⁶ | rozpakowanie 4-krotki do 2 zmiennych; złapane przez selftest |

Osobno — **błędna diagnoza, którą trzeba odnotować**: początkowo przypisałem ujemne R niezgodności
katalogowego P₃, powołując się na medianę `p3_meas/p3_cat` = 0.37 w grupie R < 0.3 wobec 1.02 w
grupie R > 0.7. To wnioskowanie było nieuprawnione: `p3_meas` liczone jest z mapy A, więc gdy ruchu
nie ma, A jest szumem i `p3_meas` traci sens — korelacja równie dobrze **wynika** z niskiego R,
zamiast je powodować. Właściwą przyczyną był znak.

---

## 11. Użycie

```julia
# pojedynczy pulsar
SpaTs.travel_test(vpmout*"J0820-1350"; max_lag=15, show_=false)

# z projekcją na deklarowaną geometrię (kierunek degradacji)
SpaTs.travel_test(vpmout*"J0820-1350"; max_lag=15, p2_template=:auto, show_=false)

# podpasmo z katalogu _16
SpaTs.travel_test(vpmout*"J0601-0527_16"; datafile="pulsar_full_debase.txt", show_=false)

# weryfikacja fundamentów metody
julia --project=. -e 'include("modules/travel.jl"); Travel.selftest()'
```

Kluczowe argumenty: `max_lag` (domyślnie 40; przy znanym P₃ dobre jest 2–3·P₃), `max_dphi`
(domyślnie połowa szerokości on-pulse), `hp_halfwin` (50), `nblocks` (4 — blok ma być krótszy niż
oczekiwany epizod dryfu), `p2_template=:auto` z `p3_template` z LRFS.

Skrypty wsadowe: `~/claude/work/scripts/travel_batch_full.jl` (533 pulsary, wznawialny),
`travel_summary.jl` (rozkład R i listy kandydatów), `travel_check.jl` / `travel_stress.jl` /
`travel_ratio.jl` (walidacja).

### Jak czytać wynik

1. **Kontrola off-pulse** — jeśli \|σ\| > 3 dla T lub T_inc, szum nie jest symetryczny w czasie i
   reszta nie znaczy nic.
2. **Detekcja** — max(T, T_inc) ≥ 5σ. Bez niej R mierzy czułość, nie fizykę pulsara.
3. **Spójność blokowa** — oba znaki oznaczają reverser; wartości ≈ 0 przy istotnym T_inc oznaczają
   niezależne mody, nie dryf (punkt 3 w §9).
4. **\|R\|** — blisko 1: koherentna modulacja wędruje; blisko 0: nie wędruje. Znak R to kierunek
   dryfu w konwencji Szary+2022 (dodatni = od wcześniejszych do późniejszych długości).
