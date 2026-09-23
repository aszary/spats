# Test „travel": dryf podpulsów czy modulacja amplitudowa, bez użycia P3

**Stan na 2026-09-23.** Metoda rozstrzygania, czy wzór podpulsów **przemieszcza się** w długości
(dryf, modulacja fazowa), czy tylko **jaśnieje i gaśnie w miejscu** (P3-only, modulacja
amplitudowa) — niezależna od kryterium Song et al. (2023) i wolna od jego głównego obciążenia.

Repozytorium: `github.com/aszary/spats`, gałąź `claude`.
Kod: `modules/travel.jl` (moduł `Travel`), `Plot.travel`, `SpaTs.travel_test`.
Skrypty: `~/claude/work/scripts/travel_*.jl`, `check_onpulse.jl`.
Dziennik roboczy: `docs/separations_analysis_log.md` (= `~/claude/work/NOTES.md`, symlink).

---

## 1. Problem

Song et al. (2023) klasyfikują cechę w 2DFS jako dryf, jeśli jej **centroida mocy** jest istotnie
przesunięta względem osi 1/P₂ = 0; brak przesunięcia przy maksimum na niezerowym 1/P₃ to P3-only.

Słabość jest w estymatorze, nie w idei. Stochastyczna zmienność kształtu pulsu wrzuca moc
skoncentrowaną **wzdłuż osi 1/P₂ = 0**. Gdy ta moc jest asymetryczna — a zwykle jest — centroida
liczona w ręcznie zakreślonym prostokącie przesuwa się od zera i przy niedoszacowanym błędzie
przekracza próg. Pozorny offset bierze się wtedy z **biasu centroidy**, nie z faktycznego
przemieszczania się podpulsów. Dodatkowo test istotności przez shuffle kolejności impulsów bada
hipotezę „to tylko szum", a nie właściwą hipotezę zerową „to mogłaby być czysta modulacja
amplitudowa".

---

## 2. Idea: separowalność i tożsamość zerowa

Niech `δI(n, φ)` będzie fluktuacją natężenia po odjęciu profilu statycznego. Dwuwymiarowa
autokorelacja:

```
K(Δ, τ) = Σ_{n,φ}  δI(n,φ) · δI(n+τ, φ+Δ)
```

**Każda** modulacja amplitudowa jest *separowalna*: wszystkie długości podążają za jedną falą
czasową, przeskalowaną rzeczywistym (także ujemnym) współczynnikiem, `δI(n,φ) = a(φ)·w(n)`. Wtedy
`K(Δ,τ) = [Σ_φ a(φ)a(φ+Δ)]·[Σ_n w(n)w(n+τ)]`, a drugi czynnik jest **dokładnie parzysty w τ**:

```
Ĉ_w(−τ) = Σ_n w(n)·w(n−τ)  = [m = n−τ] =  Σ_m w(m+τ)·w(m) = Ĉ_w(τ)
```

To ta sama suma po tych samych parach, tylko przeindeksowana — **żadnego założenia o w(n)**: może
wędrować okresem, gasnąć, przestać być okresowe. Stąd

```
A(Δ, τ) = K(Δ, τ) − K(Δ, −τ)  ≡  0     przy modulacji amplitudowej
```

**tożsamościowo**, dla zrealizowanych danych, a nie w wartości oczekiwanej. Fala bieżąca
`δI = f(φ − v·n)` daje `K = C_f(Δ − vτ)`, co parzyste w τ nie jest.

Bez zakładania separowalności `A(Δ,τ) = Σ_φ [C_{φ,φ+Δ}(τ) − C_{φ+Δ,φ}(τ)]`, czyli netto odpowiedź
na pytanie **„czy długość φ wyprzedza φ+Δ?"**. Zeruje się dla: składowych w antyfazie (a(φ)
zmieniające znak pozostaje separowalne), niezależnych modulacji o różnych P3, zapowanych impulsów
i po filtrze wysokopasmowym (ten sam dla każdej długości).

---

## 3. Związek z 2DFS

To **nie jest nowy kanał informacji** — transformata K to widmo mocy 2D, a to, co przeżywa w A, to
asymetria względem k_φ → −k_φ, formalnie ten sam kanał co klasyczne kryterium 2DFS. Różnica jest w
estymatorze:

| | Song et al. | test travel |
|---|---|---|
| część symetryczna | wchodzi do centroidy i ją obciąża | rzutowana do zera **algebraicznie** |
| obszar | ręcznie zakreślony prostokąt | cała płaszczyzna (Δ, τ) |
| założenie o P3 | trafienie w bin f₃ | żadne (P3 wchodzi dopiero w szablonie, §5) |
| model zerowy | szum (shuffle) | modulacja amplitudowa |

---

## 4. Statystyki detekcji

**`T = Σ A(Δ,τ)²`** — omnibus, bez żadnych założeń o geometrii.

**`T_inc = Σ_b Σ A_b²`** — suma niekoherentna po blokach impulsów. T liczy się na jednej globalnej
mapie, więc dryfer spędzający tyle samo czasu w każdym kierunku kasuje się. Kalibracja: zbalansowany
reverser daje **T/T_inc = 3.4·10⁻⁶**. T_inc płaci wyższym progiem szumu, więc dla dryfu stałego jest
mniej czuła — raportować obie.

---

## 5. Klasyfikator ρ

### 5.1 Szablon i dwie połowy modelu

Dla wzoru `cos(2π(φ/P₂ − n/P₃))`:

```
A(Δ,τ)/K(0,0) = (1 − τ/N)(1 − Δ/M) · 2 · sin(2πΔ/P₂) · sin(2πτ/P₃)
```

Czynnik `(1 − τ/N)(1 − Δ/M)` to **trójkątny taper korelacji liniowej** — przy opóźnieniu (Δ,τ)
sumuje się po (N−τ)(M−|Δ|) parach wobec NM w zerze. Czysta geometria; pominięcie zaniża projekcję
o 23% przy N = 600, M = 40.

Rozwinięcie cosinusa różnicy daje dwie połowy **o równych współczynnikach**:

```
cos(2π(Δ/P₂ − τ/P₃)) = cos(2πΔ/P₂)cos(2πτ/P₃)  +  sin(2πΔ/P₂)sin(2πτ/P₃)
                       └── parzysta w Δ ──┘        └── nieparzysta w Δ ──┘
```

Nieparzystą mierzy `A = antisym_map(K)`, parzystą `E = sym_map(K)`. Modulacja amplitudowa, będąc
separowalną, wkłada wszystko w parzystą i nic w nieparzystą.

### 5.2 Dlaczego ρ, a nie iloraz R

Pierwotnie używałem `R = frac_odd/frac_even` (1 dla dryfu, 0 dla AM). To **tangens** kąta w
płaszczyźnie (odd, even) i ma biegun: gdy projekcja parzysta maleje, R eksploduje (zaobserwowane
4.03 i 3.40) i wymaga strażnika „mianownik istotny", który generuje NaN-y. Lekarstwem byłyby progi
na `rank1` i spójność blokową — czyli dwa parametry swobodne.

Zamiast tego **sinus tego samego kąta**:

```
ρ = √2 · frac_odd / √(frac_odd² + frac_even²)
```

ρ = 1 dla sztywnego dryfu, 0 dla modulacji amplitudowej, **ograniczone przez √2**, zawsze
zdefiniowane, zero parametrów. Znak niesie kierunek dryfu, więc klasyfikatorem jest **|ρ|**.
`R` zostaje w wyniku, bo relacje w §5.1 są w nim sformułowane.

### 5.3 Dlaczego iloraz w ogóle — problem cyrkularności

Ocena samego `frac_odd` wymagałaby wiedzy, ile zwraca prawdziwy dryfer, a oczywista próbka
referencyjna (pulsary oznaczone `drift`) to dokładnie zbiór podejrzany o skażenie. ρ bierze miarę
**z tego samego pulsara**: siła modulacji, harmoniczne, zanik koherencji i taper mnożą obie połowy
identycznie i kasują się. Skala problemu: wśród dryferów `frac_odd` rozciąga się na 25×, a ρ na 2×.

### 5.4 Wyznaczanie geometrii

**|P₂| dopasowywane** skanem (siatka logarytmiczna, 48 punktów) maksymalizującym **projekcję
parzystą**. To istotne: parzysta mierzy koherentną modulację niezależnie od ruchu, więc wybór
geometrii **nie może wyprodukować ruchu**. Znak P₂ wychodzi ze znaku projekcji nieparzystej.

**P₃ NIE jest dopasowywane** — brane z pomiaru LRFS w `params.json`. Swobodny fit ucieka w róg
dużych P₂/P₃, gdzie szablon jest prawie stały i łapie gładkie tło mapy (J2053-7200 dopasowało 63
zamiast 3.06, ρ spadło do 0.12). Ortogonalizacja szablonu względem modelu stałego tego **nie
naprawia** — sprawdzone, tło nie leży wzdłuż kierunku stałej.

Skan liczony w postaci zamkniętej (dwa mnożenia macierzy).

---

## 6. Model zerowy i kontrole

**Surogaty.** Pod H₀ sygnał wnosi do A dokładnie zero, więc cała wariancja pochodzi z członów
szumowych. Surogat = wiodący mod SVD danych (czyli dokładnie separowalny) przeskalowany do
odszumionej mocy, plus szum **bootstrapowany z własnego off-pulse'u pulsara** (ciągły pasek tej
samej szerokości, losowe przesunięcie cykliczne w czasie). Gdy profil jest szerszy niż najdłuższy
ciągły fragment off-pulse, pasek brany cyklicznie po liście binów (flaga `offpulse_wrapped`).

**Kontrola off-pulse.** Ta sama statystyka na pasku bez sygnału — musi wyjść zgodna z zerem. Jedyna
realna podatność metody to szum niesymetryczny w czasie (dryf wzmocnienia, RFI, zła linia bazowa).
W pełnym przebiegu odrzuciła **3 z 515**.

**Spójność blokowa.** Projekcja mapy każdego bloku na sumę pozostałych (**leave-one-out** —
projekcja na mapę globalną zawierałaby człon własny i dawała ~1/√n_bloków dla szumu).

---

## 7. Dane: pełne pasmo

Katalogi `<PSR>_16/` zawierały tylko podpasma, tworzone z 16-kanałowego `pulsar.spCf16` przez
`paz -Z` (zapuje wymienione kanały, zostaje reszta):

| plik | `paz -Z` | zostaje | szerokość |
|---|---|---|---|
| `low` | `3-15` | kanały 0–2 | 3/16 |
| `mid` | `0-6 9-15` | kanały 7–8 | 2/16 |
| `high` | `0-12` | kanały 13–15 | 3/16 |

Podpasmo kosztuje czynnik ~2.25 w RMS (J0601-0527: 0.0275 pełne pasmo wobec 0.062 dla `low`,
zgodnie z √(16/3)). **Pierwszy przebieg mieszał dane** — 85 pulsarów pełnopasmowych i 430 na 3/16
pasma — co dla frakcji detekcji jest niedopuszczalne. Pełne pasmo odtwarzane jest raz i **zostaje w
katalogu `_16`**:

```
pulsar.full             pam -F -u <dir> -e full pulsar.spCf16      8.5 MB
pulsar_full.debase.gg   pmod -onpulse "<bst> <ben>" -device /NULL -debase   8.5 MB
pulsar_full_debase.txt  pdv -t -F                                   55 MB
```

~7 s na pulsara, ~38 GB dla 534 katalogów. Katalogi bez `_16` mają `pulsar.debase.txt` już
pełnopasmowe.

**Wstępne przetwarzanie**: filtr wysokopasmowy biegnącą średnią, `hp_halfwin = 50` impulsów —
świadomie łagodny. Ten sam dla każdej długości, więc separowalność (a z nią tożsamość zerowa)
zostaje nienaruszona.

---

## 8. Walidacja

### 8.1 Syntetyczna (`Travel.selftest()`)

| test | wynik |
|---|---|
| FFT z paddingiem zerowym vs suma wprost | 4.3·10⁻¹⁶ |
| pole separowalne: wędrujące P3 + nulling + skok znaku a(φ) | max\|A\| = 1.1·10⁻¹⁶ × K(0,0) |
| syntetyczna fala bieżąca | P₂ = 18.0 (prawda 18), P₃ = 12.0 (prawda 12), rank1 = 1.00 |
| projekcja matched w prawdziwej geometrii | frac = 1.000 |
| ρ dla czystego dryfu | 1.000 |
| ρ dla fali bieżącej + stojącej o tych samych okresach | 0.607 przy przewidzianym analitycznie 0.600 |
| zbalansowany reverser | T/T_inc = 3.4·10⁻⁶ |

Przewidywanie dla mieszaniny: fala stojąca to pół bieżącej w przód plus pół w tył, więc
`ρ = (a_f² − a_b²)/(a_f² + a_b²)` w wersji ilorazowej.

### 8.2 Krzywa odniesienia ρ(P₂/M)

Syntetyk z obwiednią gaussowską i szumem, sztywny dryf o rosnącym P₂:

| P₂/M | 0.25 | 0.5 | 1.0 | 1.5 | 2.0 | 3.0 | 4.0 |
|---|---|---|---|---|---|---|---|
| **ρ** | 1.00 | 1.00 | 0.98 | 0.87 | 0.85 | 0.69 | 0.64 |

Dryf o P₂ szerszym niż profil jest tłumiony **łagodnie, nie zapada się**. To jest materiał
odniesienia, nie kryterium stosowane przez kod.

### 8.3 Pełna próbka (533 pulsary, pełne pasmo)

533 policzone, 18 błędów (5 bez danych, 7 bez okna on-pulse, 6 z profilem szerszym niż off-pulse),
kontrola off-pulse odrzuciła 3 z 515, do analizy 512, z detekcją ≥5σ **387**.

Etykiety Song et al. użyte **wyłącznie jako zbiór testowy** — nic się na nich nie uczy.

W **reżimie rozdzielczym** (P₂fit ≤ M/2, czyli szablon zamyka pełny cykl w przeszukiwanym zakresie
Δ — warunek rozdzielczości, nie dobrany próg), 179 pulsarów:

| | n | mediana ρ | kwartyle |
|---|---|---|---|
| **drift** | 172 | **1.011** | 0.830–1.087 |
| p3only | 7 | 0.324 | 0.092–0.926 |

**Mediana 1.011 na 172 niezależnych pulsarach przy przewidywaniu teorii dokładnie 1 i zerowej
liczbie parametrów swobodnych** — to najmocniejszy wynik tej pracy.

Poza reżimem: drift 0.32 (n≈121), p3only 0.15 (n≈45).

Rozkład samej geometrii też rozróżnia: mediana P₂fit/M wynosi **0.39 dla drift** (63% poniżej 1)
wobec **1.39 dla p3only** (21% poniżej 1). Dryfer ma znajdywalne P₂ wewnątrz profilu, P3-only nie ma
preferowanej geometrii i fit ucieka w górę siatki.

### 8.4 ρ nie zależy od S/N

Pozorna zależność w próbce zbiorczej okazała się **efektem składu (paradoks Simpsona)**:

| | T ∈ [5,30) | [30,300) | [300,∞) |
|---|---|---|---|
| ρ przy P₂/M < 0.5 (drift) | 0.912 | 1.006 | 1.032 |
| ρ przy P₂/M ≥ 1 (drift) | 0.342 | 0.326 | 0.318 |
| udział P₂/M < 0.5 wśród dryferów | 40% | 54% | 61% |

W obu grupach ρ jest płaskie; S/N steruje tylko tym, do której grupy pulsar trafia.

---

## 9. Okna on-pulse i stabilność wyniku

### 9.1 Same okna są w porządku

Na 515 pulsarach (profile średnie z `pdv -t -F -T`, `check_onpulse.jl`):

| | mediana | kwartyle |
|---|---|---|
| M (okno z `params.json`) | 116 binów | 90–157 |
| W₃σ (zasięg emisji) | 56 binów | 37–80 |
| **W₃σ / M** | **0.49** | 0.38–0.60 |

Okna są medianowo dwukrotnie szersze niż kontur 3σ, ale szczyt leży wewnątrz okna u **512 z 515**.
To normalna praktyka — próg 3σ obcina skrzydła — a nie błąd automatu `pmod`.

### 9.2 Poszerzanie okna szkodzi, nie pomaga

Syntetyk, ta sama emisja, coraz szersze okno (P₂ prawdziwe = 150 binów):

| M | frac_odd | ρ | P₂fit |
|---|---|---|---|
| 100 (dopasowane) | 0.384 | **0.872** | **150.2** ✓ |
| 150 | 0.213 | 0.802 | 185.7 |
| 200 | 0.123 | 0.750 | 231.5 |
| 250 | 0.070 | 0.653 | **308.9** ✗ |
| 300 | — | — | błąd: on-pulse przerósł off-pulse |

Zawyżone M rozcieńcza sygnał (K(0,0) rośnie o kolumny szumu), psuje ρ, **zawyża P₂** (taper zakłada
dane sięgające M) i zjada obszar off-pulse. Zwiększanie samego zasięgu Δ przy stałym oknie też nie
pomaga (md = 0.5M vs 0.9M: przy P₂/M = 3 ρ spada z 0.688 na 0.614). Powód jest zasadniczy:
**informacja o okresie w długości pochodzi wyłącznie z obszaru, który świeci** — dokładanie szumu
nie wydłuża bazy pomiarowej.

### 9.3 Test stabilności okna — kluczowe zastrzeżenie

ρ policzone przy oknie 1.0 / 1.5 / 2.0 × W₃σ:

| pulsar | rola | ρ @1.0× | ρ @1.5× | ρ @2.0× |
|---|---|---|---|---|
| J0034-0721 | kontrola (dryfer) | 1.042 | 1.067 | 1.048 |
| J0151-0635 | kontrola (dryfer) | 0.970 | 1.011 | 1.006 |
| J1239+2453 | kand. degradacji | **0.011** | 0.374 | **0.461** |
| J2048-1616 | kand. degradacji | 0.021 | 0.003 | **0.481** |
| J1921+2003 | kand. degradacji | 0.055 | 0.084 | 0.035 |
| J1810-5338 | kand. promocji | 0.566 | 0.438 | 0.770 |

**Prawdziwe dryfery są odporne** (ρ i P₂ powtarzalne do kilku procent) — więc walidacja z §8.3 nie
zależy od wyboru okien. **Ale obiekty o niskim ρ są chwiejne**, nawet 40-krotnie. Mechanizm jest
zrozumiały: ρ = odd/hypot(odd, even), więc przy liczniku bliskim zeru drobne zmiany zawartości okna
nim rzucają.

**Konsekwencja: żadna lista kandydatów nie ma prawa iść dalej bez testu stabilności okna.** To nie
jest dostrajany parametr, tylko wymóg odporności — obiekt, który przy 1.0× daje 0.01, a przy 2.0×
daje 0.48, po prostu nie jest zmierzony. Dotyczy to obu kierunków; w szczególności J1810-5338,
najczystszy kandydat do promocji z pierwszego przejścia, **nie przeszedł** tego testu.

---

## 10. Ograniczenia i problemy otwarte

1. **Listy kandydatów nie są jeszcze wiarygodne** — patrz §9.3. Potrzebny test stabilności okna na
   każdym obiekcie (trzy przebiegi zamiast jednego, albo tylko dla obiektów skrajnych).

2. **Systematyczne ρ > 1 u dobrych dryferów.** 25 obiektów przekracza 1.1 na 3σ, w tym J1519-6106 z
   ρ = 1.315 ± 0.008 przy rank1 = 0.97 i spójności 0.91. Model przewiduje ρ ≤ 1. Mediana trzyma się
   (1.011), ale rozrzut w górę nie jest szumem i **nie jest wyjaśniony**. Do rozstrzygnięcia przed
   publikacją. Kandydaci: składowa stojąca odejmująca od kanału parzystego (zademonstrowane na
   syntetyku, R = 1.32 dla dryfu z rampą), złe uwarunkowanie przy P₃ ≈ 2 (przy P₃ = 2.05 norma
   szablonu nieparzystego jest 6× mniejsza niż parzystego), ruch nie-sztywny, bi-drifting.

3. **Gdzie metoda ma moc.** Poza reżimem rozdzielczym obie klasy zapadają się do niskiego ρ. Ale
   krzywa z §8.2 mówi, że prawdziwy dryf przy P₂/M = 1.5–3 dałby 0.69–0.87, a realne „drift" z
   P₂/M ≥ 1 dają 0.32 — więc metoda tam **jednak rozróżnia**, słabiej. Ograniczenie do P₂ ≤ M/2
   było za ostre; właściwym odniesieniem jest krzywa, nie cięcie próbki.

4. **Null jest rank-1, a pole może być rank ≥ 2.** J1907+0731 daje T_inc = 4.8σ przy T = 1.7σ i
   czystej kontroli off-pulse, przy projekcjach blokowych ≈ 0 — sygnatura losowego uporządkowania
   między niezależnymi modami, nie dryfu. Reguła robocza: **T_inc bez zgodności blokowej nie jest
   kandydatem na dryf.**

5. **Geometria dopasowywana na tych samych danych**, na których mierzona jest projekcja, zawyża
   `frac_even` i tym samym lekko zaniża ρ. Surogaty używają ustalonego szablonu.

6. **Degeneracja nieusuwalna**: „dryf" i „kontinuum składowych opóźnionych w czasie" to ten sam
   obserwabl.

7. **P₂ raportowane w binach**, a w literaturze w stopniach (`360·P₂/nbin`). Do zamiany w tabeli
   wynikowej.

8. **Kolejność kanałów niezweryfikowana** — zakładam kanał 0 = dół pasma za nazewnictwem w kodzie,
   ale `vap`/`psredit` nie wystawiają częstotliwości per kanał dla tych plików. Bez znaczenia przy
   pełnym paśmie, istotne przy porównaniu międzyczęstotliwościowym.

---

## 11. Historia poprawek

Każdy z tych błędów dawał wynik wyglądający wiarygodnie.

| błąd | objaw | przyczyna |
|---|---|---|
| brak tapera w szablonie | `frac` = 0.767 zamiast 1 dla syntetyku o znanej geometrii | korelacja liniowa sumuje po (N−τ)(M−\|Δ\|) parach wobec NM w zerze |
| projekcja blokowa z członem własnym | 0.48 dla czystego szumu | projekcja na mapę globalną zawiera wkład bloku → ~1/√n |
| znak P₂ ze zgadywanki `ridge` | ujemne R czytane jako „brak ruchu", 31 fałszywych kandydatów | kryterium musi być na module; znak jest fizyczny |
| swobodny fit P₃ | J2053-7200 dopasowało 63 zamiast 3.06, ρ = 0.12 | maksimum ucieka w róg dużych P₂/P₃; ortogonalizacja względem stałej nie wystarcza |
| iloraz zamiast kąta | R = 4.03 i 3.40, cztery NaN-y, konieczność dwóch progów | tangens ma biegun — ρ (sinus) go nie ma |
| batch mieszał dane | 85 pulsarów pełnopasmowych, 430 na 3/16 pasma | różnica czułości ~2.25× w amplitudzie |
| zmiana sygnatury `_travel_maps` | test reversera czytał 2000 zamiast 3.4·10⁻⁶ | rozpakowanie 4-krotki do 2 zmiennych; złapane przez selftest |
| sprawdzenie okna on-pulse po fakcie | 7 pulsarów z nieczytelnym `ArgumentError` z `Cmd` | `nothing` interpolowane do polecenia `pmod` |

**Błędne wnioski, które trzeba odnotować:**

1. Przypisałem ujemne R niezgodności katalogowego P₃, powołując się na medianę `p3_meas/p3_cat`
   = 0.37 w grupie R < 0.3 wobec 1.02 w grupie R > 0.7. Nieuprawnione: `p3_meas` liczone jest z mapy
   A, więc gdy ruchu nie ma, A jest szumem i `p3_meas` traci sens — korelacja **wynika** z niskiego
   R, zamiast je powodować. Właściwą przyczyną był znak.
2. Ogłosiłem, że korelacja ρ z P₂fit/M unieważnia klasyfikację. Przeceniałem: łańcuch jest
   etykieta → czy istnieje znajdywalna geometria → ρ, czyli metoda działająca. Kontrola syntetyczna
   (§8.2) pokazała, że tłumienie przy szerokim P₂ jest stopniowane i policzalne.
3. Ograniczenie do reżimu P₂ ≤ M/2 przy budowie list kandydatów — za ostre, patrz §10 pkt 3.
4. Przedstawiłem ρ > 1 jako wyjątkową kategorię na podstawie dwóch obiektów z pilota. W pełnej
   próbce jest ich 25, głównie dryferów — to zaludniony ogon, nie anomalia.

Wspólny mianownik: **sam pomiar ρ per pulsar nie zmienił się ani razu**; zmieniała się wyłącznie
ocena, kiedy wolno go interpretować. Każda korekta wyszła z testu, nie z rozumowania.

---

## 12. Użycie

```julia
# pojedynczy pulsar
SpaTs.travel_test(vpmout*"J0820-1350"; max_lag=15, show_=false)

# z dopasowaniem geometrii (|P2| skanowane, P3 z LRFS)
SpaTs.travel_test(vpmout*"J0820-1350"; max_lag=15, p2_template=:auto, show_=false)

# pelnopasmowy plik z katalogu _16
SpaTs.travel_test(vpmout*"J0601-0527_16"; datafile="pulsar_full_debase.txt", show_=false)

# weryfikacja fundamentow
julia --project=. -e 'include("modules/travel.jl"); Travel.selftest()'
```

Argumenty: `max_lag` (domyślnie 40; przy znanym P₃ dobre 2–3·P₃), `max_dphi` (domyślnie M/2),
`hp_halfwin` (50), `nblocks` (4), `p2_template=:auto` z `p3_template` z LRFS, `p2_cap_frac`,
`orth_even`.

Skrypty wsadowe (`~/claude/work/scripts/`): `travel_batch_full.jl` (533 pulsary, wznawialny, tryb
pilotażowy `--limit N --out PLIK`), `travel_rho_summary.jl`, `travel_rho_pilot.jl`,
`travel_variants.jl`, `check_onpulse.jl`, `travel_check.jl`, `travel_stress.jl`.

Wyniki: `~/output/claude/travel_batch_full.csv`, `onpulse_check.csv`, `travel_rho_all.png`,
`travel_rho_pilot.png`.

### Jak czytać wynik

1. **Kontrola off-pulse** — jeśli \|σ\| > 3 dla T lub T_inc, szum nie jest symetryczny w czasie i
   reszta nie znaczy nic.
2. **Detekcja** — max(T, T_inc) ≥ 5σ. Bez niej ρ mierzy czułość, nie pulsara.
3. **Stabilność okna** — ρ przy kilku szerokościach musi być zgodne (§9.3). Bez tego pojedynczy
   obiekt nie jest zmierzony.
4. **Spójność blokowa** — oba znaki to reverser; wartości ≈ 0 przy istotnym T_inc to niezależne
   mody, nie dryf (§10 pkt 4).
5. **\|ρ\|** — blisko 1: koherentna modulacja wędruje; blisko 0: nie wędruje; odniesienie dla
   szerokiego P₂ w §8.2. Znak ρ to kierunek dryfu (dodatni = od wcześniejszych do późniejszych
   długości, Szary+2022).
