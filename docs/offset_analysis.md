# Offsety składowych 1023 → 1523 MHz — wyniki analizy

Stan na 2026-09-09. Dane: `input/offsets.csv` (161 pulsarów), złączone z `input/psrcat.db`
po nazwie PSRJ. Wielkość analizowana to to samo, co liczy `Plot._read_offsets`: dla pulsara
o ≥ 2 składowych **zmiana separacji** między skrajnymi komponentami,
`Δsep = offset_ostatni − offset_pierwszy`, dla pulsara jednoskładnikowego surowe przesunięcie
tej jednej składowej. Offset dodatni = składowa przesuwa się w stronę późniejszych faz
przy 1523 MHz.

Wszystkie liczby poniżej pochodzą z jednego przebiegu skryptu weryfikacyjnego (opis
odtworzenia na końcu) — nie z odczytu z wykresu.

## Próbka

Po cięciu `grade ≥ 6` zostaje **159 ze 161** pulsarów (odpadają J1714−1054 i J2139+2242):

| | liczba |
|---|---|
| wielokomponentowe (Δ separacji) | 105 |
| jednokomponentowe (przesunięcie) | 54 |
| rozkład składowych | 54× n=1, 82× n=2, 22× n=3, 1× n=4 |

Zakres: `P` od 0.048 do 7.73 s, `Ṗ` od 1.1e−19 do 1.4e−13. **Zero pulsarów milisekundowych** —
cała próbka leży w głównej wyspie populacji, co ogranicza wnioski do pulsarów normalnych.

## 1. Profile zwężają się z częstotliwością

Efekt populacyjny, nie pojedyncze przypadki:

- `Δsep < 0` dla **74 ze 105**, test dwumianowy **p = 1.9e−5**
- mediana **−0.463°**, przedział bootstrapowy 68%: [−0.517, −0.358]
- Wilcoxon p = 1.4e−6

To klasyczne radius-to-frequency mapping, tu zmierzone jednorodnie na jednej próbce
dwuczęstotliwościowej zamiast składane z klasyfikacji morfologicznych.

## 2. To nie jest błąd DM ani wyrównania pasm

Rozbicie Δsep na ruch obu krawędzi profilu osobno:

| | mediana | znaki | p |
|---|---|---|---|
| składowa wiodąca | **+0.220°** [+0.163, +0.290] | 74/105 dodatnich | 3.3e−5 |
| składowa tylna | **−0.086°** [−0.143, −0.056] | 68/105 ujemnych | 3.2e−3 |
| centroid (wspólne przesunięcie) | +0.065° [+0.042, +0.089] | 61/44 | 0.12 (n.s.) |

Obie krawędzie ruszają się do środka **niezależnie**, a wspólne przesunięcie jest zgodne
z zerem. Błąd dedyspersji albo złe wyrównanie faz między pasmami przesunęłyby oba komponenty
w tę samą stronę — dane to wykluczają.

Zwężenie jest przy tym symetryczne: test Wilcoxona na `|lead| − |trail|` daje p = 0.51.
Różnica median (0.220 vs 0.086) to ogon rozkładu, nie systematyka. Wniosek uboczny: **brak
sygnału aberracji/retardacji** na poziomie ~0.05°.

## 3. Siła zwężania spada z Ė

Po zdjęciu skalowania rozmiaru wiązki (`W ∝ P^−1/2`, stąd wielkość `Δsep·√P`), kwintyle
spin-down luminosity:

| log Ė | n | mediana Δsep·√P | zwężeń >3σ | poszerzeń >3σ | nieistotnych |
|---|---|---|---|---|---|
| 28.9 – 30.9 | 21 | −1.105 | 16 | **0** | 5 |
| 31.0 – 31.4 | 21 | −0.638 | 16 | 1 | 4 |
| 31.4 – 32.1 | 21 | −0.322 | 10 | 5 | 6 |
| 32.1 – 33.1 | 21 | −0.055 | 8 | 8 | 5 |
| 33.1 – 35.9 | 21 | −0.240 | 9 | 4 | 8 |

Częściowa korelacja Spearmana `Δsep` vs `log Ė` przy kontrolowanym `log P`, na podpróbkach:

| podpróbka | n | ρ | p |
|---|---|---|---|
| wszystkie wielokomponentowe | 105 | +0.295 | 0.0024 |
| tylko >3σ | 77 | **+0.411** | 2.3e−4 |
| tylko n = 2 | 82 | +0.405 | 1.8e−4 |
| tylko grade = 10 | 75 | +0.374 | 0.0010 |
| po obcięciu \|Δsep\| < 3° | 96 | +0.207 | 0.044 |

Trend przeżywa każde cięcie i jest **najsilniejszy na najczystszych podpróbkach**, co jest
zachowaniem oczekiwanym od sygnału, a nie od artefaktu.

Kierunek zgadza się z literaturą — empiryczny model emisji Karastergiou & Johnston (2007)
wiąże zakres wysokości emisji z Ė, a pulsary o wysokim Ė mają emitować z wąskiego przedziału
wysokości, przez co RFM przestaje działać. Warto natomiast odnotować, że to **nie jest
przewidywanie trywialne**: standardowa formuła na wysokość emisji (Kijak & Gil) daje
zależność od częstotliwości praktycznie uniwersalną, `r ∝ ν^−0.26`, czyli stałą frakcjonalną
zmianę szerokości dla wszystkich pulsarów. Trend rozstrzyga między tymi klasami modeli.

### Czego natomiast literatura nie przewiduje

Przy wysokim Ė **znak się rozjeżdża, a nie amplituda maleje**. Osłabienie RFM przewiduje
`Δsep → 0`, czyli migrację pomiarów do słupka „nieistotne". Tymczasem liczba nieistotnych
rośnie tylko z 5 do 8, a pojawiają się wysokoistotne poszerzenia — w czwartym kwintylu
8 zwężeń i 8 poszerzeń. W najniższym kwintylu Ė poszerzeń nie ma **ani jednego**.

## 4. Poszerzenia

18 pulsarów pokazuje `Δsep > 0` na ponad 3σ (przy 59 zwężeniach na tym samym progu):

| pulsar | Δsep | istotność | n | Ė |
|---|---|---|---|---|
| J1822−4209 | +2.573 ± 0.269 | 9.6σ | 2 | 1.9e32 |
| J1843−0459 | +1.978 ± 0.381 | 5.2σ | 2 | 7.8e31 |
| J1757−2421 | +1.700 ± 0.228 | 7.5σ | 3 | 4.0e34 |
| J1834−0426 | +1.365 ± 0.167 | 8.2σ | 4 | 1.2e32 |
| J1557−4258 | +1.151 ± 0.188 | 6.1σ | 3 | 3.7e32 |
| J1847−0402 | +1.074 ± 0.035 | **30.3σ** | 2 | 9.6e33 |

Najgroźniejszym wyjaśnieniem było przekręcone parowanie składowych między pasmami —
`GaussianFit.component_offsets` paruje je czystym sortowaniem po μ, więc zmiana względnych
amplitud mogłaby je zamienić miejscami. Ta hipoteza jednak nie tłumaczy danych, patrz
test monotoniczności niżej.

## 5. Testy kontrolne

### 5.1. Jednokomponentowe jako próba kontrolna

Pulsary o jednej składowej nie mają separacji, którą można zmienić — jakikolwiek trend z Ė
w ich przesunięciach musiałby być systematyką pomiarową (dedyspersja, wyrównanie pasm, S/N,
zachowanie fitu):

| kanał | n | ρ (częściowa, \| log P) | p |
|---|---|---|---|
| wielokomp.: Δsep | 105 | +0.295 | 0.0024 |
| wielokomp.: wspólne przesunięcie | 105 | +0.180 | 0.067 |
| **jednokomp.: przesunięcie** | 54 | **−0.006** | **0.965** |

Kanał kontrolny jest płaski. Trend siedzi wyłącznie w separacji, nie w pozycji profilu —
to wyklucza całą klasę wyjaśnień systematycznych naraz.

### 5.2. Podłoga systematyczna

Rms przesunięć jednoskładnikowych = **0.707°**. Tyle wynosi rozrzut przy *zerowej* zmianie
separacji, więc to realistyczna miara błędu pojedynczego pomiaru. Rms Δsep = 2.027°,
F = 8.2, p = 4.3e−14 — sygnał leży 8× ponad podłogą wariancji.

Dla porównania mediana deklarowanego błędu to 0.087°, czyli **osiem razy mniej niż realna
podłoga**. Globalne χ²/dof wokół średniej ważonej wynosi 97.

### 5.3. Uporządkowanie deformacji (n ≥ 3)

Dla 23 pulsarów o ≥ 3 składowych można sprawdzić monotoniczność offsetów wzdłuż profilu —
jednorodne skalowanie wymaga, żeby offsety zmieniały się monotonicznie z numerem składowej:

- **11 monotonicznie malejących** (ściskanie), oczekiwane przypadkiem 3.7 → **p = 3.9e−4**
- 6 monotonicznie rosnących (rozciąganie)
- 6 niemonotonicznych

Razem 17/23 monotonicznych przy 7.4 oczekiwanych. Offsety to nie jest szum per składowa.

Sześć rozciągnięć wygląda tak:

| pulsar | offsety składowych | Δsep |
|---|---|---|
| J1757−2421 | −1.154, +0.267, +0.546 | +1.700 ± 0.228 |
| J1328−4921 | −1.505, −0.299, +0.038 | +1.543 ± 0.553 |
| J1527−5552 | −0.451, +0.025, +0.385 | +0.836 ± 0.136 |
| J1017−5621 | −0.606, +0.118, +0.141 | +0.747 ± 0.199 |
| J1807−0847 | −0.015, +0.165, +0.407 | +0.421 ± 0.077 |
| J1057−5226 | +0.062, +0.137, +0.139 | +0.076 ± 0.040 |

To są gładkie, monotoniczne rampy. **Złe sparowanie składowych rozsypałoby offsety, a nie
ułożyło je w rampę** — poszerzenia przeżywają swój najgroźniejszy zarzut.

### 5.4. Zlanie składowych przy wysokim Ė

Hipoteza „przy wysokim Ė fituje się substrukturę szerokiego zlanego profilu, a nie parę
stożkową" przewiduje większe błędy przy wysokim Ė. Dane pokazują odwrotność: częściowa
korelacja `log(err)` vs `log Ė` przy kontrolowanym `log P` wynosi **−0.291 (p = 0.0028)**,
czyli błędy przy wysokim Ė są *mniejsze*. Poszerzenia i zwężenia mają identyczne rozkłady
błędów (Mann-Whitney p = 1.000).

## 6. Ograniczenia

- **Wielkość jest w stopniach absolutnych, nie frakcyjna.** To główne otwarte zastrzeżenie,
  patrz sekcja 7.
- **Rozrzut populacyjny 20× przewyższa deklarowane błędy** (rms 1.73° vs mediana błędu
  0.087°, χ²/dof ≈ 97). Nie da się z tych danych rozstrzygnąć, czy to realny rozrzut
  międzyźródłowy, czy niedoszacowane błędy. W obu przypadkach średnia ważona i jej formalna
  istotność są bezużyteczne — wszystko powyżej policzone na medianach i testach rangowych.
- **Δsep dla n ≥ 3 gubi wnętrze profilu.** `_read_offsets` bierze tylko skrajne składowe,
  więc 22 pulsary z trzema i 1 z czterema składowymi oddają tylko część informacji.
- **Próg 3σ działa niespójnie między typami.** Wszystkie 54 pomiary jednoskładnikowe są >3σ
  (mediana błędu 0.031°), wśród wielokomponentowych tylko 77/105 (mediana 0.129°, bo błędy
  się dodają). Puste symbole na diagramie P–Ṗ to z definicji wyłącznie kółka.
- **Wzbogacenie poszerzeń w profile złożone** (n ≥ 3: 6/18 vs 10/59) idzie w podejrzanym
  kierunku, ale Fisher p = 0.18 — nie rozstrzyga.
- **Rozciągnięcia vs ściskania nie różnią się w Ė** (Mann-Whitney p = 0.40), ale to 6 vs 11
  obiektów — test nic nie może wykryć.

## 7. Po co ΔW/W

Mierzona wielkość faktoryzuje się na dwie rzeczy:

```
Δsep [°]  =  (ΔW/W)  ×  W
             fizyka      geometria
             ile RFM     jak szeroka wiązka
```

Trend z Ė może siedzieć w którymkolwiek z czynników, a dostępny jest tylko iloczyn.
Mnożenie przez `√P` i kontrolowanie `log P` w korelacji cząstkowej załatwia człon wiodący,
ale nie to, że przy ustalonym P szerokość profilu i tak zmienia się o czynnik kilku (stożek
vs rdzeń, różne kąty β, różne wysokości emisji). Ten rozrzut wchodzi jako nieuwzględniony
szum wprost do korelowanej wielkości.

Gorzej: **jeśli samo W koreluje z Ė przy ustalonym P**, trend może być w całości efektem
geometrycznym przy uniwersalnym RFM. Kierunek działa na korzyść obecnego wniosku — młodsze
pulsary mają raczej szersze profile, więc uniwersalne ΔW/W dałoby im *większe* Δsep
w stopniach, czyli odwrotnie niż obserwujemy — ale to argument jakościowy o wielkości,
której się nie mierzy.

Normalizacja usuwa ten problem z definicji i dokłada trzy rzeczy:

1. **Porównywalność z literaturą.** ΔW/W przelicza się na wykładnik `W ∝ ν^−a`. Uwaga na
   pułapkę: `longitude` z `component_offsets` jest środkowopasmowe, więc `frac = ΔW/W_mid`,
   a nie `ΔW/W_low`. Poprawnie `W_high/W_low = (1 + frac/2)/(1 − frac/2)`. Dla przykładowego
   pulsara z `frac = −0.176 ± 0.026` daje to `a = 0.44 ± 0.06`, przy typowych z literatury
   0.2–0.3. Naiwne potraktowanie `frac` jako `ΔW/W_low` myli wykładnik o 0.04.
2. **Modele przewidują ΔW/W, nie Δsep.** Obraz Karastergiou & Johnston mówi o frakcjonalnej
   ewolucji profilu — test staje się bezpośredni zamiast pośredniego.
3. **Kalibracja istotności.** Podłoga 0.71° jest w stopniach. Na profilu szerokim 30° to 2%,
   na profilu 5° to 14%, a dziś oba pomiary są traktowane jednakowo.

### Stan implementacji

Podpięte w commicie `184c3ba`. `Plot._offset_summary` (`modules/plot.jl`) liczy średnią
ważoną pozycji składowych, separację i ΔW/W, i dopisuje wiersz do `input/separations.csv`
(podmienia istniejący wiersz dla tego samego pulsara). Wołane z `analyse_p3folds4`
i `analyse_average_offset`; nazwa pulsara wyciągana z katalogu danych przez
`Data.psr_from_dir`. Szczegóły rachunku: longitude to `(μ_high + μ_low)/2`, więc jej błąd
to połowa błędu offsetu i te same wagi obsługują obie średnie; błędy separacji kwotowane
jako σ_ext (= rms/√n), nie σ_int.

**Pozostaje przelecieć pulsary ponownie**, żeby `input/separations.csv` się zapełniło —
dopiero wtedy zastrzeżenie z sekcji 7 da się domknąć.

## Odtworzenie

Analiza nie jest częścią pipeline'u — to jednorazowy skrypt czytający `input/offsets.csv`
i `input/psrcat.db` tak samo, jak robi to `Plot._read_offsets` (łącznie z mapą
`PSR_RENAMED`: J1402−5124 → J1402−5021). Statystyka: `scipy.stats` — testy dwumianowe
i Wilcoxona na znakach, Spearman i korelacje cząstkowe na rangach, przedziały medianowe
z bootstrapu 20000 prób.

Diagram P–Ṗ z kolorem kodującym offsety: `Plot.ppdot_offsets("output")`. Skala koloru jest
percentylowa (`offset_vmax_quantile = 0.90`, czyli ±2.2°) — użycie maksimum oddawało zakres
jednemu pulsarowi z komentarzem „need to be redone".
