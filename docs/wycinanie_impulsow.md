# Automatyczne wycinanie złych impulsów (pazur2)

Uruchamiane przez `autozap_all` w `spats.jl`. Wyniki trafiają do `/home/psr/data/OUTPUT/pazur2/<psr>/`.

## Co jest wycinane

1. **Nulle i wygaszenia** (`auto_zap_core`): odcinki bez wykrywalnej emisji. Krótkie przygaśnięcia z widocznym sygnałem zostają; tniemy je tylko wtedy, gdy trwają ≥ 30 impulsów.
2. **Zaszumione odcinki i zły początek obserwacji**.
3. **Pojedyncze śmieci**: paski i skoki szumu w jednym impulsie.
4. **RFI** (`rfi_rms`, reguła z potoku MeerKAT TPA, wersja lokalna): rms poza impulsem odstaje od mediany ~30 sąsiednich impulsów o więcej niż 3σ.

Działają przy tym dwa bezpieczniki:

- **bramka S/N** (`profile_gate`): wycięty odcinek wraca, jeśli z nim S/N profilu jest wyższe. Potem bramka oddaje dobre krawędzie odcinków, które zostały wycięte, jeśli te krawędzie mają własny sygnał (≥3σ) i podnoszą S/N. Dzięki temu kilka fatalnych impulsów nie ciągnie za sobą całego odcinka;
- **limit dla RFI**: wycięcia RFI nie mogą obniżyć S/N profilu o więcej niż 1%, liczone zarówno po wycięciu nulli, jak i na danych z samym RFI (nulle zostawione, `rfi_ref=:both`).

Off-pulse liczony jest bez emisji ze średniego profilu (interpuls, składowe poza oknem).

## Klucze w `params.json`

| klucz | co zawiera | kiedy używać |
|---|---|---|
| `zaps` | wycięcia Macieja (nigdy nie zmieniane) | — |
| `zaps2` | wszystkie moje wycięcia (RFI + nulle); wycięte w `pulsar.spCF` | profil, S/N, P3, LRFS/2DFS, foldy |
| `zaps2_rfi` | tylko RFI i śmieci | frakcja nulli, rozkład energii, wszystko, co dotyczy nullingu |
| `zaps2_null` | tylko nulle | analiza samych nulli |
| `zaps2_other` | wygaszona emisja i marginesy przy nullach (od automatu v4) | wycięte do profilu, ale ani RFI, ani null |
| `zap_cmp` | S/N profilu: bez wycięć / moje / tylko RFI / Macieja; liczby impulsów | porównania |

W kodzie odpowiadają im funkcje `Data.zaps_of(p)`, `Data.zaps_rfi_of(p)` i `Data.zaps_null_of(p)`.

**Uwaga:** `pulsar.spCF` w pazur2 ma wycięte `zaps2`, czyli również nulle. Do badań nullingu trzeba zbudować archiwum od nowa z surowych danych, z wycięciem tylko `zaps2_rfi`.

## Pliki i weryfikacja

**Pliki dla każdego pulsara:**

- `pulsar_autozap.pdf`: co i dlaczego wycięto (kolory powodów);
- `pulsar_contrast.pdf`: dane po wycięciu (wycięte impulsy białe) dla wariantów przed / automat / Maciej, z profilami.

**Pliki zbiorcze:**

- `zap_summary.csv`: zawiera kolumny `n_rfi` i `n_null`;
- `do_recznie.csv`: pulsary bez okna on-pulse lub bez P3.

**Weryfikacja:** `julia tests/zap_benchmark.jl` uruchamia benchmark na danych syntetycznych ze znaną prawdą. Porównuje automat z regułą TPA i podejściem w stylu clfd.

## Automat v4: czyste nulle

Do listy `zaps2_null` trafiają tylko impulsy bez sygnału.

- Odcinek oznaczony jako null, w którym średni sygnał jest istotny (>3σ) i wynosi co najmniej 15% typowego impulsu, to:
  - **zaszumiony odcinek** (RFI), jeśli szum w nim jest wyraźnie podwyższony; tak było w J1812+0226;
  - **wygaszona emisja** (powód 6) w pozostałych przypadkach.
- Marginesy przy nullach nie liczą się jako nulle, bo impulsy na brzegach często mają część emisji.
- Wycięcia i S/N się nie zmieniają (benchmark syntetyczny bez zmian); zmienia się tylko podział na `zaps2_rfi` / `zaps2_null` / `zaps2_other`.
- Test syntetyczny: prawdziwe nulle zostają nullami we wszystkich przypadkach; odcinek z 30% emisji w jasnym pulsarze zawsze dostaje etykietę wygaszonej emisji.

## Weryfikacja na prawdziwych danych

`zap_verify()` (`modules/zap_verify.jl`, uruchamiane z `main()`, wznawialne). Wyniki w `/home/psr/data/OUTPUT/zap_verify/`:

- `summary.txt`: podsumowanie do wklejenia;
- `nulls.pdf`: profil średni nulli vs dobrych impulsów;
- `verify.csv`: wyniki dla każdego pulsara.

Sprawdza:

- spójność z zapisanymi wynikami;
- sygnał w nullach i frakcję nulli: nasza vs rozkład energii vs literatura;
- nowy bezpiecznik RFI vs stary;
- kanały `.spCf16`: czy wycięcie złych kanałów przed sumowaniem pomaga.

Przy okazji zapisuje `pazur2/<psr>/pulsar_rfi.spCF` z wyciętym tylko RFI, do badań nullingu.

## Stan (2026-11-17, automat v4)

Zrobione: kanały przed sumowaniem (24 pulsary), weryfikacja nulli (145 ze 183 pulsarów z czystymi nullami),
test MeerGuard (tylko do sprawdzania kanałów), model frakcji nullingu (`tools/nulling_energies.jl` +
`tools/nulling_model.py`, zgodność z literaturą −0,2 pp), przegląd literatury, diagnostyka (`tools/zap_checks.jl`).
Porównanie z metodą ręczną: 262 : 6 : 256 (524 pulsary).

Dokończenie po pełnym przebiegu: `tools/finish_v4.jl` (diagnostyka pulsarów z lepszą metodą ręczną,
kanały dla długich obserwacji, raport modelu nulli z oznaczeniem „model nie rozróżnia”, porównanie).

Frakcja nullingu: brać `NF_model` z `nulling_model_v4/summary.csv` tylko przy `reliable` = prawda.
Liczba wyciętych nulli automatu nie jest frakcją nullingu (automat zostawia krótkie nulle).

Dalej: P3 i dryf (P3-fold, klasyfikator dryf/modulacja amplitudy, diagram P–Pdot).

## Do zrobienia (decyzja 2026-11-18)

- **Nie wycinać marginesu przy nullach.** Teraz każdy wycięty odcinek dostaje po 1 impulsie z każdej
  strony (`auto_zap_core(pad=1)`, powód 5 „margines”, od v4 w `zaps2_other`). Przy nullach nie jest to
  potrzebne: impuls na brzegu nulla z częścią emisji to nadal emisja. Margines zostawić tylko przy
  odcinkach RFI / zaszumionych (brzegi zakłócenia). Przed zmianą: skan wszystkich pulsarów obiema
  wersjami (jak `tools/rfi_v5.jl`) i test syntetyczny; zmieni się tylko profil (nulle bez zmian).

## Stan v6 (2026-11-18)

Przyjęty wariant **automat + P≥0.95**: nulle automatu (jak v5) plus impulsy, które model energii
(TPA XIX) uznaje za pewne nulle (P_null ≥ 0.95, tylko pulsary z wiarygodnym modelem). Wybrany w
`tools/final_v6.jl` spośród 8 wariantów (w tym jednej reguły z publikacji, P_null ≥ 0.5, która wycina
3% emisji i zniekształca profil o 3.6%). Syntetycznie: kompletność nulli 97.9%, bias profilu ≤ 0.7%;
prawdziwe dane: lepszy od metody ręcznej w 336 pulsarach, gorszy w 2. Dynesty (167 pulsarów) potwierdza
NF modelu: mediana różnicy +0.1 pp. Domyślnie w kodzie: `Data.NULL_MODE[] = :union`, `NULL_THR[] = 0.95`,
lista nulli z `NULL_MODEL_DIR` (final_v6/model/). Powrót do v5: `Data.NULL_MODE[] = :auto`.
Kopia poprzedniej wersji: `/home/psr/data/OUTPUT/pazur2_przed_v6`.

## Stan v7 (2026-11-18)

`tools/final_v7.jl`: każdy krok automatu wyłączany osobno (12 wariantów), dane syntetyczne
(`tests/zap_ablation.jl`) i wszystkie pulsary. Jedna reguła: krok zostaje, jeśli bez niego cokolwiek się
pogarsza. Zbędne okazało się tylko wypełnianie krótkich wysp między cięciami (bez zmian na obu rodzajach
danych) — wyłączone na stałe (`auto_zap_core(min_keep=0)`). Skoki szumu i pojedyncze nulle są zbędne
syntetycznie, ale potrzebne na prawdziwych danych (8 i 14 pulsarów gorzej), więc zostają.

Linia bazowa GP (jądro SE na okręgu fazy, siatka ell × λ) jako dodatkowa opcja w energiach, wybierana tą
samą regułą co średnia/Fourier (najmniejszy rozrzut E_off): wybrana w 129 z 514 pulsarów, tam rozrzut
mediana 0.92× (najlepiej 0.62×). Zgodność z TPA XIX bez zmian (|różnica| mediana 0.21 pp wobec 0.22,
w 2σ 94.8% wobec 95.3%, wspólnych wiarygodnych 364 wobec 358). Model nulli: `final_v7/model/`.

Nulle w 16 kanałach: nulle modelu (P≥0.95, 99 pulsarów) bez emisji w żadnym kanale; nulle automatu: emisja
w jednym kanale tylko w J0034−0721 (kanał 6, 11%) i J1104−6103 (kanał 7, 46%) — do sprawdzenia.

Wynik końcowy: lepszy od metody ręcznej w 337, gorszy w 2, remis 192; średni zysk S/N +18.3%, wycięte
7.6%. Kopia v6: `/home/psr/data/OUTPUT/pazur2_przed_v7`.

## Stan v8 (2026-11-25)

Nulle wycinane jedną regułą z publikacji (TPA XIX, Keith+2026): w pulsarach z wiarygodnym modelem energii null
to impuls z P_null ≥ 0,5; w pulsarach bez wiarygodnego modelu nulle wycina automat jak w v7
(`Data.NULL_MODEL_FALLBACK[] = :auto`; `:none` = bez nulli). Reszta automatu (RFI, kanały, bezpiecznik,
margines przy RFI) bez zmian. Domyślnie w kodzie: `NULL_MODE = :model`, `NULL_THR = 0.5`, model nulli z
`final_v7/model/`. Wybrane jako prostsze (jedna reguła zamiast kilku) mimo znanego kosztu: w ok. 1 na 5
pulsarów z nullami reguła wycina też słabą emisję (profil jaśniejszy o >1% w 33 z 174, najgorzej +22%).

Wynik (531 pulsarów): lepszy od metody ręcznej w 346, gorszy w 2 (J1843+2024, J2111+2106), remis 183;
średni zysk S/N +19,1%, wycięte 8,9%; pulsary z nullami 247 (czyste <2% energii: 189, domieszka 2–10%: 41,
>10%: 17, głównie po kilka nulli); NF a literatura: mediana +0,5 pp. Kopia v7: `/home/psr/data/OUTPUT/pazur2_przed_v8`.

Marginesy w v8 (decyzja 2026-11-25: zostają): przy nullach z modelu brak marginesu; przy RFI, zaszumionych
odcinkach i śmieciach margines 1 impulsu; przy nullach automatu (pulsary bez modelu) margines zostaje
(`NULL_MARGIN[] = true`). Usunięcie marginesu wszędzie (`pad=0`) pogorszyło 59 pulsarów i wykrywanie RFI (−32 pp).
Punkt „Nie wycinać marginesu przy nullach” z 2026-11-18 zamknięty.

## Stan v9 (2026-11-26): jedna zasada dla wszystkich pulsarów

Nulle wycinamy tylko tam, gdzie model energii (TPA XIX) pewnie odróżnia je od emisji: w pulsarach z wiarygodnym
modelem null to impuls z P_null ≥ 0,5; w pozostałych (model niepewny albo brak modelu) nulli nie wycinamy wcale.
RFI bez zmian (z marginesem 1 impulsu). Domyślnie w kodzie: `NULL_MODE = :model`, `NULL_THR = 0.5`,
`NULL_MODEL_FALLBACK = :none`, `NULL_MODEL_ANY = false`.

Wybór (`tools/one_rule.jl`, wszystkie pulsary, te same dane):

| opcja | wobec metody ręcznej | średni zysk S/N | wycięte | uwagi |
|---|---|---|---|---|
| v8: model pewny + automat | 345 : 2 : 184 | +20,4% | 8,1% | dwie metody |
| B: model wszędzie | 343 : 4 : 184 | +20,4% | 8,5% | w słabych pulsarach wycina prawdziwą emisję (J1549: 71% impulsów) |
| **3: model tylko pewny (przyjęta)** | 333 : 6 : 192 | +19,2% | 7,3% | jedna zasada, nic nie fałszuje; nulle zostają w 40 słabszych pulsarach (J1931 −41%) |
| automat dla wszystkich | 323 : 2 : 206 | +18,2% | 5,6% | znajduje tylko 40% nulli |
| impuls E < 3σ / seria 5 < 3σ | 231 : 220 / 276 : 117 | +6,7% / +11,3% | 43% / 26% | wycinają dużo emisji |

Margines przy RFI zostaje: bez niego 31 pulsarów traci >1% S/N (J1148 −9,4%). Budowa i przyjęcie:
`one_rule_build()`, `one_rule_adopt()`, potem `one_rule_replot()`. Kopia v8: `pazur2_przed_v9`.

### v9.1 (2026-11-26): silne RFI zawsze

Impuls z silnym RFI (z > 6: szum poza impulsem odstaje o ponad 6σ od mediany ~31 sąsiadów) wycinamy zawsze, bez
bezpiecznika S/N (`Data.RFI_HARD[] = 6.0`; `Inf` = jak v9). Powód: poziome paski (przesunięcie tła na całym impulsie)
nie psują S/N średniego profilu, więc bezpiecznik ich nie przepuszczał (J0932−3217: z ≈ 13–24 zostawione), a psują
energie, nulle i LRFS. Porównanie z listą wyciętych impulsów TPA VIII (`stats.csv`, `tools/tpa_zap_compare.jl`):
przed poprawką z 1551 impulsów TPA wycięliśmy 66%, 508 zostało (część z z ≫ 3). Słabsze RFI (3 < z ≤ 6) jak dotąd,
z bezpiecznikiem. Kopia v9: `pazur2_przed_v9_2`.
