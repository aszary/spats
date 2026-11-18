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
