# Automatyczne wycinanie złych impulsów (pazur2)

Uruchamiane przez `autozap_all` w `spats.jl`. Wyniki trafiają do `/home/psr/data/OUTPUT/pazur2/<psr>/`.

## Co jest wycinane

1. **Nulle i wygaszenia** (`auto_zap_core`): odcinki bez wykrywalnej emisji. Krótkie przygaśnięcia z widocznym sygnałem zostają; tniemy je tylko wtedy, gdy trwają ≥ 30 impulsów.
2. **Zaszumione odcinki i zły początek obserwacji**.
3. **Pojedyncze śmieci**: paski i skoki szumu w jednym impulsie.
4. **RFI** (`rfi_rms`, reguła z potoku MeerKAT TPA, wersja lokalna): rms poza impulsem odstaje od mediany ~30 sąsiednich impulsów o więcej niż 3σ.

Działają przy tym dwa bezpieczniki:

- **bramka S/N** (`profile_gate`): wycięty odcinek wraca, jeśli z nim S/N profilu jest wyższe;
- **limit dla RFI**: wycięcia RFI nie mogą obniżyć S/N profilu o więcej niż 1%.

Off-pulse liczony jest bez emisji ze średniego profilu (interpuls, składowe poza oknem).

## Klucze w `params.json`

| klucz | co zawiera | kiedy używać |
|---|---|---|
| `zaps` | wycięcia Macieja (nigdy nie zmieniane) | — |
| `zaps2` | wszystkie moje wycięcia (RFI + nulle); wycięte w `pulsar.spCF` | profil, S/N, P3, LRFS/2DFS, foldy |
| `zaps2_rfi` | tylko RFI i śmieci | frakcja nulli, rozkład energii, wszystko, co dotyczy nullingu |
| `zaps2_null` | tylko nulle | analiza samych nulli |
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
