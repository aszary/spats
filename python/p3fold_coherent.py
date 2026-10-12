#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p3fold_coherent.py — samodzielny port (Python 3, numpy/scipy/matplotlib) foldu koherentnego P3 ze spats
=========================================================================================================

Źródło (Julia, repo github.com/aszary/spats, gałąź `claude`):
  * `spats.jl`: `SpaTs.p3fold_coherent` (wersja pierwotna) i `SpaTs.p3fold_coherent_agent` (wersja bieżąca),
  * `modules/p3fold_viterbi.jl`: `coherent_fold`, `coherent_fold_jackknife`, `coherent_fold_agent`,
    `auto_cutoff_agent`, `kernel_fold_agent`, `kernel_sigma_agent` i funkcje pomocnicze,
  * `modules/data.jl` (`load_ascii`, `zap!`), `modules/tools.jl` (`p3fold`, `intensity_pulses`, `read_params`),
    `modules/plot.jl` (`p3fold_compare`).
Dokumentacja i historia testów: `docs/coherent_fold_params.md`, dziennik `docs/separations_analysis_log.md`.

Uruchomienie
------------
    python3 p3fold_coherent.py <katalog_z_params.json> [--datafile pulsar.debase.txt] [--plotdir .]
            [--name nazwa] [--legacy] [--cutoff 0.01] [--threshold-q 0.3] [--kernel-sigma 0.07]
            [--no-split-nulls] [--fold bins]

Katalog musi zawierać `params.json` (klucze: p3, p3_ybins, bin_st, bin_end, opcjonalnie zaps) i plik z pojedynczymi
pulsami w formacie ASCII PSRCHIVE (`pdv -t`: pierwsza linia nagłówek, potem kolumny: puls, kanał, bin, natężenie).
Wynik: `<plotdir>/<name>_p3fold_compare.pdf/.png` (fold koherentny | fold stałym P3 | energia pulsów | P3(n)).
W Pythonie: `res = p3fold_coherent_agent(outdir, ...)` zwraca słownik z foldem, fazą, P3(n) itd.

Zależności: numpy, scipy (signal), matplotlib. Sprawdzone ze scipy 1.8 i 1.15.

Różnice względem wersji w Julii (wyniki zgodne co do wartości, nie co do bitu)
--------------------------------------------------------------------------------
  * Generator liczb losowych: tasowania pulsów (null w `auto_cutoff`, próg |s|) używają numpy, nie MersenneTwister
    z Julii — inne permutacje, więc próg |s| i wynik CV różnią się nieznacznie (zwykle ten sam wybór f_c i σ).
  * `filtfilt`: Butterworth w postaci sekcji drugiego rzędu, odbicie nieparzyste na brzegach jak w DSP.jl;
    długość wypełnienia ustawiona jak w DSP.jl (6 × liczba sekcji), więc brzegi zbliżone, ale nie identyczne.
  * Indeksy: `bin_st`, `bin_end` w params.json są 1-bazowe i włącznie (konwencja Julii) — kod to przelicza.
  * Weryfikacja (2026-10-12): s(n) po filtrze identyczne z Julią do 4 cyfr; S/N, f_c i σ jądra takie same dla
    J1750-3503 (2.1, 1/147, 0.1) i J0152-1637 (2.7, 1/20.9, 0.1). Pokrycie P3(n) różni się (J1750 78% vs 70%,
    J0152 96% vs 90%) wyłącznie przez próg |s| z tasowania (słabość 13).

=====================================================================================================================
METODA (wariant C, `coherent_fold_agent` — domyślna)
=====================================================================================================================
Cel: fold P3 („jak wygląda cykl modulacji”) dla pulsarów, w których P3 wędruje, skacze, są nulle, a pojedynczy puls ma
za mało S/N, żeby przypisać mu fazę na oko. Zamiast stałego P3 każdy puls dostaje własną fazę cyklu.

 1. Szablon: jedna, pełnej długości FFT po pulsach w każdym binie on-pulse; bin częstotliwości k = round(N/P3)
    daje zespolony wzór L(φ) — „jak każdy bin zmienia się z okresem P3”. Przy dryfie |L| to profil modulacji,
    a arg L rośnie wzdłuż długości (pochylone pasmo).
 2. Filtr dopasowany: każdy puls (po odjęciu średniego profilu) rzutowany na conj(L): z(n) = Σ_φ X(n,φ)·L*(φ).
    Kąt z(n) to „w którym miejscu cyklu puls wygląda, że jest”, długość — jak wyraźnie.
 3. Demodulacja przy f3 = 1/P3 i filtr dolnoprzepustowy (Butterworth rzędu 6, zero-fazowy filtfilt) z częstością
    odcięcia f_c: s(n) = LP[z(n)·e^{−i·nośna(n)}]. Faza całkowita = nośna + arg s. Jest to przesuwny odpowiednik
    LRFS bez bloków (brak granic bloków, na których wędrówka P3 niszczy koherencję).
 4. f_c automatycznie (`auto_cutoff_agent`) z siatki f3·{1/16, 1/10, 1/8, 1/6, 1/4, 1/3}: walidacja krzyżowa po
    długości (faza z połowy binów, fold drugiej połowy; połówki = przeplatane bloki po 4 biny), miara ΔR² =
    R²_cv − R²_cv dla przetasowanych pulsów. Odejmowanie tasowania jest konieczne: szybki filtr sortuje pulsy wg
    położenia podpulsu i „wyjaśnia” wariancję także bez okresowości. f_c > rozrzut |1/P3_lok − f3| (żeby nadążyć
    za P3) i f_c ≲ f3/3 (powyżej przecieka składowa −f3 i modulacja natężenia/nulle).
 5. Adaptacyjna nośna: 2 dodatkowe przebiegi, w których nośna = faza z poprzedniego przebiegu (rozwinięta po
    pulsach z |s| ≥ próg, liniowo pomiędzy) — lokalne P3 daleko od P3 nie jest tłumione przez filtr.
 6. Próg |s|: kwantyl `threshold_q` (domyślnie mediana) |s| z danych o przetasowanej kolejności pulsów (poziom szumu).
    Z P3(n) wypadają pulsy poniżej progu, brzegi filtfilt (1/(2 f_c) pulsów) i — przy `split_nulls` — nulle
    wykryte z energii (średnia po 5 P, ułamek nulli = 2·frac(Ē < 0) jak u Ritchingsa 1976, epizody ≥ 2 P;
    uśrednianie konieczne — z pojedynczych pulsów J1750-3503 wychodziło fałszywie 18% nulli).
    Fold używa WSZYSTKICH pulsów (faza jest modulo 2π).
 7. P3(n) z nachylenia fazy całkowitej (regresja ważona |s|²) w oknie max(1/(2 f_c), 3·P3, 30) pulsów, osobno
    w każdym ciągłym odcinku użytych pulsów — bez rozwijania fazy przez przerwę (liczba cykli w przerwie ~P3 jest
    nieznana; J1750: wszystkie skoki > 3 P leżały przy przerwach). Wartości poza [2, 3·P3] → NaN (pojedynczy punkt
    −700 w J1919+0134). Błędy: P3(n) liczone osobno w `n_groups` = 4 rozłącznych zakresach długości (niezależny
    szum), std/√n.
 8. Fold jądrowy (domyślny): dla 64 faz siatki średnia WSZYSTKICH pulsów ważona zawiniętym Gaussem w fazie,
    w = exp(−d²/2σ²), d = różnica faz [cykle]. σ z danych (`kernel_sigma_agent`): faza z połowy binów, każdy puls
    drugiej połowy przewidywany z foldu pozostałych (leave-one-out na siatce 256 faz), siatka σ 0.01–0.3 cyklu,
    wybierane największe σ w 1% od maksimum. Alternatywnie `fold="bins"`: ŚREDNIA (nie suma) pulsów w binie fazy
    — fazy z danych obsadzają biny nierówno (J1539-4828: 38–155 pulsów/bin), suma pokazywałaby obsadzenie.
 9. S/N filtra dopasowanego: |L_on| względem szumu L w binach off-pulse — sprawdzenie, czy jest co śledzić.

Wersja pierwotna `p3fold_coherent` (`--legacy`): krok 1–3 ze stałym f_c = 1/300, bez adaptacyjnej nośnej, progu,
nulli; fold = SUMA pulsów w `ybins` binach; P3(n) z nachylenia fazy w oknie 20 P; błędy z 4 zakresów długości
(`coherent_fold_jackknife`). Przy P3 ≲ 15 f_c = 1/300 jest za wolne — fold ≈ fold stałym P3.

Typowe ustawienia: domyślne; dla słabych (S/N filtra < 2) lub małego pokrycia P3(n): threshold_q=0.2–0.3,
kernel_sigma=0.07; gdy „nulle” to w rzeczywistości modulacja natężenia: split_nulls=False.

=====================================================================================================================
TESTY, KTÓRE ZROBILIŚMY (szczegóły: docs/coherent_fold_params.md §3–7 i dziennik)
=====================================================================================================================
 * Dobór f_c (H1, 10 pulsarów): optimum 0.07–0.33·f3; stałe 1/300 dawało 2–6× mniejsze ΔR² przy P3 ≲ 15.
   Na danych auto f_c: J1750-3503 1/147, J0034-0721 1/66, J0818-3232 1/35, J1825+0004 1/142; P3(n) w 60–92% pulsów;
   2–10 s na pulsar. Wzrost ΔR² nie zawsze oznacza lepszy fold dryfu (J1946+1805: modulacja natężenia).
 * Maskowanie i adaptacyjna nośna (H2): maska |s| usuwa skoki P3(n) (przejścia |s| przez zero); adaptacyjna nośna
   zmniejsza maskowanie z 26–74% do 10–36% pulsów; pokrycie śladu sLRFS 93–100%.
 * Kalman (H3): najostrzejsze pasma, ale P3(n) bezużyteczne na danych; zostaje wariant C.
 * Benchmark syntetyczny H4 (36 przypadków, P3 6–45, wędrówka/skok, 30% nulli, S/N 1.5–5): korelacja foldu
   z prawdą 0.68 (f_c = 1/300) → 0.85 (auto); faza R 0.93; błąd P3(n) 5.6% (mediana) przy pokryciu 94%.
 * Benchmark H5 (200 losowych przypadków: P3 3–50 stałe/wędrujące/skokowe, nulle, S/N 0.8–8, 15% AM): mediana błędu
   P3(n) 1.5%, pokrycie 70%, faza R 0.87; stałe/wędrówka/skok 0.5/2.4/1.7%; krótkie nulle 3.1%; S/N < 1.5 2.1%;
   AM 3.4%; punkty z błędem > 20%: 4%; skok P3 wykryty 61/78; w ±2σ tylko 77% (błędy niedoszacowane ~1.5×).
 * Okno P3(n) ≥ 30 P: błąd dla P3 < 8 1.6% → 0.85% (bez wpływu na większe P3).
 * Fold jądrowy vs biny (6c): lepszy w 40/40 syntetyków (korelacja 0.813 vs 0.761 dla najlepszej liczby binów
   z CV i 0.677 dla ybins = 2·P3), wyższy LOO na 7/7 pulsarach — dlatego domyślny.
 * Czy pochylenie jest w danych (§5, J0255 i 5 innych): faza z połowy binów + fold drugiej połowy daje ten sam wzór
   (korelacja 0.994); lokalne foldy stałym P3 pokazują to samo przesunięcie fazy.
 * Ślepy przegląd 302 + 227 pulsarów (2026-10-10/11, §6.5): ocena wzrokowa foldu pod losowymi ID; dryf wśród
   dryferów Song+23 spada z Ė (Ė ≤ 2·10³²: 58%, 2·10³²–10³³: 16%, > 10³³: ~6%).
 * Kalibracja klasyfikacji „dryf / brak dryfu” (2026-10-11, syntetyki, 30 dryfów + 30 AM): test SVD (poniżej)
   24/30 wykrytych przy 1/30 fałszywych alarmów. Samego foldu koherentnego jako klasyfikatora NIE skalibrowano.

=====================================================================================================================
SŁABOŚCI (ważne przy interpretacji)
=====================================================================================================================
 1. FOLD NIE JEST TESTEM NA DRYF. Faza każdego pulsu pochodzi z jego własnego podobieństwa do szablonu (krok 2), więc
    uśrednianie pulsów posortowanych wg tej fazy odtwarza szablon — także pochylone pasma — nawet gdy kolejne pulsy
    nie przesuwają się w czasie. Sprawdzone (2026-10-11): fold z PRZETASOWANYCH pulsów z szablonem liczonym od nowa
    daje strukturę o amplitudzie 0.56–1.41 amplitudy danych, z ukośnymi pasmami o losowym kierunku. Kolejność pulsów
    wchodzi tylko przez filtr (krok 3). Przykład: J0152-1637 — fold pokazuje pasma, pojedyncze pulsy i test obrotu
    SVD nie pokazują dryfu. Do rozstrzygania dryf/brak potrzebny test używający kolejności pulsów.
 2. Kontrast foldu zawyżony ~2× (null ze wspólnym szablonem: 25–60% amplitudy foldu) — nie traktować jako głębokości
    modulacji.
 3. Null (fold przetasowanych pulsów) nie nadaje się do klasyfikacji: ze wspólnym szablonem odtwarza szablon razem
    z kierunkiem pochylenia (korelacja z danymi 0.75–0.96) i „zjada” prawdziwy dryf (14/18 „artefaktów” ślepego
    przeglądu miało spójny gradient fazy LRFS); z szablonem liczonym od nowa ma losową fazę i jego średnia się znosi.
 4. Falowanie „S” i zygzak przy P3 ≈ 2–3 to w większości artefakt; przy P3 ≈ 2 alias — kierunek dryfu nieokreślony.
 5. Przegapiony skok P3: gdy skok częstości > f_c, P3(n) zostaje przy starym P3 bez sygnału błędu (syntetyk #156:
    25.8 zamiast 18.9 przez ~400 P) — wynik pewny i fałszywy.
 6. Okno P3(n) zależy od f_c: przy f_c = f3/16 okno i brzegi ~8·P3, przy częstych nullach pokrycie spada do zera;
    piki P3(n) na początkach odcinków (J1137, J1655, J1056).
 7. Auto f_c optymalizuje fold, nie P3(n): przy wędrówce P3 ≳ 20% bywa za wolny (6/56 przypadków).
 8. Błędy σ niedoszacowane ~1.5× (w ±2σ 77% zamiast 95%).
 9. Nulle: skok fazy po nullu przechodzi przez filtr; głęboka AM bywa brana za nulle (syntetyk #20: 40%).
10. Słabe pulsary (S/N ≲ 2): fold koherentny bywa płaski, a fold stałym P3 pokazuje modulację (J0905-6019,
    J1742-4616). Przy silnej modulacji natężenia próg |s| z tasowania zawyżony (J1703-3241: pokrycie 2% przy S/N 11.7).
11. Złe dane (zła efemeryda, RFI, modulowana linia bazowa) dają pozorne pochylenie — zawsze oglądać panel energii;
    RFI zerować i liczyć ponownie.
13. Próg |s| z tasowania zależy od ziarna losowego: tylko 5 tasowań, a |s| po filtrze jest silnie
    autokorelowane (efektywnie kilkadziesiąt niezależnych wartości). J1750-3503: próg 0.26–0.33 dla 20 ziaren
    (Julia: 0.316), pokrycie P3(n) 70–78%. Wynik P3(n) na brzegach odcinków zależy więc od ziarna.
12. Odwrócenia dryfu krótsze niż cykl P3 (J1750-3503, epizody ~28 P przy P3 ≈ 49) nie są rozdzielane — fold uśrednia
    oba kierunki; minimum |s| nie jest wiarygodnym wskaźnikiem odwrócenia.

Uzupełnienie — test dryfu z kolejności pulsów (osobny skrypt, nie w tym pliku): `~/claude/work/scripts/svd_drift2.jl`.
SVD fluktuacji (po usunięciu składowej „cały profil jaśniej/ciemniej”), dwa wiodące wzory w długości, amplitudy
(u1, u2) kolejnych pulsów: dryf = obrót w jedną stronę, AM = brak obrotu; główny wykres — skumulowany obrót na tle
nulli z losowymi fazami Fouriera. Na syntetykach 24/30 dryfów, 1/30 fałszywych alarmów (koniec skumulowanego
obrotu poza nullami, p < 0.05); zgubione głównie przy S/N ≤ 1.4. Nie potrzebuje P3 ani szablonu.
"""

import argparse
import json
import os

import numpy as np
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ===================================================================================================================
# Wejście / wyjście (Data.load_ascii, Data.zap!, Tools.read_params, Tools.p3fold, Tools.intensity_pulses)
# ===================================================================================================================

def read_params(path):
    """params.json → dict."""
    with open(path) as f:
        return json.load(f)


def load_ascii(path):
    """ASCII PSRCHIVE (`pdv -t`): nagłówek, potem wiersze `puls kanał bin natężenie` → macierz pulsy × biny."""
    with open(path) as f:
        head = f.readline().split()
    npulse, nbin = int(head[5]), int(head[11])
    tab = np.loadtxt(path, skiprows=1, usecols=(0, 2, 3))
    data = np.zeros((npulse, nbin))
    data[tab[:, 0].astype(int), tab[:, 1].astype(int)] = tab[:, 2]
    return data


def zap(data, ranges):
    """Zeruje pulsy z `ranges` (napis „a-b c d-e”, numeracja 0-bazowa jak w PSRCHIVE). Działa w miejscu."""
    if not ranges:
        return data
    n = data.shape[0]
    for tok in str(ranges).split():
        b = tok.split("-")
        lo = int(b[0]); hi = lo if len(b) == 1 else int(b[1])
        data[max(lo, 0):min(hi, n - 1) + 1, :] = 0.0
    return data


def p3fold_const(data, p3, ybins=10, start_ind=None):
    """Fold stałym P3 (Tools.p3fold): SUMA pulsów w `ybins` binach fazy n/P3, od pulsu floor(P3/2)."""
    if start_ind is None:
        start_ind = int(np.floor(p3 / 2))
    d = data[start_ind:, :]
    out = np.zeros((ybins, d.shape[1]))
    dp3 = p3 / ybins
    for i in range(d.shape[0]):
        j = int(np.floor(np.mod((i + 1) / dp3, ybins)))
        out[j, :] += d[i, :]
    return out


def intensity_pulses(data):
    """Energia pulsów znormalizowana do [0, 1] (Tools.intensity_pulses)."""
    e = data.sum(axis=1).astype(float)
    if e.min() < 0:
        e = e - e.min()
    return e / e.max()


# ===================================================================================================================
# Wspólne elementy foldu koherentnego
# ===================================================================================================================

def _on(bin_st, bin_end):
    """Zakres on-pulse z 1-bazowych, włącznych bin_st, bin_end (konwencja params.json / Julii)."""
    return slice(bin_st - 1, bin_end)


def template_weights(X, p3):
    """Filtr dopasowany przy f3 = 1/P3: conj binu k = round(N/P3) pełnej FFT po pulsach."""
    n = X.shape[0]
    k = int(np.clip(np.round(n / p3), 1, n // 2))
    return np.conj(np.fft.fft(X, axis=0)[k, :])


def lowpass(x, cutoff, order=6):
    """Zero-fazowy Butterworth (filtfilt) części rzeczywistej i urojonej; cutoff w cyklach/puls (fs = 1)."""
    sos = signal.butter(order, cutoff, btype="low", fs=1.0, output="sos")
    pad = min(6 * sos.shape[0], len(x) - 1)
    f = lambda v: signal.sosfiltfilt(sos, v, padtype="odd", padlen=pad)
    return f(np.real(x)) + 1j * f(np.imag(x))


def build_template(data, bins, ybins):
    """SUMA pulsów w przypisanych binach fazy (0-bazowe) → ybins × N_binów."""
    out = np.zeros((ybins, data.shape[1]))
    np.add.at(out, bins, data)
    return out


def phase_bins(phase, ybins):
    return np.floor(np.mod(phase / (2 * np.pi) * ybins, ybins)).astype(int)


def matched_snr(data, p3, bin_st, bin_end):
    """S/N filtra dopasowanego: |L_on| względem szumu L w binach off-pulse."""
    n, nb = data.shape
    k = int(np.clip(np.round(n / p3), 1, n // 2))
    L = np.fft.fft(data, axis=0)[k, :]
    on = _on(bin_st, bin_end)
    off = np.r_[0:bin_st - 1, bin_end:nb]
    if off.size == 0:
        return np.inf
    s_off = np.std(np.r_[L[off].real, L[off].imag], ddof=1)
    n_on = bin_end - bin_st + 1
    return np.inf if s_off == 0 else np.sqrt(np.sum(np.abs(L[on]) ** 2)) / (s_off * np.sqrt(n_on))


# ===================================================================================================================
# Wersja pierwotna: coherent_fold, coherent_fold_jackknife (p3fold_coherent)
# ===================================================================================================================

def windowed_slope(y, window):
    """Lokalne nachylenie dy/dx z regresji w przesuwnym oknie `window` próbek (zwęża się przy brzegach)."""
    n = len(y); out = np.zeros(n); h = max(1, window // 2)
    for i in range(n):
        lo, hi = max(0, i - h), min(n - 1, i + h)
        if hi - lo < 2:
            out[i] = (y[i] - y[i - 1]) if i > 0 else ((y[1] - y[0]) if n > 1 else 0.0)
            continue
        xs = np.arange(lo, hi + 1); ys = y[lo:hi + 1]
        den = np.sum((xs - xs.mean()) ** 2)
        out[i] = 0.0 if den == 0 else np.sum((xs - xs.mean()) * (ys - ys.mean())) / den
    return out


def coherent_fold(data, p3, bin_st, bin_end, ybins=10, lowpass_cutoff=1 / 200, filter_order=4, p3_window=20):
    """Fold koherentny (wersja pierwotna): stała nośna 2πn/P3, stałe f_c, fold = suma w binach."""
    n = data.shape[0]; on = _on(bin_st, bin_end)
    X = data[:, on] - data[:, on].mean(axis=0)
    z = X @ template_weights(data[:, on], p3)          # szablon z danych bez odjęcia średniej (jak w Julii)
    nn = np.arange(1, n + 1)
    s = lowpass(z * np.exp(-1j * 2 * np.pi / p3 * nn), lowpass_cutoff, filter_order)
    phase = 2 * np.pi / p3 * nn + np.unwrap(np.angle(s))
    b = phase_bins(phase, ybins)
    return dict(folded=build_template(data, b, ybins), phase=phase, bin=b,
                p3_per_pulse=2 * np.pi / windowed_slope(phase, p3_window), snr=matched_snr(data, p3, bin_st, bin_end))


def coherent_fold_jackknife(data, p3, bin_st, bin_end, ybins=10, lowpass_cutoff=1 / 200, filter_order=4,
                            p3_window=20, n_groups=4):
    """coherent_fold + błędy P3(n) i fazy z `n_groups` rozłącznych zakresów długości (std/√n)."""
    main = coherent_fold(data, p3, bin_st, bin_end, ybins, lowpass_cutoff, filter_order, p3_window)
    n = data.shape[0]
    edges = np.round(np.linspace(bin_st, bin_end + 1, n_groups + 1)).astype(int)
    gp3 = np.full((n_groups, n), np.nan); gph = np.full((n_groups, n), np.nan)
    for g in range(n_groups):
        st, en = edges[g], edges[g + 1] - 1
        if en < st:
            continue
        r = coherent_fold(data, p3, st, en, ybins, lowpass_cutoff, filter_order, p3_window)
        gp3[g] = r["p3_per_pulse"]; gph[g] = r["phase"] - r["phase"].mean()
    main["p3_per_pulse_err"] = np.std(gp3, axis=0, ddof=1) / np.sqrt(n_groups)
    main["phase_err"] = np.std(gph, axis=0, ddof=1) / np.sqrt(n_groups)
    return main


# ===================================================================================================================
# Wersja bieżąca: wariant C (coherent_fold_agent) i jej części
# ===================================================================================================================

def carrier_track(X, p3, cutoff, filter_order=6, niter=0, threshold=None):
    """s = LP(z·e^{−i·nośna}); przy niter > 0 nośna przechodzi na fazę z poprzedniego przebiegu
    (rozwiniętą po pulsach z |s| ≥ threshold, liniowo pomiędzy). X — z odjętą średnią kolumn."""
    n = X.shape[0]; nn = np.arange(1, n + 1)
    z = X @ template_weights(X, p3)
    carrier = 2 * np.pi / p3 * nn
    s = lowpass(z * np.exp(-1j * carrier), cutoff, filter_order)
    for _ in range(niter):
        idx = np.arange(n) if threshold is None else np.flatnonzero(np.abs(s) >= threshold)
        if idx.size < 3:
            break
        th = np.unwrap(np.angle(s[idx]))
        carrier = carrier + np.interp(np.arange(n), idx, th)
        s = lowpass(z * np.exp(-1j * carrier), cutoff, filter_order)
    return s, carrier


def shuffle_level(X, p3, cutoff, filter_order=6, niter=0, threshold=None, q=0.5, nshuffle=5, seed=1):
    """Kwantyl q |s| dla danych o przetasowanej kolejności pulsów (brak okresowości) — poziom szumu |s|."""
    rng = np.random.default_rng(seed); vals = []
    for _ in range(nshuffle):
        s, _ = carrier_track(X[rng.permutation(X.shape[0])], p3, cutoff, filter_order, niter, threshold)
        vals.append(np.abs(s))
    return np.quantile(np.concatenate(vals), q)


def _fold_r2(Y, phase, ybins):
    """Ułamek wariancji fluktuacji wyjaśniony foldem w ybins binach minus obciążenie szumem (ybins−1)/(N−1)."""
    n = Y.shape[0]; b = phase_bins(phase, ybins); between = 0.0
    for y in range(ybins):
        m = b == y
        if m.any():
            between += m.sum() * np.sum(Y[m].mean(axis=0) ** 2)
    return between / np.sum(Y ** 2) - (ybins - 1) / (n - 1)


def _bin_halves(nb):
    """Przeplatane bloki po 4 biny: dwie połówki o niezależnym szumie, obie pokrywają cały profil."""
    A = np.array([j for j in range(nb) if (j // 4) % 2 == 0]); return A, np.setdiff1d(np.arange(nb), A)


def _cv_r2(X, A, B, p3, cutoff, ybins, filter_order):
    sc = 0.0
    for P, Q in ((A, B), (B, A)):
        s, carrier = carrier_track(X[:, P], p3, cutoff, filter_order)
        sc += _fold_r2(X[:, Q], carrier + np.angle(s), ybins) / 2
    return sc


def auto_cutoff(data, p3, bin_st, bin_end, ybins=10, grid=(1 / 16, 1 / 10, 1 / 8, 1 / 6, 1 / 4, 1 / 3),
                filter_order=6, nshuffle=5, seed=2):
    """f_c = argmax ΔR² (CV po długości − to samo dla przetasowanych pulsów) na siatce g·f3."""
    X = data[:, _on(bin_st, bin_end)]; X = X - X.mean(axis=0); n = X.shape[0]
    A, B = _bin_halves(X.shape[1]); cut = [g / p3 for g in grid]
    rng = np.random.default_rng(seed); perms = [rng.permutation(n) for _ in range(nshuffle)]
    score = [_cv_r2(X, A, B, p3, fc, ybins, filter_order) -
             np.mean([_cv_r2(X[pp], A, B, p3, fc, ybins, filter_order) for pp in perms]) for fc in cut]
    return dict(cutoff=cut[int(np.argmax(score))], grid=cut, score=score)


def weighted_p3(s, carrier, keep, window, bounds=(2.0, np.inf)):
    """P3(n) z nachylenia fazy (ważone |s|²) w ±window/2 pulsów, osobno w każdym ciągłym odcinku `keep`."""
    n = len(s); out = np.full(n, np.nan); wt = np.abs(s) ** 2; h = max(2, window // 2); i = 0
    while i < n:
        if not keep[i]:
            i += 1; continue
        j = i
        while j < n - 1 and keep[j + 1]:
            j += 1
        run = np.arange(i, j + 1); L = len(run)
        ph = carrier[run] + np.unwrap(np.angle(s[run]))
        for m in range(L):
            lo, hi = max(0, m - h), min(L - 1, m + h)
            if hi - lo + 1 < max(5, h // 2):
                continue
            x = run[lo:hi + 1].astype(float); y = ph[lo:hi + 1]; w = wt[run[lo:hi + 1]]
            xm = np.sum(w * x) / np.sum(w); ym = np.sum(w * y) / np.sum(w)
            v = 2 * np.pi / (np.sum(w * (x - xm) * (y - ym)) / np.sum(w * (x - xm) ** 2))
            out[run[m]] = v if bounds[0] <= v <= bounds[1] else np.nan
        i = j + 1
    return out


def energy_nulls(data, on, smooth=5, minlen=2):
    """Nulle z energii: średnia po `smooth` P, ułamek nulli 2·frac(Ē<0); < 5% → brak; epizody ≥ `minlen` P."""
    n = data.shape[0]; E = data[:, on].sum(axis=1); h = smooth // 2
    Es = np.array([E[max(0, i - h):min(n, i + h + 1)].mean() for i in range(n)])
    nf = 2 * np.mean(Es < 0); out = np.zeros(n, bool)
    if nf < 0.05:
        return out
    m = Es < np.quantile(Es, min(nf, 0.9)); i = 0
    while i < n:
        if m[i]:
            j = i
            while j < n - 1 and m[j + 1]:
                j += 1
            if j - i + 1 >= minlen:
                out[i:j + 1] = True
            i = j + 1
        else:
            i += 1
    return out


def variant_c(X, p3, cutoff, filter_order=6, niter=2, nshuffle=5, threshold_q=0.5, exclude=None):
    """Wariant C na jednym zakresie: adaptacyjna nośna, próg z tasowania, bez brzegów i `exclude` (nulli)."""
    n = X.shape[0]
    thr0 = shuffle_level(X, p3, cutoff, filter_order, nshuffle=nshuffle, q=threshold_q)
    s, carrier = carrier_track(X, p3, cutoff, filter_order, niter, thr0)
    thr = thr0 if niter == 0 else shuffle_level(X, p3, cutoff, filter_order, niter, thr0, threshold_q, nshuffle)
    nedge = max(4, int(np.round(1 / (2 * cutoff))))
    window = max(nedge, int(np.round(3 * p3)), 30)
    edge = np.zeros(n, bool); edge[:min(n, nedge)] = True; edge[max(0, n - nedge):] = True
    keep = (np.abs(s) >= thr) & ~edge
    if exclude is not None:
        keep &= ~exclude
    return dict(s=s, carrier=carrier, keep=keep, threshold=thr, window=window,
                p3=weighted_p3(s, carrier, keep, window, bounds=(2.0, 3 * p3)))


def _phase_grid(nphase):
    return np.arange(nphase) / nphase


def _kernel_weights(g, u, sigma):
    """Zawinięty Gauss między fazami siatki g i pulsów u (cykle): G × N."""
    d = g[:, None] - u[None, :]; d -= np.round(d)
    return np.exp(-d ** 2 / (2 * sigma ** 2))


def kernel_fold(data, phase, sigma, nphase=64):
    """Fold jądrowy: dla `nphase` faz średnia wszystkich pulsów ważona zawiniętym Gaussem (σ w cyklach)."""
    W = _kernel_weights(_phase_grid(nphase), np.mod(phase / (2 * np.pi), 1.0), sigma)
    return (W @ data) / W.sum(axis=1, keepdims=True)


def _loo_kernel_score(Y, u, sigma, ngrid=256):
    """LOO foldu jądrowego na siatce `ngrid` faz: 1 − Σ|Y − Ŷ|² / Σ|Y|²."""
    W = _kernel_weights(_phase_grid(ngrid), u, sigma); S = W @ Y; C = W.sum(axis=1)
    k = np.mod(np.round(u * ngrid).astype(int), ngrid); err = 0.0
    for i in range(len(u)):
        c = C[k[i]] - 1.0
        err += np.sum(Y[i] ** 2) if c <= 1e-9 else np.sum(((S[k[i]] - Y[i]) / c - Y[i]) ** 2)
    return 1 - err / np.sum(Y ** 2)


def kernel_sigma(data, p3, bin_st, bin_end, cutoff, grid=(0.01, 0.015, 0.02, 0.03, 0.045, 0.07, 0.1, 0.15, 0.2, 0.3),
                 filter_order=6, niter=2, nshuffle=5, threshold_q=0.5, exclude=None):
    """σ foldu jądrowego z LOO na połówkach binów długości; największe σ w 1% od maksimum."""
    X = data[:, _on(bin_st, bin_end)]; X = X - X.mean(axis=0)
    A, B = _bin_halves(X.shape[1]); score = np.zeros(len(grid))
    for P, Q in ((A, B), (B, A)):
        c = variant_c(X[:, P], p3, cutoff, filter_order, niter, nshuffle, threshold_q, exclude)
        u = np.mod((c["carrier"] + np.angle(c["s"])) / (2 * np.pi), 1.0)
        Y = X[:, Q] - X[:, Q].mean(axis=0)
        score += np.array([_loo_kernel_score(Y, u, sg) for sg in grid]) / 2
    best = score.max(); i = np.flatnonzero(score >= best - 0.01 * abs(best))[-1]
    return dict(sigma=grid[i], grid=list(grid), score=score)


def coherent_fold_agent(data, p3, bin_st, bin_end, ybins=10, lowpass_cutoff="auto", filter_order=6, niter=2,
                        nshuffle=5, n_groups=4, threshold_q=0.5, split_nulls=True, fold="kernel",
                        kernel_sigma_="auto", nphase=64):
    """Fold koherentny, wariant C (opis w nagłówku pliku). Zwraca dict:
    folded (wg `fold`), folded_kernel, kernel_sigma, folded_bins (średnia w binie), counts, phase [rad], bin,
    p3_per_pulse (NaN gdzie nie zmierzone), p3_per_pulse_err, used, nulls, amplitude |s|, threshold,
    lowpass_cutoff, cutoff_grid, cutoff_score, sigma_grid, sigma_score, snr."""
    if fold not in ("kernel", "bins"):
        raise ValueError("fold musi być 'kernel' albo 'bins'")
    n = data.shape[0]; on = _on(bin_st, bin_end)
    if lowpass_cutoff == "auto":
        ac = auto_cutoff(data, p3, bin_st, bin_end, ybins, filter_order=filter_order, nshuffle=nshuffle)
        fc, cgrid, cscore = ac["cutoff"], ac["grid"], ac["score"]
    else:
        fc, cgrid, cscore = float(lowpass_cutoff), None, None
    X = data[:, on] - data[:, on].mean(axis=0)
    nulls = energy_nulls(data, on) if split_nulls else np.zeros(n, bool)
    c = variant_c(X, p3, fc, filter_order, niter, nshuffle, threshold_q, nulls)
    phase = c["carrier"] + np.angle(c["s"])
    b = phase_bins(phase, ybins)
    counts = np.bincount(b, minlength=ybins)
    with np.errstate(invalid="ignore", divide="ignore"):
        folded_bins = build_template(data, b, ybins) / counts[:, None]
    folded_bins[counts == 0] = np.nan
    if kernel_sigma_ == "auto":
        ks = kernel_sigma(data, p3, bin_st, bin_end, fc, filter_order=filter_order, niter=niter, nshuffle=nshuffle,
                          threshold_q=threshold_q, exclude=nulls)
        ksig, sgrid, sscore = ks["sigma"], ks["grid"], ks["score"]
    else:
        ksig, sgrid, sscore = float(kernel_sigma_), None, None
    folded_kernel = kernel_fold(data, phase, ksig, nphase)
    edges = np.round(np.linspace(bin_st, bin_end + 1, n_groups + 1)).astype(int)
    gp3 = np.full((n_groups, n), np.nan)
    for g in range(n_groups):
        st, en = edges[g], edges[g + 1] - 1
        if en - st < 2:
            continue
        Xg = data[:, st - 1:en]; Xg = Xg - Xg.mean(axis=0)
        gp3[g] = variant_c(Xg, p3, fc, filter_order, niter, nshuffle, threshold_q, nulls)["p3"]
    err = np.full(n, np.nan)
    for i in range(n):
        v = gp3[:, i][np.isfinite(gp3[:, i])]
        if len(v) >= 2:
            err[i] = np.std(v, ddof=1) / np.sqrt(len(v))
    return dict(folded=folded_kernel if fold == "kernel" else folded_bins, fold=fold, folded_kernel=folded_kernel,
                kernel_sigma=ksig, sigma_grid=sgrid, sigma_score=sscore, folded_bins=folded_bins, counts=counts,
                phase=phase, bin=b, p3_per_pulse=c["p3"], p3_per_pulse_err=err, used=c["keep"], nulls=nulls,
                amplitude=np.abs(c["s"]), threshold=c["threshold"], lowpass_cutoff=fc, cutoff_grid=cgrid,
                cutoff_score=cscore, snr=matched_snr(data, p3, bin_st, bin_end))


# ===================================================================================================================
# Wykres (Plot.p3fold_compare) i funkcje „całościowe” (SpaTs.p3fold_coherent[_agent])
# ===================================================================================================================

def p3fold_compare(folded, folded_const, p3_per_pulse, p3_nominal, outdir, bin_st=None, bin_end=None, cmap="viridis",
                   darkness=1.0, repeat_num=4, name_mod="0", label="coherent fold", p3_per_pulse_err=None,
                   intensity=None, p3_ylim=None, figtitle=None):
    """Fold koherentny | fold stałym P3 (powtórzone `repeat_num` cykli) | energia pulsów | P3(n) z błędami."""
    nb = folded.shape[1]; bin_st = 1 if bin_st is None else bin_st; bin_end = nb if bin_end is None else bin_end
    sl = _on(bin_st, bin_end)
    dv = np.tile(folded[:, sl], (repeat_num, 1)); dc = np.tile(folded_const[:, sl], (repeat_num, 1))
    le = dv.shape[0]; ticks = [le // 4, le // 2, le * 3 // 4]
    tl = [f"{repeat_num * f}$P_3$" for f in (0.25, 0.5, 0.75)]
    plt.rc("font", size=8.); plt.rc("axes", linewidth=0.5); plt.rc("lines", linewidth=0.5)
    nrows, img, p3r, h = (3, 2, 1, 6.29921) if intensity is None else (8, 4, 3, 7.3)
    fig = plt.figure(figsize=(6.29921, h))
    plt.subplots_adjust(left=0.1, bottom=0.07, right=0.99, top=0.96 if figtitle is None else 0.93,
                        wspace=0.15, hspace=0.35)
    if figtitle is not None:
        fig.suptitle(figtitle, fontsize=11, y=0.985)
    ax = plt.subplot2grid((nrows, 2), (0, 0), rowspan=img)
    ax.imshow(dv, origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness * np.nanmax(dv))
    ax.set_yticks(ticks); ax.set_yticklabels(tl); ax.set_title(label)
    ax = plt.subplot2grid((nrows, 2), (0, 1), rowspan=img)
    ax.imshow(dc, origin="lower", cmap=cmap, interpolation="none", aspect="auto", vmax=darkness * np.nanmax(dc))
    ax.tick_params(labelleft=False); ax.set_title("constant $P_3$")
    x = np.arange(1, len(p3_per_pulse) + 1)
    if intensity is not None:
        ax = plt.subplot2grid((nrows, 2), (img, 0), colspan=2)
        ax.minorticks_on(); ax.plot(x, intensity, c="grey", lw=0.6); ax.set_xlim(1, len(x))
        ax.tick_params(labelbottom=False); ax.set_ylabel("intensity")
    ax = plt.subplot2grid((nrows, 2), (nrows - p3r, 0), rowspan=p3r, colspan=2)
    ax.minorticks_on()
    if p3_per_pulse_err is not None:
        ax.fill_between(x, p3_per_pulse - p3_per_pulse_err, p3_per_pulse + p3_per_pulse_err, color="grey",
                        alpha=0.3, linewidth=0)
    ax.plot(x, p3_per_pulse, c="grey", lw=0.8); ax.axhline(p3_nominal, c="red", ls="--", lw=0.8)
    ax.set_xlim(1, len(x))
    if p3_ylim is not None:
        ax.set_ylim(p3_ylim)
    ax.set_xlabel("pulse number"); ax.set_ylabel("$P_3$ (P)")
    for ext in ("pdf", "png"):
        path = os.path.join(outdir, f"{name_mod}_p3fold_compare.{ext}"); fig.savefig(path); print(path)
    plt.close(fig)


def _load(outdir, datafile):
    p = read_params(os.path.join(outdir, "params.json"))
    data = load_ascii(os.path.join(outdir, datafile))
    zap(data, p.get("zaps"))
    return p, data


def p3fold_coherent(outdir, ybins=None, lowpass_cutoff=1 / 300, filter_order=6, n_groups=4, darkness=1.0,
                    datafile="pulsar.debase.txt", plotdir=None, name_mod="pulsar_coherent", figtitle=None):
    """Wersja pierwotna (SpaTs.p3fold_coherent): stałe f_c, fold = suma w binach, błędy z 4 zakresów długości."""
    plotdir = outdir if plotdir is None else plotdir
    p, data = _load(outdir, datafile)
    yb = int(np.round(p["p3_ybins"])) if ybins is None else ybins
    p3 = float(p["p3"]); bst, ben = int(p["bin_st"]), int(p["bin_end"])
    r = coherent_fold_jackknife(data, p3, bst, ben, yb, lowpass_cutoff, filter_order, n_groups=n_groups)
    print(f"Matched-filter SNR: {r['snr']:.1f}")
    p3fold_compare(r["folded"], p3fold_const(data, p3, yb), r["p3_per_pulse"], p3, plotdir, bin_st=bst, bin_end=ben,
                   darkness=darkness, name_mod=name_mod, label="coherent fold",
                   p3_per_pulse_err=r["p3_per_pulse_err"], intensity=intensity_pulses(data[:, _on(bst, ben)]),
                   figtitle=figtitle)
    return r


def p3fold_coherent_agent(outdir, ybins=None, lowpass_cutoff="auto", filter_order=6, n_groups=4, threshold_q=0.5,
                          split_nulls=True, fold="kernel", kernel_sigma_="auto", nphase=64, darkness=1.0,
                          datafile="pulsar.debase.txt", plotdir=None, name_mod="pulsar_coherent_agent",
                          figtitle=None):
    """Wersja bieżąca (SpaTs.p3fold_coherent_agent): wariant C + fold jądrowy; wykres jak p3fold_compare."""
    plotdir = outdir if plotdir is None else plotdir
    p, data = _load(outdir, datafile)
    yb = int(np.round(p["p3_ybins"])) if ybins is None else ybins
    p3 = float(p["p3"]); bst, ben = int(p["bin_st"]), int(p["bin_end"])
    r = coherent_fold_agent(data, p3, bst, ben, yb, lowpass_cutoff, filter_order, n_groups=n_groups,
                            threshold_q=threshold_q, split_nulls=split_nulls, fold=fold, kernel_sigma_=kernel_sigma_,
                            nphase=nphase)
    fc = r["lowpass_cutoff"]
    print(f"Matched-filter SNR: {r['snr']:.1f}, lowpass_cutoff = 1/{1 / fc:.1f} ({fc * p3:.3f} f3), "
          f"kernel σ = {r['kernel_sigma']:.3g} cycle, nulls (energy) {100 * r['nulls'].mean():.1f}%, "
          f"P3(n) measured in {100 * np.isfinite(r['p3_per_pulse']).mean():.1f}% of pulses")
    label = (f"coherent fold, kernel σ = {r['kernel_sigma']:.3g} ($f_c$ = 1/{int(round(1 / fc))})" if fold == "kernel"
             else f"coherent fold, mean, {yb} bins ($f_c$ = 1/{int(round(1 / fc))})")
    p3fold_compare(np.nan_to_num(r["folded"]), p3fold_const(data, p3, yb), r["p3_per_pulse"], p3, plotdir,
                   bin_st=bst, bin_end=ben, darkness=darkness, name_mod=name_mod, label=label,
                   p3_per_pulse_err=r["p3_per_pulse_err"], intensity=intensity_pulses(data[:, _on(bst, ben)]),
                   figtitle=figtitle)
    return r


def main():
    ap = argparse.ArgumentParser(description="Fold koherentny P3 (port spats: p3fold_coherent[_agent]). "
                                             "Opis metody, testy i słabości w nagłówku pliku.")
    ap.add_argument("outdir", help="katalog z params.json i plikiem pojedynczych pulsów")
    ap.add_argument("--datafile", default="pulsar.debase.txt")
    ap.add_argument("--plotdir", default=None, help="gdzie zapisać wykres (domyślnie outdir)")
    ap.add_argument("--name", default=None, help="przedrostek nazwy wykresu")
    ap.add_argument("--legacy", action="store_true", help="wersja pierwotna p3fold_coherent (f_c = 1/300)")
    ap.add_argument("--cutoff", type=float, default=None, help="stałe f_c [cykle/P] zamiast auto")
    ap.add_argument("--threshold-q", type=float, default=0.5, help="kwantyl progu |s| (słabe: 0.2–0.3)")
    ap.add_argument("--kernel-sigma", type=float, default=None, help="stałe σ foldu jądrowego [cykle] (słabe: 0.07)")
    ap.add_argument("--no-split-nulls", action="store_true", help="nie wykrywać nulli z energii")
    ap.add_argument("--fold", choices=("kernel", "bins"), default="kernel")
    a = ap.parse_args()
    name = a.name or os.path.basename(os.path.normpath(a.outdir))
    if a.legacy:
        p3fold_coherent(a.outdir, lowpass_cutoff=a.cutoff or 1 / 300, datafile=a.datafile, plotdir=a.plotdir,
                        name_mod=name + "_coherent", figtitle=name)
    else:
        p3fold_coherent_agent(a.outdir, lowpass_cutoff=a.cutoff if a.cutoff else "auto", threshold_q=a.threshold_q,
                              split_nulls=not a.no_split_nulls, fold=a.fold,
                              kernel_sigma_=a.kernel_sigma if a.kernel_sigma else "auto", datafile=a.datafile,
                              plotdir=a.plotdir, name_mod=name + "_coherent_agent", figtitle=name)


if __name__ == "__main__":
    main()
