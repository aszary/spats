# Metody wyznaczania pozycji składowych profilu bez dopasowywania Gaussów

Przegląd literatury (lipiec 2026) pod kątem zastąpienia fitów Gaussów w `analyse_p3folds4`
(offset między p3foldami 1023 vs 1523 MHz, per wiersz p3foldu) metodą nieparametryczną.

## 1. Template matching w dziedzinie Fouriera — FFTFIT (Taylor 1992)

Standard całej społeczności timingowej od ponad 30 lat. Przesunięcie między dwoma profilami
mierzy się jako **liniowy gradient fazy** między ich transformatami Fouriera: jeśli
`f(φ) = g(φ − Δφ)`, to `F_k = G_k · e^{-ikΔφ}`, więc Δφ wyznacza się fitem prostej do różnicy
faz harmonik (χ² z wagami z amplitud). Zalety kluczowe:

- **żadnych założeń o kształcie** — "templatem" może być po prostu profil z niższej
  częstotliwości, a mierzy się, o ile trzeba go przesunąć, żeby pasował do wyższej;
- precyzja znacznie poniżej jednego binu, z **analitycznym wzorem na błąd**
  (Taylor 1992, App. A) — nie potrzeba błędów z kowariancji fitu Gaussów;
- implementacje istnieją w PSRCHIVE/PulsePortraiture, algorytm to ~30 linii kodu w Julii
  (FFT + ważona regresja fazy).

Jedyny istotny caveat, dobrze opisany w [Hassall et al. 2012 (LOFAR)](https://www.aanda.org/articles/aa/pdf/2012/07/aa18970-12.pdf):
metoda zakłada, że kształt jest identyczny z dokładnością do skali i przesunięcia. Gdy kształt
ewoluuje z częstotliwością (komponenty przesuwają się różnie), globalny fit daje "średni"
offset obciążony przez zmiany kształtu. Rozwiązanie: **stosować FFTFIT w oknach wokół
poszczególnych składowych** — wtedy dostaje się offset per komponent, dokładnie jak z Gaussów,
tylko bez modelu kształtu.

## 2. Cross-korelacja (CCF) + interpolacja piku

Najstarsza i najczęstsza metoda w analizie pojedynczych pulsów i dryfujących subpulsów — np.
klasyczna analiza dryfu przez CCF w
[Prószyński & Wolszczan 1986 (B0809+74, B1237+25, B1919+21)](https://ui.adsabs.harvard.edu/abs/1986ApJ...307..540P/abstract),
a P2 standardowo liczy się z ACF/CCF w kolejnych longitude'ach
([Weltevrede et al. 2003](https://www.aanda.org/articles/aa/pdf/2003/31/aah4435.pdf)).
Procedura: policz dyskretną CCF między profilem low i high (w oknie komponentu), znajdź
maksimum, a pozycję sub-binową uzyskaj przez interpolację **parabolą** lub — jak w wariancie
"Gaussian interpolation shift" używanym w timingu
([Wang et al. 2024](https://arxiv.org/html/2405.08629)) — dopasowaniem krzywej do kilku punktów
CCF wokół piku (to fit do CCF, nie do profilu, więc nie zakłada kształtu składowej).
Błąd: z krzywizny CCF albo bootstrapem (dodawanie realizacji szumu o zmierzonym RMS off-pulse
i powtarzanie pomiaru). Matematycznie CCF i FFTFIT to niemal to samo; FFTFIT jest po prostu
elegantszą, sub-binową wersją z gotowym błędem.

## 3. Punkty odniesienia (fiducial points) — metody momentów i progów

Rodzina metod używana głównie do wyrównywania profili przy pomiarach DM i ewolucji profilu
z częstotliwością ([Ahuja et al. 2007](https://arxiv.org/pdf/astro-ph/0702440)):

- **pik przez interpolację paraboliczną** 3–5 binów wokół maksimum;
- **centroid ważony natężeniem** (pierwszy moment) w oknie komponentu — bardzo odporny na szum,
  ale czuły na asymetrię i wybór okna;
- **mediana strumienia skumulowanego** (punkt 50% energii w oknie) — jeszcze odporniejszy na
  szum niż pik;
- **punkty przecięcia na x% maksimum** na zboczach składowej i ich środek — Kramer i in.
  używali środka profilu na poziomie 10% mocy jako punktu odniesienia; dla profili
  dwuskładnikowych klasycznie brano środek między zewnętrznymi pikami.

Te metody są trywialne w implementacji i bezmodelowe, ale offset per komponent wymaga dobrze
dobranego okna, a asymetryczna ewolucja kształtu wewnątrz okna przesuwa centroid inaczej niż pik.

## 4. Wygładzanie nieparametryczne + zera pochodnej

Zamiast modelować składowe, wygładza się profil metodą nie zakładającą kształtu i szuka
ekstremów po znaku pochodnej:

- **odszumianie falkowe** — dokładnie to robi `psrsmooth` w PSRCHIVE (falki Daubechies,
  progowanie współczynników) do budowy templatów
  ([Pulsar Timing Techniques, Lorimer](https://ar5iv.labs.arxiv.org/html/1309.1767));
- **regresja procesem Gaussa (GP)** — pozycję piku i jej błąd dostaje się z próbek posteriora;
  coraz częstsza w nowszych pracach;
- **B-splajny + PCA** — tak buduje szerokopasmowe templaty
  [PulsePortraiture (Pennucci)](https://github.com/pennucci/PulsePortraiture/blob/master/ppspline.py):
  PCA na profilach per kanał, wygładzenie falkowe eigenprofili, splajn w częstotliwości
  ([Pennucci 2014, "Elementary Wideband Timing"](https://iopscience.iop.org/article/10.1088/0004-637X/790/2/93)).
  To najbliższy literaturowy odpowiednik problemu "jak profil zmienia się między
  częstotliwościami bez zakładania Gaussów".

## 5. Kontekst A/R (pozycje składowych względem rdzenia)

W pracach o aberracji/retardacji
([Gangadhara & Gupta 2001, B0329+54](https://iopscience.iop.org/article/10.1086/321439);
[Krzeszowski et al. 2009](https://academic.oup.com/mnras/article/393/4/1617/1008138))
słabe składowe wykrywano techniką **window-threshold** (uśrednianie tylko pulsów
przekraczających próg w danym oknie longitude), ale finalne pozycje pików i tak brano
z Gaussów — więc ta gałąź literatury nie pomaga, poza samym window-threshold do wzmacniania
słabych składowych.

## 6. Specyficznie dla p3foldów

W pracach o dryfie na dwóch częstotliwościach (np.
[Bhattacharyya et al. 2008, B0826−34](https://academic.oup.com/mnras/article/383/4/1538/1747098),
[Maan 2019 o sygnaturach karuzeli w wielu częstotliwościach](https://arxiv.org/abs/1812.01010))
offset drift-bandów między częstotliwościami ogląda się właśnie po zfoldowaniu modulo P3,
a mierzy korelacyjnie. Dwie opcje wykraczające poza analizę per wiersz:

- **2D cross-korelacja całych map p3fold** — daje jednocześnie offset w longitude i w fazie P3
  jednym pomiarem, uśredniając po całym cyklu;
- **subpulse phase track** (Weltevrede, PSRSALSA): faza pierwszej harmoniki P3 jako funkcja
  longitude — porównanie phase tracków na dwóch częstotliwościach mierzy offset niezależnie
  od kształtu składowych.

## Rekomendacja

Dla `analyse_p3folds4` najlepszy stosunek wiarygodności do wysiłku ma **FFTFIT per okno
komponentu, per wiersz p3foldu**: okna wokół składowych definiuje się raz (na profilu średnim),
w każdym oknie mierzy się przesunięcie low→high z gradientu fazy w Fourierze, błąd z wzoru
Taylora albo bootstrapem z szumu off-pulse. Znika cały problem zbieżności fitów, zamiany
komponentów miejscami i interaktywnego odrzucania punktów — a "longitude" punktu bierze się
jako centroid okna zamiast średniej μ. Jako sanity check warto dorzucić 2D CCF całych map —
jeden globalny offset do porównania ze średnią ważoną per komponent.

## Źródła

- [Taylor 1992 / przegląd technik timingu (Lorimer)](https://ar5iv.labs.arxiv.org/html/1309.1767)
- [Hassall et al. 2012](https://www.aanda.org/articles/aa/pdf/2012/07/aa18970-12.pdf)
- [Pennucci 2014](https://iopscience.iop.org/article/10.1088/0004-637X/790/2/93)
- [PulsePortraiture / ppspline](https://github.com/pennucci/PulsePortraiture/blob/master/ppspline.py)
- [Wang et al. 2024](https://arxiv.org/html/2405.08629)
- [Ahuja et al. 2007](https://arxiv.org/pdf/astro-ph/0702440)
- [Gangadhara & Gupta 2001](https://iopscience.iop.org/article/10.1086/321439)
- [Krzeszowski et al. 2009](https://academic.oup.com/mnras/article/393/4/1617/1008138)
- [Prószyński & Wolszczan 1986](https://ui.adsabs.harvard.edu/abs/1986ApJ...307..540P/abstract)
- [Weltevrede et al. 2003](https://www.aanda.org/articles/aa/pdf/2003/31/aah4435.pdf)
- [Bhattacharyya et al. 2008 (B0826−34)](https://academic.oup.com/mnras/article/383/4/1538/1747098)
- [Maan 2019](https://arxiv.org/abs/1812.01010)
