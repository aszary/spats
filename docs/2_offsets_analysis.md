# Do pulsar profiles get narrower at higher frequency? (short version)

*2026-10-16. Short, plain-language version of [`1_offsets_analysis.md`](1_offsets_analysis.md),
which has all the details, numbers and remaining figures.*

---

## The question

A pulsar sends out a radio beam. Each time it rotates we see a short pulse, and averaging many
pulses gives a **profile**. Many profiles have two humps (**components**). The distance between
the outer humps is the profile width **W**.

A popular idea, called *radius-to-frequency mapping*, says that higher radio frequencies come from
lower in the magnetosphere, where the beam is narrower. If it is true, **W should be smaller at
higher frequency.**

We tested this on MeerKAT data by comparing each pulsar at **1023 MHz** and at **1523 MHz**.

## How we measured it

![Example](figures/offsets_01_example_profile.png)

1. For each pulsar we took the profile in the low band (blue) and in the high band (orange).
2. We fitted a bell-shaped curve (a Gaussian) to each hump and found its centre.
3. We measured W in both bands and computed the relative change

   **ΔW/W = (W at 1523 MHz − W at 1023 MHz) / W**

   A negative value means the profile is narrower at the higher frequency. In the example above
   ΔW/W = −0.08, i.e. the profile is 8% narrower.

Dividing by W lets us compare wide and narrow profiles fairly.

## Which pulsars we kept

We started with 533 pulsars and ended with **193 good measurements**. Most losses came from
one problem: in 184 pulsars the fitted humps did not match between the two bands.

![Fit failures](figures/offsets_05_mismatch_vs_profile_shift.png)

We checked why. The two bands of the data are lined up almost perfectly (horizontal axis, below 1 bin).
So **the data are fine and the fitting method is what fails**: it places the curves on different
features in the two bands (orange points, above the dashed line). A better fitting method could
probably rescue many of these pulsars.

## Result 1: most profiles do get narrower

![Histogram](figures/offsets_06_frac_histogram.png)

- **124 of 193 pulsars (64%)** are narrower at the higher frequency.
- The chance of getting this by luck is about 1 in 10 000.
- The typical change is small: **about 2.4% narrower** over this frequency range.
- Pulsars differ a lot from each other, and some get wider instead.

## Result 2: it is a real effect, not a mistake in the data

A common error in this kind of work is a wrong correction for the delay of radio waves in space
(the dispersion measure, DM). Such an error would **shift the whole profile** to one side,
with both humps moving the same way.

![Edges](figures/offsets_10_edge_offsets.png)

That is not what we see. The **left hump moves right and the right hump moves left**, by the same
amount on average (the cloud of points leans towards the lower right). Both edges move in towards the
centre, so the profile really is narrowing.

## Result 3: an interesting hint (not proven)

Song et al. (2023) divided these pulsars into two types: **drifters**, whose sub-pulses march
across the profile, and **P3-only** pulsars, which only change brightness.

![Drifters vs P3-only](figures/offsets_12_frac_drift_vs_p3only.png)

- Drifters: **69% narrower**.
- P3-only: **only 38% narrower**.

That looks like a clear difference, but there is a catch. P3-only pulsars also lose energy faster
(higher spin-down power Ė), and Ė itself seems to matter a little. When we account for this,
the difference becomes borderline (about a 1 in 20 chance of being luck). We only have 29 P3-only
pulsars, so **we need more data to decide.**

## Things we checked that were less interesting

- **Dependence on Ė:** a weak trend (faster-spinning-down pulsars narrow less). It was seen in the
  first group of pulsars, but **not confirmed** in the second, independent group.
- **Dependence on period P and magnetic field:** none.
- **Dependence on the drift period P3:** none.
- **Width vs period:** slower pulsars have narrower profiles, as expected. This is a sanity check
  that the measurement makes physical sense.
- **Size of the effect:** expressed as a power law in frequency, the index is only ~0.06, a weak
  dependence, as expected above 1 GHz.
- **Noisy measurements** make the effect look weaker, so the true fraction of narrowing pulsars
  is probably above 64%.

## What next

1. **Rescue the 184 failed pulsars** with a better fitting method. This would give ~300 pulsars and
   should settle the question about Ė.
2. **Get more P3-only pulsars** to test Result 3.
3. Check two doubtful pulsars (J1717−3425, J1819+1305), which are left out for now.
