# Review protocol (PRISMA 2020)

A protocol is written **before** the searches and not changed afterwards; anything decided
later goes in a dated amendment at the end. That is what makes the review reproducible, and
it is the difference between a systematic review and a bibliography.

Draft of 2026-10-01. Not yet registered.

## 1. The two papers this feeds

| | Paper A — software | Paper B — method |
|---|---|---|
| What it is | a description of the code, its validation and its availability | the measurements and what they mean for the method |
| Candidate venues | SoftwareX, Journal of Open Source Software, Computers & Geosciences | Canadian Geotechnical Journal, Acta Geotechnica, Computers and Geotechnics, Géotechnique Letters |
| Needs the review? | a short related-work section only | **yes**, as the state of the art |
| Ready when | the repository has a DOI and a tagged release | the review and the remaining measurements are done |
| Risk | low; mostly written already | the contributions have to be positioned against what others have done |

They share the validation data, the figures and the measurements, so the review is done once
and serves both. Paper A can go out first and be cited by Paper B.

## 2. Review questions

**RQ1 (broad).** How have image-based full-field deformation measurements been used in
geotechnical physical modelling with large displacements, and what do those methods do about
finite strain, rigid rotation, boundary effects and measurement uncertainty?

**RQ2 (narrow, citation review).** What have the works citing Pinyol & Alvarado (2017) done
with the numerical-particle method: applied it, extended it, criticised it, or replaced it —
and which of their contributions could be folded into this code?

RQ2 is a subset of RQ1 in subject but not in search strategy: one is a database search, the
other a citation search. PRISMA 2020 reports those as separate streams feeding one flow
diagram, which is exactly how they will be reported.

**The anchor records for RQ2.** Three, because the method this code implements and the
moisture measurement it was joined to were published separately, and a work extending either
may cite only one of them:

| | Record | What it anchors |
|---|---|---|
| **A1** | Pinyol, N.M. & Alvarado, M. (2017). *Novel analysis for large strains based on particle image velocimetry*. Canadian Geotechnical Journal 54(7): 933-944. doi:[10.1139/cgj-2016-0327](https://doi.org/10.1139/cgj-2016-0327) | the numerical-particle method this code implements |
| **A2** | Parera, F., Pinyol, N.M. & Alonso, E.E. (2021). *Massive, continuous, and non-invasive surface measurement of degree of saturation by shortwave infrared images*. Canadian Geotechnical Journal 58(6): 749-762. doi:[10.1139/cgj-2019-0051](https://doi.org/10.1139/cgj-2019-0051) | saturation from SWIR images, which `MOISTER=2` implements |
| **A3** | Morales, G., Pinyol, N.M., Alvarado, M. & Alonso, E.E. (2025). *Digital image-based measurement of degree of saturation on moving soil*. Canadian Geotechnical Journal 62: 1-7. doi:[10.1139/cgj-2023-0760](https://doi.org/10.1139/cgj-2023-0760) | the two joined, and the test `examples/dam-break-swir` reproduces |

A3 being recent matters for the reading of the counts: it will have few citations yet, and a
low number there is a fact about its age, not about its uptake. That has to be said in the
results rather than left for a reviewer to point out.

## 3. Sources, and what each is for

| Source | Role in PRISMA | Why |
|---|---|---|
| **Scopus** | database | main source; best citation tracking in engineering |
| **Web of Science** | database | overlaps Scopus but not completely; both are reported |
| **Google Scholar** | database, first 200 hits only | catches theses, reports and conference papers the other two miss; too noisy to screen in full, so the cap is declared in advance |
| **Citation search** on the 2017 paper | *other methods* | the whole of RQ2 |
| **Reference lists** of the included records | *other methods* | backward chaining, one generation |
| **Perplexity / LLM search** | *other methods*, scoping only | see below |

### On Perplexity

Worth using, but **not as a database**, and this has to be declared honestly or the review is
not reproducible. A language-model search does not return the same results twice, gives no
reproducible hit count, and cannot be re-run by a reviewer — three things PRISMA requires of
a database search.

Where it is legitimately useful:

* **before the protocol is fixed**, to find the vocabulary the field actually uses, which
  makes the real search strings better. Terms it surfaces go into the strings below, and that
  is a contribution to the method, not to the results.
* **as a supplementary source** under "other methods", on the same footing as asking a
  colleague. Anything it finds must then be located in Scopus or WoS and cited from there; if
  it cannot be found in a real database, it does not go in.
* **never** for the counts in the flow diagram, and never as the only place a record came
  from.

The same applies to any other AI search tool, this one included.

## 4. Search strings

To be run on title, abstract and keywords, with **no date limit** at the lower end and
records up to the search date. Run them, export the full record set including abstracts, and
keep the result count and the exact date of each run — the flow diagram needs both.

### S1 — method × domain (RQ1 core)

**Scopus**

```
TITLE-ABS-KEY (
  ( "particle image velocimetry" OR PIV OR "digital image correlation" OR DIC
    OR "image-based deformation" OR "image based deformation" OR "optical flow"
    OR photogrammetr* )
  AND
  ( geotechn* OR "soil mechanics" OR slope OR landslide OR embankment OR "retaining wall"
    OR "physical model*" OR centrifuge OR "dam break" OR tailings OR "granular material"
    OR sand OR clay )
)
```

**Web of Science**

```
TS=(
  ("particle image velocimetry" OR PIV OR "digital image correlation" OR DIC
   OR "image-based deformation" OR "optical flow" OR photogrammetr*)
  AND
  (geotechn* OR "soil mechanics" OR slope OR landslide OR embankment OR "retaining wall"
   OR "physical model*" OR centrifuge OR "dam break" OR tailings OR "granular material"
   OR sand OR clay)
)
```

### S2 — S1 narrowed to large deformation (the subset Paper B argues in)

Add to S1:

```
AND TITLE-ABS-KEY (
  "large displacement*" OR "large deformation*" OR "large strain*" OR "finite strain"
  OR Lagrangian OR "material point method" OR "shear band*" OR "post-failure"
  OR "runout" OR "progressive failure"
)
```

### S3 — water content from images (the moisture half of the method)

```
TITLE-ABS-KEY (
  ( "water content" OR moisture OR "degree of saturation" OR suction )
  AND
  ( "image analysis" OR "digital image" OR photograph* OR spectral OR "near infrared"
    OR SWIR OR "short-wave infrared" OR hyperspectral OR colorimetr* )
  AND
  ( soil OR geotechn* OR sand OR silt OR clay OR "porous media" )
)
```

### C1 — citation search (the whole of RQ2)

Not a string but a navigation, repeated for **each of the three anchors** on **each of the
three sources**: nine exports in all.

| Anchor | DOI |
|---|---|
| A1 | `10.1139/cgj-2016-0327` |
| A2 | `10.1139/cgj-2019-0051` |
| A3 | `10.1139/cgj-2023-0760` |

* **Scopus**: search the DOI, open the record, *Cited by*, export all.
* **Web of Science**: search the DOI, *Times Cited*, export all.
* **Google Scholar**: find the record, *Cited by*, export the first 200.

Record the count from each source **and the date of the search**, because these grow. The
step-by-step for running them is in [`review-step-1.md`](review-step-1.md).

## 5. Inclusion and exclusion

**Include** a record when all of these hold:

1. it measures a displacement or strain field from images, or it uses one that was measured
   that way;
2. the material is soil, rock or a granular medium, in the laboratory, in a centrifuge, or
   in the field;
3. it reports something about **how** the measurement is made or interpreted — algorithm,
   validation, uncertainty, strain measure, boundary handling — and not only what was found
   with it;
4. it is a journal article, conference paper, thesis or technical report with a retrievable
   full text;
5. it is in English or Spanish.

**Exclude** when any of these hold:

1. PIV applied to fluids with no granular or soil phase;
2. medical, biological or industrial DIC with no geotechnical bearing;
3. purely numerical work with no image-based measurement;
4. abstract-only records, posters, or anything whose full text cannot be obtained;
5. a duplicate of a record already included, keeping the journal version over the preprint.

For **RQ2** the citing works are screened differently: everything that cites the 2017 paper
is read, and classified rather than excluded (see the extraction form), because a citation
that merely mentions the method in passing is itself a finding about its uptake.

## 6. Screening

Two stages, both recorded in the flow diagram with reasons for exclusion at the second:

1. **Title and abstract**, against the criteria above.
2. **Full text**, for everything that survives.

The user screens; records that are not clear-cut are set aside and decided together, and the
number of those is reported. Where a second screener is available, a sample of at least 20 %
is screened twice and the agreement is reported.

## 7. What is extracted from each included record

Fields chosen so that they answer the two questions and feed the papers directly:

| Field | Why it is there |
|---|---|
| reference, year, venue, DOI | the obvious |
| measurement technique and software | what the field actually uses, and whether PIVlab dominates |
| test type and scale | centrifuge, 1g model, element test, field |
| typical displacement per image pair, in **pixels** | the resolution question this work ran into; most papers do not report it, and that absence is itself a result |
| strain measure: incremental or finite | the core of Paper B |
| whether rigid rotation is discussed at all | the artifact this work quantified |
| what is done at the boundary of the material | where every filter tried here went wrong |
| uncertainty: reported, quantified, or not mentioned | |
| water content from images: yes/no, and how calibrated | the moisture half |
| **for RQ2 only**: how the 2017 method is used — applied as published, extended, modified, criticised, or cited in passing | the uptake question |
| **for RQ2 only**: any extension that could be folded into this code | the practical payoff |

Extraction goes into one spreadsheet, one row per record, with the screening decision and
reason in the same row so the flow diagram can be generated from it rather than counted by
hand.

## 8. How it will be reported

The PRISMA 2020 flow diagram with two identification streams (databases, and other methods),
the PRISMA 2020 checklist, the full search strings with their dates and hit counts, and the
extraction table as supplementary material.

## 9. What the review has to position

Written down now so the review is read with the right questions in mind, and so that a
contribution is not claimed if the review finds it already exists:

1. a legacy code re-engineered with **byte-for-byte equivalence** demonstrated, as a way of
   migrating scientific software without losing continuity;
2. that summing linear strain increments reports **a volumetric contraction that never
   happened** for a rigid rotation — 1.5 % on a 50° turn — and that this reads as densification
   in soil mechanics;
3. the **kinematic vorticity number** brought from structural geology to separate shear bands
   from rotating blocks in a geotechnical PIV field;
4. the boundary correction compared by **hiding nodes that do have data**, which gives a
   ground truth where the boundary normally has none;
5. the sensitivity of **water content from images**, where the width of the calibration band
   dominates everything else by an order of magnitude, and where 82 % of the values in a
   published analysis turn out to be bounds rather than measurements;
6. the negative result that **spatial filtering of the velocity field costs more than it
   gives**, with all of the bias traced to the material boundary;
7. that in the test analysed the soil moves **0.16 px per image pair while PIV resolves
   0.1 px**, which bounds what any post-processing can achieve.

Items 2, 5 and 6 are partly a critique of the authors' own earlier work. That is a choice to
make deliberately rather than by accident, and it is easier to defend when it is measured.
