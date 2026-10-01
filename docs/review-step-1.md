# Step 1 — the citation search

The first thing to run, because it is bounded (a few hundred records at most), it answers
the question that matters most — what others have done with this method — and it needs no
judgement calls while collecting.

Three papers × three sources = **nine exports**. Allow about an hour.

Everything is recorded on the way: the counts are part of the result, not bookkeeping. A
PRISMA flow diagram cannot be reconstructed afterwards from the files alone.

## Before you start

Make a folder for it, outside the repository:

```
Prueba\revision\
├── scopus\
├── wos\
├── scholar\
└── counts.csv
```

And create `counts.csv` with this first line, which you will fill in as you go:

```
date,source,anchor,doi,records_found,file
```

## The three papers

| Anchor | Short name | DOI |
|---|---|---|
| **A1** | numerical particles, 2017 | `10.1139/cgj-2016-0327` |
| **A2** | SWIR saturation, 2021 | `10.1139/cgj-2019-0051` |
| **A3** | the two joined, 2025 | `10.1139/cgj-2023-0760` |

---

## 1. Scopus

For **each** of the three DOIs:

1. In the search box choose **Document search**, set the field to **DOI**, paste the DOI,
   search.
2. Open the record. On the right there is **Cited by _N_ documents**. **Write that N down
   now**, with today's date, into `counts.csv`.
3. Click it. You now have the list of citing documents.
4. Select all (the checkbox at the top selects the whole page — make sure it says *All* and
   not just the 20 shown).
5. **Export → CSV**. In the field chooser tick:
   - *Citation information* (authors, title, year, source, volume, pages, DOI, citation count)
   - *Bibliographical information* (affiliations, publisher)
   - *Abstract & keywords* ← **this one matters most**, the screening is done on abstracts
   - *Include references* — leave **off**, it makes the file unusable in a spreadsheet
6. Save as `scopus\A1-cited-by.csv` (and `A2-`, `A3-`).

If a DOI returns nothing, try the title in quotes; CGJ records are sometimes indexed without
the DOI field populated. Say so if it happens — it is itself worth reporting.

## 2. Web of Science

For **each** DOI:

1. **Advanced search**, query `DO=(10.1139/cgj-2016-0327)` with the DOI of the anchor.
2. Open the record, click **Times Cited**. **Write that number down**, with the date.
3. **Export → Excel** (or *Tab delimited file*), and for *Record Content* choose
   **Full Record** — not "Author, Title, Source", which omits the abstract.
4. Save as `wos\A1-cited-by.xls`.

WoS caps an export at 1000 records; you will be far below that, but if it ever asks, say so.

## 3. Google Scholar

Scholar has no bulk export, so this one is more manual and noisier. It is included because it
catches theses, reports and conference papers the other two miss.

1. Search the exact title in quotes.
2. Click **Cited by _N_**. **Write N down**, with the date.
3. Tick **Include citations** off, so you get documents rather than mentions.
4. If you have **Publish or Perish** (free, Harzing), point it at the same "cited by" list and
   export to CSV — that is much the easiest route.
5. If not: for each page of results click the quotation-mark icon → **BibTeX**, and paste
   them into one file `scholar\A1-cited-by.bib`. Stop at **200 records** — that cap is in the
   protocol and declaring it in advance is what makes it defensible.

## 4. What to send me

Just the folder, or the files. For each of the nine, I need to know:

* the **count** the source reported, and the **date** you searched
* the export file itself, **with abstracts**

If an export has no abstracts, the screening cannot be done from it and I will have to ask
you to redo it — so it is worth checking one file before doing all nine.

## 5. What I do with it

1. **Merge and deduplicate** across the nine, keeping track of which source each record came
   from — PRISMA wants duplicates counted and reported, not silently dropped.
2. Build the **extraction table** with the fields in the protocol, one row per record.
3. Give you back a **first screening**: which records look relevant to RQ2, which do not, and
   why — for you to check rather than to accept.
4. Report the numbers for the flow diagram: identified, duplicates removed, screened,
   excluded with reasons, included.

I will not decide inclusion on my own. I will propose, with the reason in each row, and
anything not clear-cut gets flagged for the two of us to decide — the protocol says that
number has to be reported, so it has to be real.

---

## A warning about what this step cannot tell us

A citation count is a measure of **visibility**, not of use. A paper that cites A1 in a list
of six PIV references has told us nothing about the method. That is why the extraction form
classifies *how* each work uses it — applied, extended, modified, criticised, or cited in
passing — rather than just counting.

Expect most of them to be in passing. That is a normal and reportable finding, and it is
better to say so than to inflate it.

A3 is from 2025 and will have very few citations. That is its age, not its reception, and the
results have to say so.
