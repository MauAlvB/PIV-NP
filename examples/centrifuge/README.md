# Example: slope in a centrifuge

Data of the case from the paper (`caso paper centrifuga`): 149 PIVlab fields of 60 × 35
points, cells of 0.212115 m and Δt = 0.8 s.

The original `zapatak.PAR` of that folder is in the old format (the one of the 2022
`PIVNP.exe`, with 4 values in block 3). The one in this folder is the same case in the
current format (9 values in block 3, plus a block 4).

The data (`datos (1).txt` … `datos (149).txt`, 13 MB) is not in the repository. To run it,
copy the files here and launch:

```bash
pivnp examples/centrifuge
```

It produces `zapatak.POST.MSH`, `zapatak.POST.RES` (≈ 770 MB) and `zapatak.REC`.
