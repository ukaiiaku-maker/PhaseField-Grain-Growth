# Qiu lineage graph

```text
QIU_LEGACY_FORENSIC seed5101
  55930486: step 0 --------------------------> 9984
                    \ 55948258 replay: 9000 -> 9984
  overlap is an exact diagnostic replay, not a continuation edge

FFT_EIGENSTRAIN_V2 seed5101
  55932457: dt=0.04, step 0 -> 8379, terminal
  55950433: dt=0.02, step 0 -> running, separate refinement

QIU_SI_REFERENCE [100] 4ref
  20260912T180419Z-nogit-a25f0c: prepared, unsubmitted
```

No pair of jobs satisfies all ten continuation requirements. The authoritative
isotropic baseline is the pristine native `QIU_SI_REFERENCE`; it has no result
yet. `QIU_LEGACY_FORENSIC` and `FFT_EIGENSTRAIN_V2` remain separate identities.
