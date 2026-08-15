# External sources and provenance

## Qiu et al.

- C. Qiu et al., “Why grain growth is not curvature flow,” *PNAS* 122 (2025), DOI `10.1073/pnas.2500707122`.
- Data and code: Zenodo record `10.5281/zenodo.15120372`.
- The Zenodo API reports `CC-BY-4.0` for the record and the local `PF_Codes.zip` matches the published MD5. The journal version is separately marked CC BY-NC-ND 4.0; equations are facts and are independently reimplemented with attribution.
- The pristine archive is retained in the workspace and extracted only under ignored `.external/qiu/`. No Qiu source file is copied into `src/`.

## Numerical and physical references

- D. Fan and L.-Q. Chen, “Computer simulation of grain growth using a continuum field model,” *Acta Materialia* 45 (1997) 611–622, DOI `10.1016/S1359-6454(96)00200-5`: independent continuum-field benchmark reporting the conventional two-dimensional scaling exponent (m=2). It was consulted when the original constrained independent-well implementation produced (n=1); the required filling constraint and the directive's Qiu reference favored correcting the pairwise MPF formulation rather than switching model families.
- I. Steinbach et al., “A phase field concept for multiphase systems,” *Physica D* 134 (1999) 385–393, DOI `10.1016/S0167-2789(99)00129-3`: constrained multiphase evolution.
- J. W. Cahn and J. E. Taylor, “A unified approach to motion of grain boundaries, relative tangential translation along grain boundaries, and grain rotation,” *Acta Materialia* 52 (2004) 4887–4898, DOI `10.1016/j.actamat.2004.02.048`: shear coupling kinematics.
- J. Han, S. L. Thomas, and D. J. Srolovitz, “Grain-boundary kinetics: A unified approach,” *Progress in Materials Science* 98 (2018) 386–476, DOI `10.1016/j.pmatsci.2018.05.004`: disconnection modes and thermally activated kinetics.
- NumPy/SciPy documentation: array finite differences, FFT, nonlinear least squares, special functions, and statistical tests.

Significant modeling choices are stated as closures, not attributed as established crystallographic laws. In particular the isotropic mode barrier spectrum, local shear-memory energy, free-volume storage energy, and coarse-grained event packet are falsifiable surrogate assumptions.
