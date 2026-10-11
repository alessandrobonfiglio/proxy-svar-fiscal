# Proxy SVAR for US Fiscal Policy

MATLAB code for the narrative proxy SVAR of Mertens & Ravn (2014) applied to
US fiscal data (Caldara & Kamps dataset, 1950Q1–2006Q4). The code estimates the
dynamic effects of tax and government spending shocks on output and compares
two bootstrap methods for inference: the wild bootstrap and the moving block
bootstrap.

Macroeconometrics exam project, MSc in Economics, University of Bologna.

## How to run

1. Download the repository (**Code → Download ZIP**) and unzip it.
2. Open `Draft.m` in MATLAB and press **Run**.

No paths need to be edited: all files are located relative to `Draft.m`.
The random seed is fixed (`rng(1234)`), so results are identical on every run.

## Repository structure

```
Draft.m                     main script
CK_RESTUD_DATASET.xlsx      data (sheet FINAL)
additional files/
    doPVAR.m                VAR estimation and narrative proxy identification
    doPVARbs.m              wild bootstrap
    doPVARbs_MBB.m          moving block bootstrap
    testResiduals.m         residual diagnostics
```

## Model

- **Variables:** federal tax revenues, government spending, GDP
  (logs, real, per capita).
- **Specification:** VAR(4) with constant, linear and quadratic trend,
  and a dummy for 1975Q2.
- **Identification:**
  - The tax shock is identified using the Mertens & Ravn (2012) narrative
    series of unanticipated tax changes as an external instrument.
  - The spending shock is identified recursively: government spending
    does not respond to output within the quarter.
- **Normalization:** both shocks equal 1% of GDP. The tax shock is a tax cut.
- **Inference:** 95% bands from 1,000 bootstrap replications.
  - Wild bootstrap: Rademacher weights, applied jointly to residuals
    and proxy.
  - Moving block bootstrap: block length 8 quarters, with position-in-block
    centering. Draws are discarded if the resampled proxy has fewer than
    20 non-zero observations, if the first-stage F-statistic is below 10,
    or if estimation fails. The number of discards is printed.

## Output

- First-stage diagnostics for the proxy: R² and F-statistic.
- Estimates of the structural parameters: θ_G, θ_Y, γ_T, ζ_T, ζ_G.
- Residual diagnostics: Ljung-Box, Breusch-Godfrey, portmanteau, ACF and
  Q-Q plots.
- Impulse responses of output, government spending and tax revenues to a tax
  shock and to a spending shock, with wild and MBB bands.

## References

- Caldara, D. and Kamps, C. (2017). The Analytics of SVARs: A Unified Framework
  to Measure Fiscal Multipliers. *Review of Economic Studies*, 84(3).
- Jentsch, C. and Lunsford, K. G. (2019). The Dynamic Effects of Personal and
  Corporate Income Tax Changes in the United States: Comment.
  *American Economic Review*, 109(7).
- Mertens, K. and Ravn, M. O. (2012). Empirical Evidence on the Aggregate
  Effects of Anticipated and Unanticipated US Tax Policy Shocks.
  *American Economic Journal: Economic Policy*, 4(2).
- Mertens, K. and Ravn, M. O. (2014). A Reconciliation of SVAR and Narrative
  Estimates of Tax Multipliers. *Journal of Monetary Economics*, 68.

## Authors

Cristiano Celot, Alessandro Bonfiglio
