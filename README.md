EVM reproducibility package
==========================================

Manuscript: "An Extended van Genuchten-Mualem Model for Predicting the Unsaturated Hydraulic
Conductivity from Soil Water Characteristic Data" (Ghorbani, Sadeghi, Zhang, van Genuchten, Jones).

These MATLAB files reproduce every table, figure and number of the manuscript and of its
Supporting Information (Texts S1-S7, Figures S1-S12, Tables S1-S9). Each script writes a text
report with all the numbers it produces and prints a line ending in _DONE when it has finished.


1. Requirements
---------------
- MATLAB R2025b (tested). The code uses functions introduced up to R2020b (e.g. tiledlayout
  legends, exportgraphics); older releases were not tested.
- Optimization Toolbox (lsqcurvefit).
- Parallel Computing Toolbox (optional): with it, soils are fitted in parallel; without it,
  the same parfor loops run serially and take about 24 times longer.
- About 60 MB of free disk space for the outputs.


2. Data (place in the same folder as the scripts)
------------------------------------------------
RetMeas.csv, CondMeas.csv
    German soil hydraulic database (572 soils from TUDD, UFZ, FSUJ, KIT and TUBS):
    Hohenbrink, T.L., et al. (2023). Soil hydraulic characteristics in a wide range of
    saturation and soil properties. GFZ Data Services. https://doi.org/10.5880/fidgeo.2023.012
    Columns: Sample_ID, MeasType, pF, theta (cm3 cm-3) or k (cm d-1). Rows with MeasType KSAT
    (pF = -2) hold the measured saturated conductivity (409 soils), used as Ks,c.

Table1_literature_RetMeas.csv, Table1_literature_CondMeas.csv   (Data Set S2)
    The five literature soils of Table 1, in the same format (MeasType "Literature"):
        Berlin_Sand   Berlin sand (Peters, 2013)                      Ks = 673 cm d-1
        Pachepsky_SL  Pachepsky sandy loam (Pachepsky et al., 1984)   Ks = 12.8 cm d-1
        Gilat_Loam    Gilat loam (Mualem, 1976b)                      Ks = 17.28 cm d-1
        Pachepsky_CL  Pachepsky clay loam (Pachepsky et al., 1984)    Ks = 28.94 cm d-1
        Pachapa_FSC   Pachapa fine sandy clay (Mualem, 1976b)          Ks = 12.1 cm d-1
    The first four come from Data Sets S11, S12, S14 and S15 of Rudiyanto et al. (2015), Water
    Resour. Res. 51, 8757-8772, https://doi.org/10.1002/2015WR017703; the clay loam comes from the
    authors' data file clay_loam_Pachepsky_1984.csv (suction in m and K in m/s, converted to cm
    and cm/d). Ks (the KSAT rows) is the capillary Ks of Table 1. Scripts S5, S6 and S8 use these
    files. EVM_literature_data(srcDir) rebuilds them from the downloaded files in srcDir.


3. How to run
-------------
In MATLAB, change to this folder and type

    EVM_run_all              % all eight scripts, about 1 h with 24 parallel workers
    EVM_run_all(pwd, 20)     % quick test: the first 20 German soils in every script

or run any script on its own, e.g.  R = EVM_reproduce_manuscript;
Every script accepts (dataDir, nSoilsMax); nSoilsMax may also be a list of Sample_IDs, e.g.
EVM_VG_comparison(pwd, ["TUBS_205" "FSUJ_015"]). Script S6 compares its fits with the outputs of
Scripts S2 and S4 (if present), Script S7 with those of Script S6, and Script S8 with those of
Scripts S1, S6 and S7, so run them after them. EVM_beta_closure_analysis(pwd, 'figures'),
EVM_reproduce_manuscript(pwd, 'figures'), EVM_original_code_recheck(pwd, 'figures') and
EVM_beta_closure_recheck(pwd, 'figures') redraw their figures from the saved results without
refitting.


4. Files and what they reproduce
--------------------------------
EVM_models.m
    Reference implementation of the van Genuchten (VG), Peters-Durner-Iden (PDI) and extended
    van Genuchten-Mualem (EVM) models: theta(h), C(h) = -dtheta/dh, Sc and Snc, and
    K(h) = Kc + Knc + Kv. The exact equations are listed in its header. The scripts below contain
    the same functions (identical to rounding error; checked by Script S6).

Script S1  EVM_beta_closure_analysis.m           -> beta_closure_results\
    Monotonicity and smoothness of Eq. (6) and the soil-specific beta (Text S1).
        Figure S1  FigS1_Snc_monotonic_smoothness.png
        Figure S2  FigS2_soil_specific_beta.png
        EVM_soil_specific_beta.csv (soil-specific beta of all German soils), EVM_closure_scan.csv,
        EVM_closure_summary.csv, EVM_beta_closure_report.txt, EVM_beta_closure_results.mat
Script S2  EVM_positivity_constraint_analysis.m  -> positivity_results\
    EVM fits with theta_s > theta_r >= 0; positivity, endpoint and continuity diagnostics (Text S3).
        Figure S6  FigS6_endpoints_continuity.png
        Figure S7  FigS7_capacity_all_soils.png
        EVM_positivity_fits.csv (every soil; column e0 = epsilon), EVM_positivity_groups.csv,
        EVM_TUBS205_curves.csv, EVM_positivity_report.txt, EVM_positivity_results.mat
Script S3  EVM_suction_range_validation.m        -> range_results\
    Errors of VG, PDI and EVM by suction range, data coverage, overlap of the 12 Table 1 soils with
    the database (Table 4; Text S4; Tables S5 and S6).
        Figure S8  FigS8_range_coverage.png
        Figure S9  FigS9_range_errors.png
        Figure S10 FigS10_range_examples.png
        EVM_range_by_soil.csv, EVM_range_summary.csv, EVM_range_coverage.csv (Table S5),
        EVM_range_report.txt, EVM_range_results.mat
Script S4  EVM_VG_comparison.m                   -> vg_results\
    Conventional VG compared with PDI and EVM under the same constraints; endpoint behaviour;
    paired statistics (Table 5; Text S5; Table S7).
        Figure S11 FigS11_VG_comparison.png
        Figure S12 FigS12_VG_endpoints.png
        EVM_VG_by_soil.csv, EVM_VG_paired.csv (Table S7), EVM_VG_report.txt, EVM_VG_results.mat
Script S5  EVM_noncapillary_conductivity_check.m -> ncond_results\
    Derivation and accuracy of the non-capillary conductivity, Eq. (20) (Appendix B; Text S2):
    the power law of Eq. (20), the normalized integral with Eq. (19) (M19) and the exact inversion
    (Mex), checked against numerical integration, for the German soils and the five literature
    soils of Table 1 (Table S3).
        Figure S5  FigS5_noncapillary_conductivity.png
        EVM_ncond_suction_errors.csv (Table S3 for three values of R), EVM_ncond_by_soil.csv,
        EVM_ncond_grid.csv, EVM_ncond_report.txt, EVM_ncond_results.mat
Script S6  EVM_reproduce_manuscript.m            -> manuscript_results\
    Fits of VG, PDI and EVM with theta_s > theta_r >= 0 for the 572 German soils and the 12 Table 1
    soils: Tables 2 and 3 and Figure 1.
        Table2_Table1_soils.csv          Table 2 (with the values of the earlier version)
        Table3_groups.csv                Table 3
        EVM_parameter_sets.csv           parameters, errors and capacity diagnostics of every fit
        EVM_processed_observations.csv   Data Set S1: every observation used, with its role
        Fig1_Table1_soils_part1.png, Fig1_Table1_soils_part2.png (Figure 1)
        EVM_reproduction_report.txt, EVM_reproduction_results.mat
Script S7  EVM_original_code_recheck.m           -> recheck_results\
    The fits with the procedure of the earlier version of the manuscript (SWRFIT_PDI_EVM.m: one
    lsqcurvefit run per soil from x0 = [alpha theta_r n theta_s] = [0.01 0.001 1.5 0.01], independent
    bounds 0 <= alpha, theta_r, theta_s <= 1, 1 <= n <= 10), reproduced verbatim; they recover all
    published parameters of Table 2 and all entries of Table 3 of that version. Compares them with
    the fits with theta_s > theta_r >= 0 (Table 6; Figures 2 and 3; Texts S3, S6 and S7; Tables S4,
    S8 and S9). The output files label the two procedures "previous" and "constrained".
        EVM_recheck_fits.csv (every fit; Tables S4 and S8), EVM_recheck_table6.csv (Table 6),
        EVM_recheck_table3.csv (Table S9), EVM_recheck_summary.csv,
        Fig2_previous_vs_constrained_part1.png, Fig2_previous_vs_constrained_part2.png (Figure 2),
        Fig3_TUBS205_previous_vs_constrained.png (Figure 3),
        EVM_recheck_report.txt, EVM_recheck_results.mat
Script S8  EVM_beta_closure_recheck.m            -> beta_recheck_results\
    Closure beta = 1 + c/n (Text S1): for every c (0-100 and c -> Inf), the German soils and the
    12 Table 1 soils are refitted with independent bounds and one lsqcurvefit run (the procedure of
    the earlier version) and with theta_s > theta_r >= 0; retention error by suction range, sign and
    suction range of C < 0, RMSE of log10 K, soil-specific beta (Tables S1 and S2). The output files
    label the two procedures single_start and constrained.
        Figure S3  FigS3_closure_sensitivity.png
        Figure S4  FigS4_Table1_soils.png
        EVM_beta_recheck_scan.csv (Table S2, every fit), EVM_beta_recheck_development_set.csv
        (Table S1), EVM_beta_recheck_summary.csv, EVM_beta_recheck_limits.csv,
        EVM_beta_recheck_soil_beta.csv, EVM_beta_recheck_report.txt, EVM_beta_recheck_results.mat
EVM_run_all.m
    Runs Scripts S1-S8 in this order and lists the run times.
EVM_literature_data.m
    Builds Table1_literature_RetMeas.csv and Table1_literature_CondMeas.csv (Data Set S2).


5. Conventions
--------------
- h is the suction (the magnitude of the matric head) in cm, h >= 0; pF = log10 h; oven dryness
  h0 = 1e7 cm.
- Retention is fitted by unweighted least squares (Eq. 22) with lsqcurvefit, with
  theta_s > theta_r >= 0 enforced: theta_r = f theta_s, 0 <= f <= 0.9999, 0.01 <= theta_s <= 1,
  10^-5 <= alpha <= 1 cm-1, 1.0001 <= n <= 10, from 30 starting points (110 where the constraint is
  binding), with an interior-point check of the constrained optimum. A search with independent
  bounds aligns the constrained optimum inside Scripts S1-S4 and S6; it is not reported.
- Conductivity is predicted, not fitted: Ks,c = measured Ksat (counted once, Eq. 14),
  Ks,nc = 86.2 theta(1e5 cm) alpha^1.5 cm d-1, Knc by the power law of Eq. (20) with l = 0
  (EVM_models option KncForm = 'powerlaw', the default), Kv at 20 C (Appendix A).
- RMSE_theta is in cm3 cm-3. RMSE_K is the RMSE of log10 K (Eq. 24) and is dimensionless
  (log10 units, i.e. orders of magnitude). It is computed for the soils with measured Ksat; the
  Ksat value itself is not an observation in RMSE_K.
- Group results are means of the per-soil RMSE values.
- The earlier version of the manuscript used one lsqcurvefit run per soil from
  x0 = [0.01 0.001 1.5 0.01] with independent bounds (Script S7). Its Tables 2 and 3 computed
  RMSE_K with Ks = the first conductivity row of each soil (the wettest unsaturated K for the 163
  soils without KSAT), included the Ks point in RMSE_K, used Snc^0.5 (l = 0.5 in Eq. 20) in the EVM
  non-capillary conductivity, and truncated the printed values. The revised results use the
  procedure described above (l = 0).
