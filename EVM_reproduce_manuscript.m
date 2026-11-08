function R = EVM_reproduce_manuscript(dataDir, nSoilsMax)
%EVM_REPRODUCE_MANUSCRIPT  Reproduces Tables 2 and 3 and Figure 1 of the manuscript with the
% revised fitting procedure (theta_s > theta_r >= 0) and writes the parameter sets of every fit and
% the processed observations (Text S7 of the Supporting Information). Table 6 and Figures 2 and 3,
% which compare the previous and the revised fits, are produced by Script S7.
%
%   R = EVM_reproduce_manuscript              % data files next to this file
%   R = EVM_reproduce_manuscript(dataDir)     % data files in dataDir
%   R = EVM_reproduce_manuscript(dataDir, 30) % quick test on the first 30 German soils
%   R = EVM_reproduce_manuscript(dataDir, 'figures')   % redraw the figures from the saved results
%
% Data
%   RetMeas.csv, CondMeas.csv: German database of Hohenbrink et al. (2023), 572 soils.
%   Table1_literature_RetMeas.csv, Table1_literature_CondMeas.csv: the Table 1 literature soils
%   (same format; Data Set S2, written by EVM_literature_data.m): Berlin sand (Peters, 2013), Gilat
%   loam and Pachapa fine sandy clay (Mualem, 1976b) and Pachepsky sandy loam (Pachepsky et al.,
%   1984), from Data Sets S11, S12, S14 and S15 of Rudiyanto et al. (2015, Water Resour. Res., 51,
%   8757-8772, doi:10.1002/2015WR017703), and the Pachepsky clay loam (Pachepsky et al., 1984; the
%   authors' data file); KSAT = capillary Ks of Table 1.
% Models (exact implementations, identical to EVM_models.m):
%   VG : theta = thr + (ths - thr) G, G = [1 + (alpha h)^n]^-m, m = 1 - 1/n; Mualem K (capillary)
%   PDI: theta = (ths - thr) Sc + thr Snc, Sc = (G - G0)/(1 - G0); Snc of Iden and Durner (2014),
%        b = 0.1 + 0.2/n^2 [1 - exp(-(thr/(ths - thr))^2)]; K = Kc (Eq. 14) + Knc (Eqs. 16-17) + Kv
%   EVM: Eq. (7) with beta = 1 + 0.1/n; K = Kc (Eq. 14) + Knc (Eq. 20) + Kv (Appendix A)
%   h0 = 1e7 cm, ha = 1/alpha, Ks,c = KSAT, Ks,nc = 86.2 theta(1e5 cm) ha^-1.5, Kv at 20 C.
%   S.KncForm = 'powerlaw' uses Eq. (20) for EVM; 'integral' the normalized integral with Eq. (19) (Eq. B2).
% Fitting (Eq. 22, unweighted): theta_s > theta_r >= 0 enforced (P2: theta_r = f theta_s, f <= 0.9999),
%   log10(alpha) in [-5, 0], n in [1.0001, 10]; 30 starting points (110 where the constraint is
%   active), aligned internally with an independent-bounds search and checked with the
%   linear-inequality form of the constraint (Script S2).
% Negative water capacity: C = -dtheta/dh evaluated analytically on 20,000 suctions from
%   1e-12 cm to h0, plus the sign of the leading term as h -> 0+.
%
% Outputs (folder manuscript_results next to the data):
%   EVM_reproduction_report.txt       every number of Text S7 (Tables 2 and 3)
%   Table2_Table1_soils.csv           Table 2 (PDI, EVM and VG parameters and errors; previous values)
%   Table3_groups.csv                 Table 3 (mean errors by soil group)
%   EVM_parameter_sets.csv            parameters, errors and capacity diagnostics of every fit
%                                     (soil x model)
%   EVM_processed_observations.csv    Data Set S1: all observations used, with their role in the analysis
%   Fig1_Table1_soils_part1.png, Fig1_Table1_soils_part2.png, .mat
% Runtime: about 7 min with 24 parallel workers (R2025b, Optimization and Parallel
% Computing Toolboxes; parfor runs serially without a pool).

if nargin < 1 || isempty(dataDir)
    dataDir = fileparts(mfilename('fullpath'));
end
outDir = fullfile(dataDir, 'manuscript_results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
tStart = tic;
if nargin > 1 && (ischar(nSoilsMax) || isstring(nSoilsMax)) && any(strcmpi(nSoilsMax, 'figures'))
    L = load(fullfile(outDir, 'EVM_reproduction_results.mat'), 'R');      % redraw the figures only
    make_figures(outDir, L.R.settings, L.R.bySoil, L.R.data, L.R.ids);
    R = L.R; fprintf('Figures redrawn in %s.\n', outDir);
    return
end

%% ---------------- settings ----------------
S.h0 = 1e7; S.c = 0.1; S.fMax = 0.9999; S.nMin = 1.0001; S.hmin = 1e-2;
S.KncForm = 'powerlaw';                                   % Eq. (20); 'integral' = Eq. (B2)
S.models = ["VG" "PDI" "EVM"]; S.procs = ["P1" "P2"];
S.hg = logspace(-12, log10(S.h0) - 1e-9, 20000)';
S.groups = ["TUDD" "UFZ" "FSUJ" "KIT" "TUBS"];
S.T1 = table(["Sand"; "Loamy Sand"; "Sandy Loam"; "Loam"; "Silt"; "Silt Loam"; "Sandy Clay Loam"; "Clay Loam"; ...
    "Silty Clay Loam"; "Sandy Clay"; "Silty Clay"; "Clay"], ...
    ["Berlin_Sand"; "FSUJ_023"; "Pachepsky_SL"; "Gilat_Loam"; "FSUJ_035"; "TUBS_027"; "TUBS_128"; "Pachepsky_CL"; ...
    "KIT_022"; "Pachapa_FSC"; "KIT_009"; "TUBS_106"], ...
    ["Berlin"; "FSUJ-023"; "Pachepsky"; "Gilat"; "FSUJ-035"; "TUBS-027"; "TUBS-128"; "Pachepsky"; "KIT-022"; ...
    "Pachapa Fine"; "KIT-009"; "TUBS-106"], 'VariableNames', {'Texture', 'ID', 'Label'});
% manuscript Table 2 (theta_r, theta_s, alpha, n, RMSE_theta, RMSE_K), PDI and EVM, Table 1 order
S.ms2.PDI = [0.029 0.305 0.053 4.424 0.0052 0.13; 0.171 0.460 0.019 1.726 0.0037 0.35; 0.042 0.418 0.010 1.631 0.0073 0.39;
    0.152 0.437 0.017 3.353 0.0059 0.24; 9e-5 0.589 0.003 1.378 0.0117 0.93; 0.081 0.400 0.003 1.291 0.0024 0.22;
    0.245 0.383 0.015 1.464 0.0020 0.28; 0.015 0.494 0.007 1.161 0.0083 0.30; 0.301 0.550 0.080 1.662 0.0042 0.27;
    0.103 0.325 0.010 2.605 0.0068 0.31; 0.278 0.658 0.045 1.401 0.0068 0.42; 0.486 0.601 0.020 1.680 0.0027 0.38];
S.ms2.EVM = [0.029 0.306 0.053 4.413 0.0052 0.08; 0.173 0.461 0.017 1.765 0.0038 0.28; 0.044 0.418 0.010 1.644 0.0073 0.41;
    0.153 0.438 0.016 3.350 0.0063 0.31; 5e-5 0.589 0.003 1.378 0.0117 0.93; 0.097 0.400 0.003 1.330 0.0026 0.11;
    0.249 0.383 0.011 1.581 0.0019 0.48; 0.050 0.494 0.007 1.175 0.0085 0.33; 0.303 0.551 0.070 1.725 0.0041 0.33;
    0.104 0.325 0.010 2.629 0.0068 0.32; 0.304 0.658 0.038 1.478 0.0068 0.53; 0.485 0.603 0.014 1.995 0.0030 0.37];
% manuscript Table 3 (RMSE_theta, RMSE_K by group)
S.ms3.PDI = [0.0096 0.0050 0.0085 0.0061 0.0048; 1.40 0.56 0.89 1.00 0.87];
S.ms3.EVM = [0.0096 0.0051 0.0085 0.0061 0.0048; 1.43 0.62 0.95 1.07 0.94];
opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxFunctionEvaluations', 4000, 'MaxIterations', 1000, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-10);
optsLC = optimoptions('lsqcurvefit', 'Algorithm', 'interior-point', 'Display', 'off', ...
    'MaxFunctionEvaluations', 8000, 'MaxIterations', 2000, 'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-10, 'ConstraintTolerance', 1e-9);

%% ---------------- data ----------------
ret  = readtable(fullfile(dataDir, 'RetMeas.csv'),  'TextType', 'string');
cond = readtable(fullfile(dataDir, 'CondMeas.csv'), 'TextType', 'string');
ret.Source = repmat("Hohenbrink et al. (2023)", height(ret), 1); cond.Source = repmat("Hohenbrink et al. (2023)", height(cond), 1);
gid = unique(ret.Sample_ID, 'stable');
if nargin > 1 && ~isempty(nSoilsMax)
    if isnumeric(nSoilsMax), gid = gid(1:min(nSoilsMax, numel(gid)));
    else,                    gid = gid(ismember(gid, string(nSoilsMax)));
    end
    gid = unique([gid; intersect([S.T1.ID; "TUBS_205"], unique(ret.Sample_ID), 'stable')], 'stable');   % Table 1, Fig. 3
end
fr = fullfile(dataDir, 'Table1_literature_RetMeas.csv'); fk = fullfile(dataDir, 'Table1_literature_CondMeas.csv');
if exist(fr, 'file') && exist(fk, 'file')
    lr = readtable(fr, 'TextType', 'string'); lk = readtable(fk, 'TextType', 'string');
    lr.Source = repmat("Rudiyanto et al. (2015) SI", height(lr), 1); lk.Source = repmat("Rudiyanto et al. (2015) SI", height(lk), 1);
    cl = "Pachepsky et al. (1984)";                          % not among the Rudiyanto et al. data sets
    lr.Source(lr.Sample_ID == "Pachepsky_CL") = cl; lk.Source(lk.Sample_ID == "Pachepsky_CL") = cl;
    ret = [ret; lr]; cond = [cond; lk];
    lid = unique(lr.Sample_ID, 'stable');
else
    lid = strings(0, 1);
    warning('Literature data files not found: Table 2 and Figure 1 contain the German soils only.');
end
ret = ret(isfinite(ret.pF) & isfinite(ret.theta), :);
ids = [gid; lid]; nS = numel(ids);
[~, ir] = ismember(ret.Sample_ID, ids);
[~, ic] = ismember(cond.Sample_ID, ids);
D = struct('id', cellstr(ids), 'pF', [], 'th', [], 'mt', [], 'pFK', [], 'K', [], 'mtK', [], 'Ksat', NaN, 'src', "");
for i = 1:nS
    s = ir == i;
    D(i).pF = ret.pF(s); D(i).th = ret.theta(s); D(i).mt = ret.MeasType(s); D(i).src = ret.Source(find(s, 1));
    isKs  = ic == i & cond.MeasType == "KSAT";
    unsat = ic == i & cond.MeasType ~= "KSAT" & cond.k > 0;
    D(i).pFK = cond.pF(unsat); D(i).K = cond.k(unsat); D(i).mtK = cond.MeasType(unsat);
    if any(isKs), D(i).Ksat = cond.k(find(isKs, 1)); end
end
grp = extractBefore(ids, "_"); grp(ismember(ids, lid)) = "LIT";
fprintf('Loaded %d German and %d literature soils, %d with KSAT.\n', numel(gid), numel(lid), sum(isfinite([D.Ksat])));

if isempty(gcp('nocreate'))          % parfor runs serially when no pool can start
    try
        pc = parcluster('Processes');
        pc.NumWorkers = max(pc.NumWorkers, feature('numcores'));
        parpool(pc, feature('numcores'));
    catch
        try, parpool('Processes'); catch, end %#ok<NOCOM>
    end
end

%% ---------------- fits and diagnostics ----------------
out = cell(nS, 1);
parfor i = 1:nS
    out{i} = analyse_soil(D(i), S, opts, optsLC);
end
fprintf('Fits finished (%.1f min).\n', toc(tStart)/60);
T = struct2table(cell2mat(cellfun(@(o) o.row, out, 'UniformOutput', false)));
hMinObs = arrayfun(@(d) 10^min(d.pF), D);
T = [table(ids, grp, [D.src]', [D.Ksat]', hMinObs, 'VariableNames', {'Sample_ID', 'Group', 'Source', 'Ksat_cm_d', ...
    'hMinObs_cm'}), T];
germ = T.Group ~= "LIT";

% long table of the parameter sets (EVM_parameter_sets.csv)
P = table();
for M = S.models
    for pr = S.procs
        px = char(M + "_" + pr + "_");
        tt = table(T.Sample_ID, T.Group, T.Source, repmat(M, nS, 1), repmat(pr, nS, 1), T.([px 'alpha']), T.([px 'n']), ...
            T.([px 'theta_r']), T.([px 'theta_s']), T.([px 'RMSEth']), T.([px 'RMSEK']), T.Ksat_cm_d, T.([px 'thrGtThs']), ...
            ~T.([px 'strictPos']), T.([px 'negMeas']), T.([px 'neg1cm']), T.([px 'negObs']), T.hMinObs_cm, ...
            T.([px 'hNegMin']), T.([px 'hNegMax']), T.([px 'Cmin']), T.([px 'Cmax']), T.([px 'thetaExcess']), ...
            T.([px 'thetaExcessObs']), T.([px 'nObsNeg']), ...
            'VariableNames', {'Sample_ID', 'Group', 'Source', 'Model', 'Procedure', 'alpha_per_cm', 'n', 'theta_r', 'theta_s', ...
            'RMSE_theta', 'RMSE_log10K', 'Ksat_cm_d', 'theta_r_gt_theta_s', 'C_negative_somewhere', 'C_negative_h_ge_0p01cm', ...
            'C_negative_h_ge_1cm', 'C_negative_in_measured_range', 'h_wettest_retention_obs_cm', 'hNeg_min_cm', ...
            'hNeg_max_cm', 'C_min_per_cm', 'C_max_per_cm', 'theta_max_minus_theta_s', ...
            'theta_max_minus_theta_s_measured_range', 'n_retention_obs_with_C_negative'});
        P = [P; tt]; %#ok<AGROW>
    end
end
% Table 3 by group (German soils)
rows = cell(0, 7);
for pr = S.procs
    for g = [S.groups "All"]
        sel = germ & (T.Group == g | g == "All");
        for M = S.models
            px = char(M + "_" + pr + "_");
            rows(end + 1, :) = {pr, g, M, sum(sel), sum(sel & isfinite(T.([px 'RMSEK']))), mean(T.([px 'RMSEth'])(sel)), ...
                mean(T.([px 'RMSEK'])(sel), 'omitnan')}; %#ok<AGROW>
        end
    end
end
T3 = cell2table(rows, 'VariableNames', {'Procedure', 'Group', 'Model', 'Nsoils', 'NsoilsKsat', 'meanRMSE_theta', 'meanRMSE_log10K'});
for v = ["Procedure" "Group" "Model"], T3.(v) = string(T3.(v)); end

% Table 2 (Table 1 soils, P2 = revised procedure; P1 and the manuscript values for comparison)
rows = cell(0, 19);
for k = 1:height(S.T1)
    i = find(T.Sample_ID == S.T1.ID(k), 1);
    for M = S.models
        v = nan(1, 6); v1 = nan(1, 2); act = false;
        if ~isempty(i)
            px = char(M + "_P2_");
            v = [T.([px 'theta_r'])(i) T.([px 'theta_s'])(i) T.([px 'alpha'])(i) T.([px 'n'])(i) T.([px 'RMSEth'])(i) T.([px 'RMSEK'])(i)];
            v1 = [T.(char(M + "_P1_RMSEth"))(i) T.(char(M + "_P1_RMSEK"))(i)];
            act = T.(char(M + "_active"))(i);
        end
        if M == "VG", m2 = nan(1, 6); else, m2 = S.ms2.(char(M))(k, :); end
        rows(end + 1, :) = [{S.T1.Texture(k), S.T1.ID(k), M, ~isempty(i), act}, num2cell(v), num2cell(v1), num2cell(m2)]; %#ok<AGROW>
    end
end
T2 = cell2table(rows, 'VariableNames', {'Texture', 'Soil', 'Model', 'data_available', 'constraint_active', 'theta_r', ...
    'theta_s', 'alpha_per_cm', 'n', 'RMSE_theta', 'RMSE_log10K', 'P1_RMSE_theta', 'P1_RMSE_log10K', 'ms_theta_r', ...
    'ms_theta_s', 'ms_alpha', 'ms_n', 'ms_RMSE_theta', 'ms_RMSE_K'});
for v = ["Texture" "Soil" "Model"], T2.(v) = string(T2.(v)); end

% processed observations, with their role in the analysis
vn = {'Sample_ID', 'Source', 'Group', 'Variable', 'MeasType', 'pF', 'h_cm', 'value', 'Use'};
Oc = cell(nS, 1);
for i = 1:nS
    d = D(i); id = string(d.id); hasK = isfinite(d.Ksat); nr = numel(d.pF); nk = numel(d.K);
    use = "not used (no measured Ks,c)"; if hasK, use = "RMSE of log10 K (Eq. 24)"; end
    t1 = table(repmat(id, nr, 1), repmat(d.src, nr, 1), repmat(grp(i), nr, 1), repmat("theta", nr, 1), d.mt, d.pF, ...
        10.^d.pF, d.th, repmat("retention fit (Eq. 22)", nr, 1), 'VariableNames', vn);
    t2 = table(repmat(id, nk, 1), repmat(d.src, nk, 1), repmat(grp(i), nk, 1), repmat("K_cm_per_d", nk, 1), d.mtK, ...
        d.pFK, 10.^d.pFK, d.K, repmat(use, nk, 1), 'VariableNames', vn);
    t3 = table(id, d.src, grp(i), "Ksat_cm_per_d", "KSAT", -2, 1e-2, d.Ksat, "Ks,c in Eq. (14)", 'VariableNames', vn);
    if ~hasK, t3 = t3([], :); end
    Oc{i} = [t1; t2; t3];
end
O = vertcat(Oc{:});

% checks against Scripts S2-S4 and the stand-alone implementation
chk = consistency_checks(dataDir, T, S);

%% ---------------- outputs (constrained fits only; P1 is internal) ----------------
T2w = T2(:, ~ismember(T2.Properties.VariableNames, {'constraint_active', 'P1_RMSE_theta', 'P1_RMSE_log10K'}));
T3w = T3(T3.Procedure == "P2", :);
Pw  = P(P.Procedure == "P2", :);
vn  = T.Properties.VariableNames;
Tw  = T(:, ~contains(vn, '_P1_') & ~endsWith(vn, '_active') & ~endsWith(vn, '_nFix'));
writetable(T2w, fullfile(outDir, 'Table2_Table1_soils.csv'));
writetable(T3w, fullfile(outDir, 'Table3_groups.csv'));
writetable(Pw,  fullfile(outDir, 'EVM_parameter_sets.csv'));
writetable(O,   fullfile(outDir, 'EVM_processed_observations.csv'));
fid = fopen(fullfile(outDir, 'EVM_reproduction_report.txt'), 'w');
if fid < 0, error('Cannot write the report in %s', outDir); end
for f = [1, fid]
    print_report(f, S, T, P, T2, T3, chk, gid, lid);
end
fclose(fid);
R = struct('settings', S, 'bySoil', Tw, 'parameterSets', Pw, 'table2', T2w, 'table3', T3w, 'checks', chk, ...
    'data', D, 'ids', ids);
save(fullfile(outDir, 'EVM_reproduction_results.mat'), 'R');
make_figures(outDir, S, Tw, D, ids);
fprintf('Outputs written to %s (%.1f min).\n', outDir, toc(tStart)/60);
fprintf('EVM_REPRODUCTION_DONE\n');
end

%% =====================================================================
%% local functions
%% =====================================================================
function v = med(x), if isempty(x), v = NaN; else, v = median(x); end, end
function v = mx(x),  if isempty(x), v = NaN; else, v = max(x);    end, end
function v = mn(x),  if isempty(x), v = NaN; else, v = min(x);    end, end

function o = analyse_soil(d, S, opts, optsLC)
% VG, PDI and EVM fits with independent bounds (P1) and theta_s > theta_r >= 0 (P2)
t0 = tic;
h = 10.^d.pF; hasK = isfinite(d.Ksat) && ~isempty(d.K);
for M = S.models
    mdl = char(M);
    F = fit_procedures(h, d.th, mdl, S, opts, optsLC);
    r.([mdl '_active']) = F.active; r.([mdl '_nFix']) = F.nFix;
    for pr = S.procs
        p = F.(pr).p; px = [mdl '_' char(pr) '_'];
        r.([px 'alpha']) = p(1); r.([px 'n']) = p(2); r.([px 'theta_r']) = p(3); r.([px 'theta_s']) = p(4);
        r.([px 'RMSEth']) = sqrt(mean((shp_theta(h, p, mdl, S) - d.th).^2));
        r.([px 'RMSEK']) = NaN;
        if hasK, r.([px 'RMSEK']) = sqrt(mean((shp_logK(10.^d.pFK, p, mdl, d.Ksat, S) - log10(d.K)).^2)); end
        c = cap_diag(p, mdl, S, h);
        for fn = fieldnames(c)'
            r.([px fn{1}]) = c.(fn{1});
        end
    end
end
r.tFit_s = toc(t0);
o.row = r;
end

function c = cap_diag(p, mdl, S, hObs)
% sign, size and suction range of negative water capacity on 1e-12 cm <= h < h0
C = shp_capacity(S.hg, p, mdl, S); th = shp_theta(S.hg, p, mdl, S);
neg = C < 0; ths = p(4);
c.thrGtThs = p(3) > ths;
c.sign0 = lead_sign(p, mdl, S);
c.strictPos = ~any(neg) && c.sign0 > 0;
c.negMeas = any(neg & S.hg >= S.hmin); c.neg1cm = any(neg & S.hg >= 1);
c.negObs = any(neg & S.hg >= min(hObs));                 % theta rises inside the measured suction range
c.hNegMin = NaN; c.hNegMax = NaN; c.Cmin = NaN;
if any(neg), c.hNegMin = min(S.hg(neg)); c.hNegMax = max(S.hg(neg)); c.Cmin = min(C(neg)); end
if c.sign0 < 0 && ~any(neg), c.hNegMin = 0; c.hNegMax = S.hg(1); end      % only below 1e-12 cm
c.Cmax = max(C);
c.thetaExcess = max(max(th) - ths, 0); if c.thetaExcess < 1e-12, c.thetaExcess = 0; end
m = S.hg >= min(hObs);                                   % theta above theta_s at measured suctions
c.thetaExcessObs = max(max(th(m)) - ths, 0); if c.thetaExcessObs < 1e-12, c.thetaExcessObs = 0; end
c.nObsNeg = sum(shp_capacity(hObs, p, mdl, S) < 0);
end

function s0 = lead_sign(p, mdl, S)
% sign of C as h -> 0+: C ~ A (alpha h)^ec + B (alpha h)^enc, the smaller exponent dominates
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n;
if strcmp(mdl, 'VG'), s0 = sign(ths - thr); return, end
G0 = (1 + (al*S.h0)^n)^(-m);
A = (ths - thr)*m*n/(1 - G0); ec = n - 1;
if strcmp(mdl, 'PDI')
    if ths == thr, r2 = Inf; else, r2 = (thr/(ths - thr))^2; end
    b = 0.1 + 0.2/n^2*(1 - exp(-r2));
    B = thr/(log(10)*log10(al*S.h0)); enc = 1/(b*log(10)) - 1;
else
    beta = 1 + S.c/n; B = thr*beta/log1p((al*S.h0)^beta); enc = beta - 1;
end
if A >= 0,        s0 = sign(A + B);
elseif enc < ec,  s0 = sign(B);
elseif enc > ec,  s0 = sign(A);
else,             s0 = sign(A + B);
end
end

function F = fit_procedures(h, th, mdl, S, opts, optsLC)
% P1 and P2 aligned as in Script S2 (the P2 set lies inside the P1 set)
[p1, s1] = fit_shp(h, th, mdl, 'P1', start_set('P1', th, false), S, opts, optsLC);
[p2, s2] = fit_shp(h, th, mdl, 'P2', start_set('P2', th, false), S, opts, optsLC);
if p1(3) > S.fMax*p1(4)
    [pb, sb] = fit_shp(h, th, mdl, 'P1', start_set('P1', th, true), S, opts, optsLC);
    if sb < s1, p1 = pb; s1 = sb; end
    [pb, sb] = fit_shp(h, th, mdl, 'P2', start_set('P2', th, true), S, opts, optsLC);
    if sb < s2, p2 = pb; s2 = sb; end
end
F.nFix = 0;
if s2 < s1*(1 - 1e-9), p1 = p2; s1 = s2; F.nFix = F.nFix + 1; end
F.active = p1(3) > S.fMax*p1(4);
if ~F.active && s1 < s2*(1 - 1e-9), p2 = p1; s2 = s1; F.nFix = F.nFix + 1; end
if F.active
    [px, sx] = fit_shp(h, th, mdl, 'P2x', start_set('P2x', th, false), S, opts, optsLC);
    if sx < s2*(1 - 1e-9) && px(3) < px(4), p2 = px; F.nFix = F.nFix + 1; end
end
F.P1.p = p1; F.P2.p = p2;
end

function [p, sseBest] = fit_shp(h, th, mdl, pr, Q0, S, opts, optsLC)
% Eq. (22) with lsqcurvefit; q = [log10(alpha) n theta_r theta_s] (P1, P2x) or [log10(alpha) n theta_s f] (P2)
switch pr
    case 'P1',  lb = [-5 S.nMin 0 0];    ub = [0 10 1 1];
    case 'P2',  lb = [-5 S.nMin 0.01 0]; ub = [0 10 1 S.fMax];
    case 'P2x', lb = [-5 S.nMin 0 0.01]; ub = [0 10 1 1];
end
fun = @(q, hh) shp_theta(hh, q2p(q, pr), mdl, S);
ws = warning; warning('off', 'MATLAB:nearlySingularMatrix'); warning('off', 'MATLAB:singularMatrix');
sseBest = Inf; qBest = Q0(1, :);
for j = 1:size(Q0, 1)
    q0 = min(max(Q0(j, :), lb + 1e-9), ub - 1e-9);
    try
        if strcmp(pr, 'P2x')
            [q, sse] = lsqcurvefit(fun, q0, h, th, lb, ub, [0 0 1 -S.fMax], 0, [], [], [], optsLC);
            if q(3) >= q(4), continue, end
        else
            [q, sse] = lsqcurvefit(fun, q0, h, th, lb, ub, opts);
        end
    catch
        continue
    end
    if sse < sseBest, sseBest = sse; qBest = q; end
end
warning(ws);
p = q2p(qBest, pr);
end

function Q0 = start_set(pr, th, dense)
ths0 = min(1.02*max(th), 0.99);
[la, nn, rr] = ndgrid([-3 -2 -1.3], [1.15 1.4 2.0], [0.05 0.3]);
[la2, nn2, rr2] = ndgrid([-2.5 -1], [1.02 4], [0.01 0.6 0.97]);
la = [la(:); la2(:)]; nn = [nn(:); nn2(:)]; rr = [rr(:); rr2(:)];
if dense
    [la3, nn3, rr3] = ndgrid([-3.5 -2.5 -1.5 -0.5], [1.02 1.1 1.3 1.8 3], [0.01 0.3 0.7 0.97]);
    la = [la; la3(:)]; nn = [nn; nn3(:)]; rr = [rr; rr3(:)];
end
t0 = repmat(ths0, numel(la), 1);
switch pr
    case 'P1',  Q0 = [la nn rr t0];
    case 'P2',  Q0 = [la nn t0 rr];
    case 'P2x', Q0 = [la nn rr.*t0 t0];
end
end

function p = q2p(q, pr)
if strcmp(pr, 'P2'), p = [10^q(1) q(2) q(4)*q(3) q(3)];
else,                p = [10^q(1) q(2) q(3) q(4)];
end
end

function th = shp_theta(h, p, mdl, S)
% VG (van Genuchten, 1980), PDI (Iden and Durner, 2014) and EVM (Eq. 7)
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n;
G = (1 + (al*h).^n).^(-m);
if strcmp(mdl, 'VG')
    th = thr + (ths - thr)*G;
    return
end
G0 = (1 + (al*S.h0)^n)^(-m);
Sc = (G - G0)/(1 - G0);
if strcmp(mdl, 'PDI')
    Snc = pdi_snc(h, p, S.h0);
else
    beta = 1 + S.c/n;
    Snc = 1 - log1p((al*h).^beta)/log1p((al*S.h0)^beta);
end
th = (ths - thr)*Sc + thr*Snc;
end

function Snc = pdi_snc(h, p, h0)
% Snc = 1 + [x - xa + b ln(1 + exp((xa - x)/b))]/(xa - x0), x = log10 h, xa = log10 ha, x0 = log10 h0
al = p(1); n = p(2); thr = p(3); ths = p(4);
if ths == thr, r2 = Inf; else, r2 = (thr/(ths - thr))^2; end
b = 0.1 + 0.2/n^2*(1 - exp(-r2));
u = log10(al*h);
Snc = 1 - (max(u, 0) + b*log1p(exp(-abs(u)/b)))/log10(al*h0);
end

function C = shp_capacity(h, p, mdl, S)
% C = -dtheta/dh (cm-1) of the three models
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n; x = al*h;
dG = m*n*al*x.^(n - 1)./(1 + x.^n).^(m + 1);
if strcmp(mdl, 'VG'), C = (ths - thr)*dG; return, end
G0 = (1 + (al*S.h0)^n)^(-m);
if strcmp(mdl, 'PDI')
    if ths == thr, r2 = Inf; else, r2 = (thr/(ths - thr))^2; end
    b = 0.1 + 0.2/n^2*(1 - exp(-r2));
    dSnc = 1./(1 + exp(-log10(x)/b))./(h*log(10)*log10(al*S.h0));
else
    beta = 1 + S.c/n;
    dSnc = beta*al*x.^(beta - 1)./((1 + x.^beta)*log1p((al*S.h0)^beta));
end
C = (ths - thr)*dG/(1 - G0) + thr*dSnc;
end

function [lK, Kc, Knc, Kv] = shp_logK(h, p, mdl, Ksc, S)
% log10 K (cm/d); capillary term in logarithms (no underflow as n -> 1)
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n;
lG = -m*log1p((al*h).^n);
if strcmp(mdl, 'VG')
    lG0 = -Inf; lSc = lG;
else
    lG0 = -m*log1p((al*S.h0)^n);
    lSc = lG + log1p(-exp(lG0 - lG)) - log1p(-exp(lG0));
end
le = lG/m; ld = lG0/m;
u  = (exp(le) - exp(ld))/(1 - exp(ld));
lB = log(-expm1(m*log1p(-u)));
t  = u < 1e-10;
lB(t) = log(m) + le(t) + log1p(-exp(ld - le(t))) - log1p(-exp(ld));
lKc = log(Ksc) + 0.5*lSc + 2*lB;
Kc = exp(lKc);
if strcmp(mdl, 'VG')
    lK = lKc/log(10); Knc = zeros(size(h)); Kv = zeros(size(h));
    return
end
G = exp(lG); G0 = exp(lG0);
Sc = max((G - G0)/(1 - G0), 0);
if strcmp(mdl, 'PDI')
    Snc = max(pdi_snc(h, p, S.h0), 0);
else
    beta = 1 + S.c/n;
    Snc = max(1 - log1p((al*h).^beta)/log1p((al*S.h0)^beta), 0);
end
Ksnc = 86.2*shp_theta(1e5, p, mdl, S)*al^1.5;
Rv = al*S.h0;
if strcmp(mdl, 'EVM') && strcmpi(S.KncForm, 'integral')
    Knc = Ksnc*(expm1(0.75*Snc*log(Rv))/expm1(0.75*log(Rv))).^2;     % Eq. (B2), integral with Eq. (19)
else
    Knc = Ksnc*Rv.^(-1.5*(1 - Snc));                                  % Eq. (16) / Eq. (20)
end
Kv = vapor_K(h, (ths - thr)*Sc + thr*Snc, ths);
lK = log10(Kc + Knc + Kv);
end

function Kv = vapor_K(h, th, ths)
T = 293.15; Mw = 0.018015; g = 9.81; Rg = 8.314;
rho_sv = 1e-3*exp(31.3716 - 6014.79/T - 7.92495e-3*T)/T;
rho_w  = 998.2;
Da  = 2.12e-5*(T/273.15)^2*1e4*86400;
cH  = Mw*g/(Rg*T)/100;
tau = max(ths - th, 0).^(10/3)/ths^2;
Kv  = rho_sv/rho_w*tau*Da*cH.*exp(-cH*h);
end

function chk = consistency_checks(dataDir, T, S)
% P2 fits against Script S4, EVM P2 against Script S2, local functions against EVM_models.m
chk = struct('nS4', 0, 'dth', NaN, 'dK_VG_PDI', NaN, 'dK_EVM', NaN, 'nS2', 0, 'dS2', NaN, 'models_file', NaN, ...
    'negS4', NaN, 'negS6', NaN, 'negBoth', NaN);
f4 = fullfile(dataDir, 'vg_results', 'EVM_VG_by_soil.csv');
if exist(f4, 'file')
    T4 = readtable(f4, 'TextType', 'string');
    [tf, j] = ismember(T.Sample_ID, T4.Sample_ID);
    if any(tf)
        chk.nS4 = sum(tf); chk.dth = 0; chk.dK_VG_PDI = 0;
        for M = S.models
            chk.dth = max(chk.dth, max(abs(T.(char(M + "_P2_RMSEth"))(tf) - T4.(char(M + "_P2_RMSEth_all"))(j(tf)))));
            a = T.(char(M + "_P2_RMSEK"))(tf); b = T4.(char(M + "_P2_RMSEK_all"))(j(tf));
            if M == "EVM", chk.dK_EVM = max(abs(a - b), [], 'omitnan');
            else, chk.dK_VG_PDI = max(chk.dK_VG_PDI, max(abs(a - b), [], 'omitnan'));
            end
        end
    end
end
f2 = fullfile(dataDir, 'positivity_results', 'EVM_positivity_fits.csv');
if exist(f2, 'file')
    F2 = readtable(f2, 'TextType', 'string');
    [tf, j] = ismember(T.Sample_ID, F2.Sample_ID);
    chk.nS2 = sum(tf);
    chk.dS2 = max(abs(T.EVM_P2_RMSEth(tf) - F2.P2_RMSE_theta(j(tf))));
end
if exist('EVM_models', 'file') == 2
    h = logspace(-6, log10(S.h0) - 1e-6, 300)'; d = 0;
    for p = [0.02 1.5 0.1 0.45; 0.005 1.15 0.25 0.5; 0.08 3 0.05 0.35]'
        for M = S.models
            d = max(d, max(abs(EVM_models(M, 'theta', h, p') - shp_theta(h, p', char(M), S))));
            c1 = EVM_models(M, 'capacity', h, p'); c2 = shp_capacity(h, p', char(M), S);
            d = max(d, max(abs(c1 - c2)./max(abs(c2), 1e-300)));
            k1 = EVM_models(M, 'K', h, p', 'Ksc', 50, 'KncForm', S.KncForm); k2 = 10.^shp_logK(h, p', char(M), 50, S);
            d = max(d, max(abs(k1 - k2)./k2));
        end
    end
    chk.models_file = d;
end
end

function print_report(f, S, T, P, T2, T3, chk, gid, lid)
fprintf(f, '\n========== Reproduction of Tables 2 and 3 and Figure 1 with the revised fitting procedure ==========\n');
fprintf(f, 'Soils: %d German (Hohenbrink et al., 2023) + %d Table 1 literature soils (%s); EVM conductivity: %s form\n', ...
    numel(gid), numel(lid), strjoin(lid', ', '), S.KncForm);
fprintf(f, 'Fits: theta_s > theta_r >= 0 enforced (P2). C checked on %d suctions, 1e-12 cm to h0\n', numel(S.hg));
g = T(T.Group ~= "LIT", :);
fprintf(f, '\n--- Admissibility and mean errors, German database (%d soils) ---\n', numel(gid));
fprintf(f, '   fits with C < 0 somewhere: VG %d, PDI %d, EVM %d\n', sum(~g.VG_P2_strictPos), sum(~g.PDI_P2_strictPos), ...
    sum(~g.EVM_P2_strictPos));
fprintf(f, '   mean RMSE_theta: VG %.6f, PDI %.6f, EVM %.6f | mean RMSE_K: VG %.4f, PDI %.4f, EVM %.4f\n', mean(g.VG_P2_RMSEth), ...
    mean(g.PDI_P2_RMSEth), mean(g.EVM_P2_RMSEth), mean(g.VG_P2_RMSEK, 'omitnan'), mean(g.PDI_P2_RMSEK, 'omitnan'), ...
    mean(g.EVM_P2_RMSEK, 'omitnan'));

fprintf(f, '\n--- Table 3 (German soils): mean RMSE_theta / RMSE_K (manuscript PDI / EVM) ---\n');
for gg = [S.groups "All"]
    fprintf(f, '   %-5s', gg);
    for M = S.models
        a2 = T3(T3.Procedure == "P2" & T3.Group == gg & T3.Model == M, :);
        fprintf(f, ' | %s %.4f / %.3f', M, a2.meanRMSE_theta, a2.meanRMSE_log10K);
    end
    k = find(S.groups == gg);
    if ~isempty(k)
        fprintf(f, ' | ms PDI %.4f / %.2f, EVM %.4f / %.2f', S.ms3.PDI(1, k), S.ms3.PDI(2, k), S.ms3.EVM(1, k), S.ms3.EVM(2, k));
    end
    fprintf(f, ' | N %d (%d)\n', a2.Nsoils, a2.NsoilsKsat);
end

fprintf(f, '\n--- Table 2 (Table 1 soils): theta_r theta_s alpha n RMSE_theta RMSE_K | manuscript: same six ---\n');
for i = 1:height(T2)
    t = T2(i, :);
    if ~t.data_available
        if t.Model == "VG", fprintf(f, '   %-16s %-12s data not available\n', t.Texture, t.Soil); end
        continue
    end
    fprintf(f, '   %-16s %-12s %-3s %.3f %.3f %.4f %.3f %.4f %.2f | %.3f %.3f %.3f %.3f %.4f %.2f\n', t.Texture, t.Soil, t.Model, ...
        t.theta_r, t.theta_s, t.alpha_per_cm, t.n, t.RMSE_theta, t.RMSE_log10K, t.ms_theta_r, t.ms_theta_s, t.ms_alpha, ...
        t.ms_n, t.ms_RMSE_theta, t.ms_RMSE_K);
end
k = T2.data_available & T2.Model ~= "VG";
fprintf(f, '   Table 2 vs manuscript (PDI, EVM; %d rows): max |d RMSE_theta| %.4f, max |d RMSE_K| %.2f, max |d theta_r| %.3f, max |d theta_s| %.3f, max |d n| %.3f\n', ...
    sum(k), max(abs(T2.RMSE_theta(k) - T2.ms_RMSE_theta(k))), max(abs(T2.RMSE_log10K(k) - T2.ms_RMSE_K(k))), ...
    max(abs(T2.theta_r(k) - T2.ms_theta_r(k))), max(abs(T2.theta_s(k) - T2.ms_theta_s(k))), max(abs(T2.n(k) - T2.ms_n(k))));
for M = ["PDI" "EVM" "VG"]
    q = T2(T2.data_available & T2.Model == M, :);
    fprintf(f, '   %-3s mean over the %d soils: RMSE_theta %.4f, RMSE_K %.3f (manuscript %.4f, %.3f)\n', M, height(q), ...
        mean(q.RMSE_theta), mean(q.RMSE_log10K), mean(q.ms_RMSE_theta), mean(q.ms_RMSE_K));
end

fprintf(f, '\n--- Checks ---\n');
if chk.nS4 > 0
    fprintf(f, '   vs Script S4 (%d soils): max |d RMSE_theta| %.1e (three models) | max |d RMSE_K| VG/PDI %.1e, EVM %.1e (EVM form: %s)\n', ...
        chk.nS4, chk.dth, chk.dK_VG_PDI, chk.dK_EVM, S.KncForm);
end
if chk.nS2 > 0, fprintf(f, '   vs Script S2 (%d soils): EVM P2 max |d RMSE_theta| %.1e\n', chk.nS2, chk.dS2); end
if isfinite(chk.models_file), fprintf(f, '   local functions vs EVM_models.m: max relative difference %.1e\n', chk.models_file); end
fprintf(f, '   fitting time per soil: median %.1f s, max %.1f s\n', median(T.tFit_s), max(T.tFit_s));
end

function make_figures(outDir, S, T, D, ids)
% Figure 1. Every figure is smaller than the screen (1707 x 1067 px here):
% larger off-screen figures lose the top of their axes contents in exportgraphics.
blue = [0 0.2 0.9]; red = [0.85 0 0];
hh = logspace(-2, log10(S.h0) - 1e-9, 500)';
% Figure 1: Table 1 soils, retention fitted and conductivity predicted with theta_s > theta_r >= 0
for part = 1:2
    f1 = figure('Visible', 'off', 'Position', [10 10 1000 990]);
    tl = tiledlayout(f1, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    H = [];
    for j = 1:6
        k = (part - 1)*6 + j; i = find(ids == S.T1.ID(k), 1);
        ax = nexttile(tl);
        ttl = sprintf('%s: %s', S.T1.Texture(k), S.T1.Label(k));
        if isempty(i)
            text(ax, 0.5, 0.5, {'data not available'}, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 9);
            set(ax, 'XTick', [], 'YTick', []); box(ax, 'on'); title(ax, ttl, 'FontSize', 10); continue
        end
        d = D(i); pP = pget(T, "PDI", "P2", i); pE = pget(T, "EVM", "P2", i);
        yyaxis(ax, 'left'); hold(ax, 'on');
        hd = plot(ax, 10.^d.pF, d.th, 'o', 'Color', blue, 'MarkerSize', 3);
        hp = plot(ax, hh, shp_theta(hh, pP, 'PDI', S), ':', 'Color', blue, 'LineWidth', 1.8);
        he = plot(ax, hh, shp_theta(hh, pE, 'EVM', S), '-', 'Color', blue, 'LineWidth', 1.0);
        ylim(ax, [0 ceil(10*1.15*max([d.th; pP(4); pE(4)]))/10]);
        ylabel(ax, 'Water content, \theta (cm^3 cm^{-3})');
        yyaxis(ax, 'right'); hold(ax, 'on');
        Hj = [hd hp he];
        if isfinite(d.Ksat)
            hk  = plot(ax, 10.^d.pFK, d.K, 's', 'Color', red, 'MarkerSize', 4);
            hs  = plot(ax, 1e-2, d.Ksat, 's', 'Color', red, 'MarkerFaceColor', red, 'MarkerSize', 5);
            hkp = plot(ax, hh, 10.^shp_logK(hh, pP, 'PDI', d.Ksat, S), ':', 'Color', red, 'LineWidth', 1.8);
            hke = plot(ax, hh, 10.^shp_logK(hh, pE, 'EVM', d.Ksat, S), '-', 'Color', red, 'LineWidth', 1.0);
            set(ax, 'YScale', 'log'); ylim(ax, [1e-12 10^ceil(log10(max(d.Ksat, 1)) + 1)]);
            Hj = [Hj hk hs hkp hke]; %#ok<AGROW>
        end
        if numel(Hj) > numel(H), H = Hj; end
        ylabel(ax, 'Hydraulic conductivity, K (cm d^{-1})');
        ax.YAxis(1).Color = blue; ax.YAxis(2).Color = red;
        set(ax, 'XScale', 'log', 'XTick', 10.^(-2:7), 'FontSize', 8); xlim(ax, [1e-2 S.h0]); grid(ax, 'on'); box(ax, 'on');
        if j > 4, xlabel(ax, 'Suction, h (cm)'); end
        title(ax, ttl, 'FontSize', 10);
    end
    lg = legend(H, {'\theta data', '\theta fitted, PDI', '\theta fitted, EVM', 'K data', 'K_s measured (plotted at h = 0.01 cm)', ...
        'K predicted, PDI', 'K predicted, EVM'}, 'Orientation', 'horizontal', 'NumColumns', 4, 'FontSize', 8);
    lg.Layout.Tile = 'south';
    exportgraphics(f1, fullfile(outDir, sprintf('Fig1_Table1_soils_part%d.png', part)), 'Resolution', 300); close(f1);
end
end

function p = pget(T, M, pr, i)
p = [T.(char(M + "_" + pr + "_alpha"))(i) T.(char(M + "_" + pr + "_n"))(i) ...
     T.(char(M + "_" + pr + "_theta_r"))(i) T.(char(M + "_" + pr + "_theta_s"))(i)];
end
