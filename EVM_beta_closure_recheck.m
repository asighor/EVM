function R = EVM_beta_closure_recheck(dataDir, nSoilsMax)
%EVM_BETA_CLOSURE_RECHECK  Re-check of the closure beta = 1 + c/n (Script S1) with the original fitting
% code of the manuscript, SWRFIT_PDI_EVM.m, and with the revised fitting procedure (Reviewer #2,
% comment on beta; Reviewer #4, comment 2).
%
%   R = EVM_beta_closure_recheck              % data files next to this file
%   R = EVM_beta_closure_recheck(dataDir)
%   R = EVM_beta_closure_recheck(dataDir, 30) % quick test on the first 30 German soils (or a list of IDs)
%   R = EVM_beta_closure_recheck(dataDir, 'figures')   % redraw the figures from the saved results
%
% For every c of the grid (beta = 1 + c/n; c = 0-100 and the limit c -> Inf, where Eq. 6 becomes a
% broken line in log h with a kink at h = 1/alpha), alpha, n, theta_r and theta_s of every soil are
% refitted to theta(h) (Eq. 22, unweighted) with two procedures:
%   single_start  independent bounds on theta_r and theta_s: ONE lsqcurvefit run (default options) from
%        x0 = [alpha theta_r n theta_s x5] = [0.01 0.001 1.5 0.01 0.01], 0 <= alpha, theta_r, theta_s <= 1,
%        1 <= n <= 10 (line 199 of SWRFIT_PDI_EVM.m with 0.1 replaced by c; for c > 20 the same function is
%        evaluated without overflow). At c = 0.1 these are the published fits (Script S7).
%   constrained   theta_s > theta_r >= 0 enforced as theta_r = f theta_s, 0 <= f <= 0.9999, 0.01 <= theta_s <= 1;
%        18 starting points (30 at c = 0 and 0.1) plus continuation in c. The search also fits the same
%        starting points with independent bounds and keeps such a solution when it has a lower SSE and
%        satisfies the constraint (the constrained parameter set lies inside the unconstrained one, as in
%        Script S1); these unconstrained fits are not reported.
% Internally the two procedures are labelled OB1 and P2 (P1 = the unconstrained part of the search).
% Data: German database (572 soils) and the development set, i.e. the 12 soils of Table 1
% (7 German soils and the 5 literature soils of Data Set S2).
% For every fit: RMSE_theta (all data; wet pF < 2, intermediate 2 <= pF < 4.2, dry pF >= 4.2, as in
% Script S3), the sign of C = -dtheta/dh on 1e-12 cm <= h < h0 and as h -> 0+, and the RMSE of
% log10 K (Eq. 20 with l = 0, Ks,c = measured KSAT, unsaturated K data; as in Script S6).
% Soil-specific beta (constrained, 1 <= beta <= 10) of every soil, with the BIC against the closure.
%
% Outputs (folder beta_recheck_results next to the data): EVM_beta_recheck_report.txt (every number of
% Text S1), EVM_beta_recheck_scan.csv (Table S2: every fit), EVM_beta_recheck_summary.csv (per procedure
% and c), EVM_beta_recheck_limits.csv (per soil: largest c with C > 0 under independent bounds),
% EVM_beta_recheck_development_set.csv, EVM_beta_recheck_soil_beta.csv, FigS3_closure_sensitivity.png,
% FigS4_Table1_soils.png, EVM_beta_recheck_results.mat.
% Runtime: about 9 min with 24 parallel workers (R2025b, Optimization and Parallel Computing Toolboxes).

if nargin < 1 || isempty(dataDir)
    dataDir = fileparts(mfilename('fullpath'));
end
outDir = fullfile(dataDir, 'beta_recheck_results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
tStart = tic;
if nargin > 1 && (ischar(nSoilsMax) || isstring(nSoilsMax)) && any(strcmpi(nSoilsMax, 'figures'))
    L = load(fullfile(outDir, 'EVM_beta_recheck_results.mat'), 'R');       % redraw the figures only
    make_figures(outDir, L.R.settings, L.R.summary, L.R.development, L.R.devSoils);
    R = L.R; fprintf('Figures redrawn in %s.\n', outDir);
    return
end

%% ---------------- settings ----------------
S.h0 = 1e7; S.hmin = 1e-2; S.hplot = 1e-6; S.fMax = 0.9999; S.nMin = 1.0001;
S.cRef = 0.10; S.betaMax = 10; S.cVerbatim = 20;
S.cGrid = [0 0.01 0.02 0.05 0.08 0.09 0.10 0.11 0.12 0.15 0.20 0.30 0.50 1 2 5 10 20 50 100 Inf];
S.hg = logspace(-12, log10(S.h0) - 1e-9, 20000)';      % 1e-12 cm <= h < h0
S.pFedge = [2 4.2];                                     % wet | intermediate | dry (Script S3)
S.procs = ["OB1" "P2"];          % reported: single start (independent bounds) and constrained
S.x0 = [0.01; 0.001; 1.5; 0.01; 0.01]; S.lb0 = [0; 0; 1; 0; 0]; S.ub0 = [1; 1; 10; 1; 1];
S.T1 = ["Berlin_Sand" "FSUJ_023" "Pachepsky_SL" "Gilat_Loam" "FSUJ_035" "TUBS_027" "TUBS_128" ...
    "Pachepsky_CL" "KIT_022" "Pachapa_FSC" "KIT_009" "TUBS_106"];     % Table 1 order
S.keep = ["FSUJ_015" "TUBS_011" "TUBS_205"];          % kept in quick tests
opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxFunctionEvaluations', 4000, 'MaxIterations', 1000, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-10);
nc = numel(S.cGrid);
kRef = find(abs(S.cGrid - S.cRef) < 1e-12);
k5 = find(S.cGrid == 5);

%% ---------------- data ----------------
ret  = readtable(fullfile(dataDir, 'RetMeas.csv'),  'TextType', 'string');
cond = readtable(fullfile(dataDir, 'CondMeas.csv'), 'TextType', 'string');
gid = unique(ret.Sample_ID, 'stable');
if nargin > 1 && ~isempty(nSoilsMax)
    if isnumeric(nSoilsMax), gid = gid(1:min(nSoilsMax, numel(gid)));
    else,                    gid = gid(ismember(gid, string(nSoilsMax)));
    end
    gid = unique([gid; intersect([S.T1 S.keep]', unique(ret.Sample_ID), 'stable')], 'stable');
end
lid = strings(0, 1);
fr = fullfile(dataDir, 'Table1_literature_RetMeas.csv'); fk = fullfile(dataDir, 'Table1_literature_CondMeas.csv');
if exist(fr, 'file') && exist(fk, 'file')
    lr = readtable(fr, 'TextType', 'string'); lk = readtable(fk, 'TextType', 'string');
    ret = [ret(:, {'Sample_ID', 'MeasType', 'pF', 'theta'}); lr(:, {'Sample_ID', 'MeasType', 'pF', 'theta'})];
    cond = [cond(:, {'Sample_ID', 'MeasType', 'pF', 'k'}); lk(:, {'Sample_ID', 'MeasType', 'pF', 'k'})];
    lid = unique(lr.Sample_ID, 'stable');
else
    warning('Literature data files not found: the development set contains the German soils only.');
end
ret = ret(isfinite(ret.pF) & isfinite(ret.theta), :);
ids = [gid; lid]; nS = numel(ids);
[~, ir] = ismember(ret.Sample_ID, ids);
[~, ic] = ismember(cond.Sample_ID, ids);
D = struct('id', cellstr(ids), 'german', true, 'pF', [], 'h', [], 'th', [], 'hK', [], 'K', [], 'Ksat', NaN);
for i = 1:nS
    s = ir == i;
    D(i).german = ~ismember(ids(i), lid);
    D(i).pF = ret.pF(s); D(i).h = 10.^D(i).pF; D(i).th = ret.theta(s);
    isKs  = ic == i & cond.MeasType == "KSAT";
    unsat = ic == i & cond.MeasType ~= "KSAT" & cond.k > 0;
    D(i).hK = 10.^cond.pF(unsat); D(i).K = cond.k(unsat);
    if any(isKs), D(i).Ksat = cond.k(find(isKs, 1)); end
end
isG = [D.german]'; isG = isG(:);
grp = extractBefore(ids, "_"); grp(~isG) = "LIT";
hasK = isfinite([D.Ksat]'); hasK = hasK(:) & reshape(arrayfun(@(d) ~isempty(d.K), D), [], 1);
npts = reshape(arrayfun(@(d) numel(d.th), D), [], 1);
fprintf('Loaded %d German and %d literature soils (%d German soils with KSAT).\n', sum(isG), sum(~isG), sum(hasK & isG));

if isempty(gcp('nocreate'))          % parfor runs serially when no pool can start
    try
        pc = parcluster('Processes');
        pc.NumWorkers = max(pc.NumWorkers, feature('numcores'));
        parpool(pc, feature('numcores'));
    catch
        try, parpool('Processes'); catch, end %#ok<NOCOM>
    end
end

%% ---------------- fits ----------------
out = cell(nS, 1);
parfor i = 1:nS
    out{i} = analyse_soil(D(i), S, opts);
end
fprintf('Fits finished (%.1f min).\n', toc(tStart)/60);

% nS x nc matrices of every diagnostic, by procedure
fl = {'alpha', 'n', 'thr', 'ths', 'sse', 'rmse', 'rmseW', 'rmseI', 'rmseD', 'thrGtThs', 'negAny', 'negPlot', ...
    'negMeas', 'neg1cm', 'negObs', 'hNegMax', 'rmseK', 'hStar', 'ef'};
A = struct();
for pr = S.procs
    for v = fl
        A.(pr).(v{1}) = cell2mat(cellfun(@(o) o.(pr).(v{1})(:)', out, 'UniformOutput', false));
    end
end

%% ---------------- summary by procedure and c (German soils) ----------------
rows = cell(0, 25);
for pr = S.procs
    a = A.(pr);
    for k = 1:nc
        ok = isG & isfinite(a.rmse(:, k)); kk = isG & hasK & isfinite(a.rmseK(:, k));
        o2 = ok & isfinite(a.rmse(:, kRef)); k2 = kk & isfinite(a.rmseK(:, kRef));   % changes relative to c = 0.1
        ng = ok & a.negMeas(:, k) == 1;
        rows(end + 1, :) = {pr, S.cGrid(k), sum(ok), mean(a.rmse(ok, k)), 100*(mean(a.rmse(o2, k))/mean(a.rmse(o2, kRef)) - 1), ...
            sqrt(sum(a.sse(ok, k))/sum(npts(ok))), mean(a.rmseW(ok, k), 'omitnan'), mean(a.rmseI(ok, k), 'omitnan'), ...
            mean(a.rmseD(ok, k), 'omitnan'), sum(a.thrGtThs(ok, k) == 1), sum(a.negAny(ok, k) == 1), sum(a.negPlot(ok, k) == 1), ...
            sum(ng), sum(a.neg1cm(ok, k) == 1), sum(a.negObs(ok, k) == 1), mx(a.hNegMax(ng, k)), sum(kk), ...
            mean(a.rmseK(kk, k)), 100*(mean(a.rmseK(k2, k))/mean(a.rmseK(k2, kRef)) - 1), ...
            median(a.rmseK(kk, k)), sum(a.ef(isG, k) == 0), md(a.hStar(ok & a.thr(:, k) >= 0.05, k)), ...
            sum(a.negMeas(ok, k) == 1 & a.thrGtThs(ok, k) == 1), sum(isG) - sum(ok), ...
            mx(a.hNegMax(ok & a.negAny(:, k) == 1, k))}; %#ok<AGROW>
    end
end
SUM = cell2table(rows, 'VariableNames', {'Procedure', 'c', 'Nfits', 'meanRMSE_theta', 'dMeanRMSE_theta_pct', ...
    'pooledRMSE_theta', 'meanRMSE_theta_wet', 'meanRMSE_theta_int', 'meanRMSE_theta_dry', 'N_theta_r_gt_theta_s', ...
    'N_C_negative_anywhere', 'N_C_negative_h_ge_1e_6cm', 'N_C_negative_h_ge_0p01cm', 'N_C_negative_h_ge_1cm', ...
    'N_C_negative_measured_range', 'max_hNeg_cm', 'NsoilsK', 'meanRMSE_log10K', 'dMeanRMSE_log10K_pct', ...
    'medianRMSE_log10K', 'N_exitflag0', 'median_hStar_cm_thr_ge_0p05', 'N_C_negative_with_theta_r_gt_theta_s', ...
    'N_failed', 'max_hNeg_any_cm'});
SUM.Procedure = string(SUM.Procedure);

% largest admissible c: C > 0 in every soil for every c' <= c (C < 0 anywhere includes h -> 0+)
crit = {'N_C_negative_anywhere', 'N_C_negative_h_ge_1e_6cm', 'N_C_negative_h_ge_0p01cm', 'N_C_negative_h_ge_1cm'};
LIMc = nan(numel(S.procs), numel(crit));
for ip = 1:numel(S.procs)
    q = SUM(SUM.Procedure == S.procs(ip), :);
    for j = 1:numel(crit)
        okc = cumprod(q.(crit{j}) == 0) == 1;
        if any(okc), LIMc(ip, j) = max(q.c(okc)); else, LIMc(ip, j) = -Inf; end
    end
end

% per soil: largest grid c with C > 0 for h >= 0.01 cm for all c' <= c (independent bounds)
LIM = table(ids, grp, 'VariableNames', {'Sample_ID', 'Group'});
for pr = "OB1"
    a = A.(pr); cm = inf(nS, 1); bad = nan(nS, 6);
    for i = 1:nS
        kb = find(a.negMeas(i, :) == 1, 1);
        if isempty(kb), continue, end
        if kb == 1, cm(i) = -Inf; else, cm(i) = S.cGrid(kb - 1); end
        bad(i, :) = [S.cGrid(kb) a.n(i, kb) a.thr(i, kb) a.ths(i, kb) a.hNegMax(i, kb) a.rmse(i, kb)];
    end
    LIM.(char("cmax_" + pr)) = cm;
    LIM.(char("first_violating_c_" + pr)) = bad(:, 1); LIM.(char("n_" + pr)) = bad(:, 2);
    LIM.(char("theta_r_" + pr)) = bad(:, 3); LIM.(char("theta_s_" + pr)) = bad(:, 4);
    LIM.(char("hNeg_max_cm_" + pr)) = bad(:, 5); LIM.(char("RMSE_theta_" + pr)) = bad(:, 6);
end

% P2: optimum of the pooled SSE over the grid and the best grid c of each soil
okP2 = isG & all(isfinite(A.P2.sse), 2);
[~, kOpt] = min(sum(A.P2.sse(okP2, :), 1));
[~, kBest] = min(A.P2.sse(okP2, :), [], 2);
bestCount = accumarray(kBest, 1, [nc 1])';
% the change of the fits with c (P2): soils whose RMSE_theta varies by less than 1e-4 over the grid
rangeP2 = max(A.P2.rmse(okP2, :), [], 2) - min(A.P2.rmse(okP2, :), [], 2);

%% ---------------- soil-specific beta (P2) and BIC ----------------
SB = struct2table(cell2mat(cellfun(@(o) o.free, out, 'UniformOutput', false)));
SB = [table(ids, grp, 'VariableNames', {'Sample_ID', 'Group'}), SB];
SB.RMSE_closure = A.P2.rmse(:, kRef); SB.RMSE_beta1 = A.P2.rmse(:, 1);
SB.n_closure = A.P2.n(:, kRef); SB.theta_r_closure = A.P2.thr(:, kRef);
SB.beta_closure = 1 + S.cRef./SB.n_closure;
sseC = A.P2.sse(:, kRef); sseF = min(SB.SSE_free, sseC);
SB.c_free = SB.n_free.*(SB.beta_opt - 1);                  % c that reproduces beta* with the fitted n
SB.dBIC = npts.*log(sseC./sseF) - log(npts);               % > 0: the fifth parameter is supported
SB.gain_pct = 100*(SB.RMSE_closure - SB.RMSE_free)./SB.RMSE_closure;
SB.beta_no_effect = SB.theta_r_free < 0.02;
SB.atLowerBound = SB.beta_opt <= 1.01; SB.atUpperBound = SB.beta_opt >= 0.99*S.betaMax;

%% ---------------- development set (Table 1 soils with data) ----------------
[tf, dv] = ismember(S.T1, ids); dv = dv(tf); devIds = S.T1(tf);
DEV = struct('c', S.cGrid);
for pr = S.procs
    a = A.(pr);
    DEV.pooled.(pr) = sqrt(sum(a.sse(dv, :), 1)./sum(npts(dv)));
    DEV.mean.(pr) = mean(a.rmse(dv, :), 1);
    DEV.negMeas.(pr) = sum(a.negMeas(dv, :) == 1, 1);
    DEV.negAny.(pr) = sum(a.negAny(dv, :) == 1, 1);
    DEV.thrGt.(pr) = sum(a.thrGtThs(dv, :) == 1, 1);
    [~, DEV.kOpt.(pr)] = min(DEV.pooled.(pr));
end
DS = table(devIds', grp(dv), A.OB1.n(dv, kRef), A.OB1.thr(dv, kRef), A.OB1.ths(dv, kRef), A.OB1.rmse(dv, kRef), ...
    A.P2.n(dv, kRef), A.P2.thr(dv, kRef), A.P2.ths(dv, kRef), A.P2.rmse(dv, kRef), A.P2.rmse(dv, 1), ...
    A.P2.rmse(dv, end), SB.n_free(dv), SB.theta_r_free(dv), SB.beta_opt(dv), SB.c_free(dv), SB.RMSE_free(dv), ...
    SB.dBIC(dv), LIM.cmax_OB1(dv), ...
    'VariableNames', {'Sample_ID', 'Group', 'n_OB1', 'theta_r_OB1', 'theta_s_OB1', 'RMSE_theta_OB1', 'n_P2', ...
    'theta_r_P2', 'theta_s_P2', 'RMSE_theta_P2', 'RMSE_theta_P2_beta1', 'RMSE_theta_P2_cInf', 'n_free', ...
    'theta_r_free', 'beta_opt', 'c_free', 'RMSE_theta_free_beta', 'dBIC', 'cmax_OB1'});
devReg = ols_line(1./SB.n_free(dv), SB.beta_opt(dv) - 1);
selR = dv(SB.theta_r_free(dv) >= 0.02);
devReg2 = ols_line(1./SB.n_free(selR), SB.beta_opt(selR) - 1);
gReg = ols_line(1./SB.n_free(isG), SB.beta_opt(isG) - 1);

%% ---------------- consistency checks ----------------
chk = struct('nS7', 0, 'dS7_rmse', NaN, 'dS7_par', NaN, 'nS1', 0, 'dS1_rmse', NaN, 'dS1_q', NaN, 'nS1b', 0, ...
    'dS1_free', NaN, 'nS6', 0, 'dS6_rmse', NaN, 'nS6diff', 0, 'dS6_K', NaN, 'dVerbatim', NaN);
f7 = fullfile(dataDir, 'recheck_results', 'EVM_recheck_fits.csv');
if exist(f7, 'file')
    F7 = readtable(f7, 'TextType', 'string');
    F7 = F7(F7.Model == "EVM" & F7.Procedure == "previous", :);
    [t7, j7] = ismember(ids, F7.Sample_ID); t7 = t7 & isG;
    chk.nS7 = sum(t7);
    chk.dS7_rmse = max(abs(A.OB1.rmse(t7, kRef) - F7.RMSEth(j7(t7))));
    chk.dS7_par = max(max(abs([A.OB1.alpha(t7, kRef) A.OB1.n(t7, kRef) A.OB1.thr(t7, kRef) A.OB1.ths(t7, kRef)] - ...
        [F7.alpha(j7(t7)) F7.n(j7(t7)) F7.theta_r(j7(t7)) F7.theta_s(j7(t7))])));
end
f1 = fullfile(dataDir, 'beta_closure_results', 'EVM_closure_scan.csv');
if exist(f1, 'file')
    F1 = readtable(f1, 'TextType', 'string'); d = [];
    for pr = "P2"
        for k = find(S.cGrid <= 5)
            b = F1(F1.Procedure == pr & abs(F1.c - S.cGrid(k)) < 1e-12, :);
            [t1, j1] = ismember(ids, b.Sample_ID); t1 = t1 & isG;
            d = [d; abs(A.(pr).rmse(t1, k) - b.RMSE_theta(j1(t1)))]; %#ok<AGROW>
            chk.nS1 = chk.nS1 + sum(t1);
        end
    end
    chk.dS1_rmse = max(d); chk.dS1_q = quantile(d, 0.99);
end
fb = fullfile(dataDir, 'beta_closure_results', 'EVM_soil_specific_beta.csv');
if exist(fb, 'file')
    B1 = readtable(fb, 'TextType', 'string');
    [tb, jb] = ismember(ids, B1.Sample_ID); tb = tb & isG;
    chk.nS1b = sum(tb); chk.dS1_free = max(abs(SB.RMSE_free(tb) - B1.RMSE_free(jb(tb))));
end
f6 = fullfile(dataDir, 'manuscript_results', 'EVM_parameter_sets.csv');
if exist(f6, 'file')
    P6 = readtable(f6, 'TextType', 'string');
    P6 = P6(P6.Model == "EVM" & P6.Procedure == "P2", :);
    [t6, j6] = ismember(ids, P6.Sample_ID);
    chk.nS6 = sum(t6);
    chk.dS6_rmse = max(abs(A.P2.rmse(t6, kRef) - P6.RMSE_theta(j6(t6))));
    chk.nS6diff = sum(abs(A.P2.rmse(t6, kRef) - P6.RMSE_theta(j6(t6))) > 1e-6);
    k6 = t6 & isfinite(A.P2.rmseK(:, kRef));
    chk.dS6_K = max(abs(A.P2.rmseK(k6, kRef) - P6.RMSE_log10K(j6(k6))));
end
% the expression of the original code against the overflow-free function used beyond c = 20
hh = logspace(-2, 6.8, 200)'; dv0 = 0;
for c = S.cGrid(S.cGrid <= S.cVerbatim)
    for x = [0.01 0.05 1.05 0.40; 0.003 0.6 1.3 0.45; 0.1 0.2 2.5 0.35]'
        fv = orig_model(c, S); xx = [x; 0.01];
        dv0 = max(dv0, max(abs(fv(xx, hh) - evm_theta(hh, [x(1) x(3) x(2) x(4)], 1 + c/x(3), S.h0))));
    end
end
chk.dVerbatim = dv0;

%% ---------------- outputs ----------------
rows = cell(nS*nc*numel(S.procs), 24); r = 0;
for i = 1:nS
    for pr = S.procs
        a = A.(pr);
        for k = 1:nc
            r = r + 1;
            rows(r, :) = {ids(i), grp(i), pr, S.cGrid(k), a.alpha(i, k), a.n(i, k), a.thr(i, k), a.ths(i, k), ...
                1 + S.cGrid(k)/a.n(i, k), a.rmse(i, k), a.rmseW(i, k), a.rmseI(i, k), a.rmseD(i, k), a.sse(i, k), npts(i), ...
                a.thrGtThs(i, k), a.negAny(i, k), a.negPlot(i, k), a.negMeas(i, k), a.neg1cm(i, k), a.negObs(i, k), ...
                a.hNegMax(i, k), a.rmseK(i, k), a.ef(i, k)};
        end
    end
end
TL = cell2table(rows, 'VariableNames', {'Sample_ID', 'Group', 'Procedure', 'c', 'alpha', 'n', 'theta_r', 'theta_s', ...
    'beta', 'RMSE_theta', 'RMSE_theta_wet', 'RMSE_theta_int', 'RMSE_theta_dry', 'SSE', 'Npts', 'theta_r_gt_theta_s', ...
    'C_negative_anywhere', 'C_negative_h_ge_1e_6cm', 'C_negative_h_ge_0p01cm', 'C_negative_h_ge_1cm', ...
    'C_negative_measured_range', 'hNeg_max_cm', 'RMSE_log10K', 'exitflag'});
TL.Procedure = string(TL.Procedure);
TLw = TL; TLw.Procedure = plabel(TLw.Procedure); SUMw = SUM; SUMw.Procedure = plabel(SUMw.Procedure);
LIMw = LIM; LIMw.Properties.VariableNames = cellstr(plabel(string(LIM.Properties.VariableNames)));
DSw = DS; DSw.Properties.VariableNames = cellstr(plabel(string(DS.Properties.VariableNames)));
writetable(TLw,  fullfile(outDir, 'EVM_beta_recheck_scan.csv'));
writetable(SUMw, fullfile(outDir, 'EVM_beta_recheck_summary.csv'));
writetable(LIMw, fullfile(outDir, 'EVM_beta_recheck_limits.csv'));
writetable(DSw,  fullfile(outDir, 'EVM_beta_recheck_development_set.csv'));
writetable(SB,  fullfile(outDir, 'EVM_beta_recheck_soil_beta.csv'));

fid = fopen(fullfile(outDir, 'EVM_beta_recheck_report.txt'), 'w');
if fid < 0, error('Cannot write the report in %s', outDir); end
for f = [1, fid]
    print_report(f, S, SUM, LIMc, crit, LIM, A, isG, kRef, kOpt, bestCount, rangeP2, okP2, SB, DEV, DS, devReg, ...
        devReg2, gReg, chk, ids, npts);
end
fclose(fid);
R = struct('settings', S, 'summary', SUM, 'limitsByProcedure', LIMc, 'limitsBySoil', LIM, 'soilBeta', SB, ...
    'development', DEV, 'devSoils', DS, 'scanTable', TL, 'checks', chk, 'ids', ids);
save(fullfile(outDir, 'EVM_beta_recheck_results.mat'), 'R');
make_figures(outDir, S, SUM, DEV, DS);
fprintf('Outputs written to %s (%.1f min).\n', outDir, toc(tStart)/60);
fprintf('EVM_BETA_RECHECK_DONE\n');
end

%% =====================================================================
%% local functions
%% =====================================================================
function v = mx(x), x = x(isfinite(x)); if isempty(x), v = NaN; else, v = max(x); end, end
function s = plabel(s), s = replace(string(s), ["OB1" "P2"], ["single_start" "constrained"]); end
function v = md(x), x = x(isfinite(x)); if isempty(x), v = NaN; else, v = median(x); end, end
function y = sp(z), y = max(z, 0) + log1p(exp(-abs(z))); end      % log(1 + exp(z)) without overflow

function o = analyse_soil(d, S, opts)
% OB1, P1 and P2 for every c of the grid, their diagnostics, and the soil-specific beta (P2)
nc = numel(S.cGrid);
ws = warning('off', 'all');
P = struct('OB1', nan(nc, 4)); ef = nan(nc, 1);
op = optimoptions('lsqcurvefit', 'Display', 'off');      % the original code: default options
for k = 1:nc
    try
        [x, ~, ~, ef(k)] = lsqcurvefit(orig_model(S.cGrid(k), S), S.x0, d.h, d.th, S.lb0, S.ub0, op);
        P.OB1(k, :) = [x(1) x(3) x(2) x(4)];
    catch
    end
end
[P.P1, P.P2, Q2] = scan_multistart(d.h, d.th, S, opts);
for pr = S.procs
    r = struct();
    for k = 1:nc
        g = diag_fit(P.(pr)(k, :), S.cGrid(k), d, S);
        for fn = fieldnames(g)', r.(fn{1})(k, 1) = g.(fn{1}); end
    end
    r.ef = nan(nc, 1); if pr == "OB1", r.ef = ef; end
    o.(pr) = r;
end
% soil-specific beta: alpha, n, theta_r, theta_s and beta (1-10) fitted together (P2; Script S1)
kRef = find(abs(S.cGrid - S.cRef) < 1e-12); k5 = find(S.cGrid == 5);
Q0 = [repmat(Q2(kRef, :), 5, 1), [1.05; 1.5; 3; 6; 9]; Q2(k5, :), min(1 + 5/P.P2(k5, 2), S.betaMax - 0.01)];
[p, q, sse] = fit_evm(d.h, d.th, NaN, 'P2', Q0, S, opts, true);
o.free = struct('alpha_free', p(1), 'n_free', p(2), 'theta_r_free', p(3), 'theta_s_free', p(4), 'beta_opt', q(5), ...
    'SSE_free', sse, 'RMSE_free', sqrt(sse/numel(d.th)));
warning(ws);
end

function g = diag_fit(p, c, d, S)
% errors, capacity sign and suction range of C < 0, and conductivity error of one fit
fn = {'alpha', 'n', 'thr', 'ths', 'sse', 'rmse', 'rmseW', 'rmseI', 'rmseD', 'thrGtThs', 'negAny', 'negPlot', ...
    'negMeas', 'neg1cm', 'negObs', 'hNegMax', 'rmseK', 'hStar'};
g = cell2struct(num2cell(nan(numel(fn), 1)), fn, 1);
if any(~isfinite(p)), return, end
beta = 1 + c/p(2);
g.alpha = p(1); g.n = p(2); g.thr = p(3); g.ths = p(4);
e = evm_theta(d.h, p, beta, S.h0) - d.th;
g.sse = sum(e.^2); g.rmse = sqrt(mean(e.^2));
cls = 1 + (d.pF >= S.pFedge(1)) + (d.pF >= S.pFedge(2));        % 1 wet, 2 intermediate, 3 dry
rr = nan(1, 3);
for j = 1:3, if any(cls == j), rr(j) = sqrt(mean(e(cls == j).^2)); end, end
g.rmseW = rr(1); g.rmseI = rr(2); g.rmseD = rr(3);
C = evm_capacity(S.hg, p, beta, S.h0); neg = C < 0;
s0 = lead_sign(p, beta, S);
g.thrGtThs = double(p(3) > p(4));
g.negAny = double(any(neg) || s0 < 0);
g.negPlot = double(any(neg & S.hg >= S.hplot));
g.negMeas = double(any(neg & S.hg >= S.hmin));
g.neg1cm = double(any(neg & S.hg >= 1));
g.negObs = double(any(neg & S.hg >= min(d.h)));
if any(neg), g.hNegMax = max(S.hg(neg)); elseif s0 < 0, g.hNegMax = S.hg(1); end
if isfinite(d.Ksat) && ~isempty(d.K)
    g.rmseK = sqrt(mean((evm_logK(d.hK, p, beta, d.Ksat, S) - log10(d.K)).^2));
end
g.hStar = (beta - 1)^(1/beta)/p(1);                 % suction of the maximum non-capillary capacity
end

function f = orig_model(c, S)
% EVM retention of SWRFIT_PDI_EVM.m (line 199) with 0.1 replaced by c, x = [alpha theta_r n theta_s x5];
% for c > S.cVerbatim the same function is evaluated without overflow
if c <= S.cVerbatim
    f = @(vGx,h) (vGx(4)-vGx(2))*((1./(1+(vGx(1)*h).^(vGx(3)))).^(1-1/vGx(3))-(1./(1+(vGx(1)*10^7).^(vGx(3)))).^(1-1/vGx(3)))/(1-(1./(1+(vGx(1)*10^7).^(vGx(3)))).^(1-1/vGx(3)))...
        +vGx(2)*(1-log(1+(vGx(1)*h).^(1+c/vGx(3)))/log(1+(vGx(1)*10^7).^(1+c/vGx(3))));
else
    f = @(vGx,h) evm_theta(h, [vGx(1) vGx(3) vGx(2) vGx(4)], 1 + c/vGx(3), S.h0);
end
end

function [p1, p2, q2] = scan_multistart(h, th, S, opts)
% P1 and P2 for every c (Script S1): 18 starting points plus the solution at the previous c, 12 more
% at c = 0 and c = 0.1 (forward sweep), then a backward sweep; the P2 set lies inside the P1 set, so a
% lower P2 SSE is also a P1 solution and a P1 optimum with theta_r < theta_s is also the P2 optimum
nc = numel(S.cGrid);
kRef = find(abs(S.cGrid - S.cRef) < 1e-12);
Q = struct();
for pr = {'P1', 'P2'}
    Q0 = start_set(pr{1}, th); Qx = start_extra(pr{1}, th);
    qAll = NaN(nc, 4); sAll = Inf(nc, 1); qPrev = [];
    for k = 1:nc
        Qk = [qPrev; Q0];
        if k == 1 || k == kRef, Qk = [Qk; Qx]; end %#ok<AGROW>
        [~, qAll(k, :), sAll(k)] = fit_evm(h, th, S.cGrid(k), pr{1}, Qk, S, opts, false);
        qPrev = qAll(k, :);
    end
    for k = nc-1:-1:1
        [~, q, sse] = fit_evm(h, th, S.cGrid(k), pr{1}, qAll(k+1, :), S, opts, false);
        if sse < sAll(k), qAll(k, :) = q; sAll(k) = sse; end
    end
    Q.(pr{1}).q = qAll; Q.(pr{1}).sse = sAll;
end
for k = 1:nc
    a = q2p(Q.P1.q(k, :), 'P1'); b = q2p(Q.P2.q(k, :), 'P2');
    if Q.P2.sse(k) < Q.P1.sse(k)*(1 - 1e-12)
        Q.P1.q(k, :) = [log10(b(1)) b(2) b(3) b(4)]; Q.P1.sse(k) = Q.P2.sse(k);
    elseif a(3) <= S.fMax*a(4) && a(4) >= 0.01 && Q.P1.sse(k) < Q.P2.sse(k)*(1 - 1e-12)
        Q.P2.q(k, :) = [log10(a(1)) a(2) a(4) a(3)/a(4)]; Q.P2.sse(k) = Q.P1.sse(k);
    end
end
p1 = NaN(nc, 4); p2 = NaN(nc, 4);
for k = 1:nc, p1(k, :) = q2p(Q.P1.q(k, :), 'P1'); p2(k, :) = q2p(Q.P2.q(k, :), 'P2'); end
q2 = Q.P2.q;
end

function [p, qBest, sseBest] = fit_evm(h, th, c, pr, Q0, S, opts, free)
% Eq. (22) with lsqcurvefit; q = [log10(alpha) n theta_r theta_s] (P1) or [log10(alpha) n theta_s f]
% (P2), plus beta as q(5) when free
if strcmp(pr, 'P1'), lb = [-5 S.nMin 0 0];    ub = [0 10 1 1];
else,                lb = [-5 S.nMin 0.01 0]; ub = [0 10 1 S.fMax];
end
if free
    lb = [lb 1]; ub = [ub S.betaMax];
    fun = @(q, hh) evm_theta(hh, q2p(q, pr), q(5), S.h0);
else
    fun = @(q, hh) evm_theta(hh, q2p(q, pr), 1 + c/q(2), S.h0);
end
sseBest = Inf; qBest = Q0(1, :);
for j = 1:size(Q0, 1)
    q0 = min(max(Q0(j, :), lb + 1e-9), ub - 1e-9);
    try
        [q, sse] = lsqcurvefit(fun, q0, h, th, lb, ub, opts);
    catch
        continue
    end
    if sse < sseBest, sseBest = sse; qBest = q; end
end
p = q2p(qBest, pr);
end

function Q0 = start_set(pr, th)
ths0 = min(1.02*max(th), 0.99);
[la, nn, rr] = ndgrid([-3 -2 -1.3], [1.15 1.4 2.0], [0.05 0.3]);
if strcmp(pr, 'P1'), Q0 = [la(:) nn(:) rr(:) repmat(ths0, numel(la), 1)];
else,                Q0 = [la(:) nn(:) repmat(ths0, numel(la), 1) rr(:)];
end
end

function Qx = start_extra(pr, th)
% 12 extra starts with n -> 1, large n and theta_r (P1) or f (P2) -> 1
ths0 = min(1.02*max(th), 0.99);
[la, nn, rr] = ndgrid([-2.5 -1], [1.02 4], [0.01 0.6 0.97]);
if strcmp(pr, 'P1'), Qx = [la(:) nn(:) rr(:) repmat(ths0, numel(la), 1)];
else,                Qx = [la(:) nn(:) repmat(ths0, numel(la), 1) rr(:)];
end
end

function p = q2p(q, pr)
if strcmp(pr, 'P1'), p = [10^q(1) q(2) q(3) q(4)];
else,                p = [10^q(1) q(2) q(4)*q(3) q(3)];
end
end

function th = evm_theta(h, p, beta, h0)
% Eq. (7): theta = (ths - thr) Sc + thr Snc, with Eqs. (3), (4) and (6)
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n;
G  = (1 + (al*h).^n).^(-m);
G0 = (1 + (al*h0)^n)^(-m);
Sc = (G - G0)/(1 - G0);
th = (ths - thr)*Sc + thr*evm_snc(al*h, al*h0, beta);
end

function Snc = evm_snc(x, x0, beta)
% Eq. (6), Snc = 1 - ln(1 + x^beta)/ln(1 + x0^beta) with x = alpha h; for large beta as softplus (no
% overflow); beta = Inf: the limit 1 - max(ln x, 0)/ln x0 (broken line with a kink at h = 1/alpha)
if isinf(beta)
    if x0 <= 1, Snc = ones(size(x)); else, Snc = 1 - max(log(x), 0)/log(x0); end
elseif beta*log(x0) < 700
    Snc = 1 - log1p(x.^beta)/log1p(x0^beta);
else
    Snc = 1 - sp(beta*log(x))/sp(beta*log(x0));
end
end

function C = evm_capacity(h, p, beta, h0)
% Eq. (9) with h > 0: C = (ths - thr)(-dSc/dh) + thr(-dSnc/dh), Eqs. (10)-(11)
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n; x = al*h;
G0  = (1 + (al*h0)^n)^(-m);
dSc = m*n*al*x.^(n - 1)./(1 + x.^n).^(m + 1)/(1 - G0);
if isinf(beta)
    dSnc = (x > 1)./(h*log(al*h0));
else
    dSnc = beta./h./(1 + exp(-beta*log(x)))/sp(beta*log(al*h0));   % = beta al x^(beta-1)/((1 + x^beta) ln(1 + x0^beta))
end
C = (ths - thr)*dSc + thr*dSnc;
end

function s0 = lead_sign(p, beta, S)
% sign of C as h -> 0+: C ~ A (alpha h)^(n-1) + B (alpha h)^(beta-1), the smaller exponent dominates
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n;
G0 = (1 + (al*S.h0)^n)^(-m);
A = (ths - thr)*m*n/(1 - G0); ec = n - 1;
if isinf(beta), s0 = sign(A); return, end            % Snc = 1 for alpha h < 1
B = thr*beta/sp(beta*log(al*S.h0)); enc = beta - 1;
if A >= 0,        s0 = sign(A + B);
elseif enc < ec,  s0 = sign(B);
elseif enc > ec,  s0 = sign(A);
else,             s0 = sign(A + B);
end
end

function lK = evm_logK(h, p, beta, Ksc, S)
% log10 K (cm/d) = log10(Kc + Knc + Kv): Eq. (14) with Ks,c counted once (capillary term in logarithms,
% no underflow as n -> 1), Eq. (20) with l = 0, Ks,nc of Eq. (17), Kv of Appendix A at 20 C
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n;
lG  = -m*log1p((al*h).^n);
lG0 = -m*log1p((al*S.h0)^n);
lSc = lG + log1p(-exp(lG0 - lG)) - log1p(-exp(lG0));
le = lG/m; ld = lG0/m;
u  = (exp(le) - exp(ld))/(1 - exp(ld));
lB = log(-expm1(m*log1p(-u)));
t  = u < 1e-10;
lB(t) = log(m) + le(t) + log1p(-exp(ld - le(t))) - log1p(-exp(ld));
Kc = exp(log(Ksc) + 0.5*lSc + 2*lB);
G = exp(lG); G0 = exp(lG0);
Sc = max((G - G0)/(1 - G0), 0);
Snc = max(evm_snc(al*h, al*S.h0, beta), 0);
Ksnc = 86.2*evm_theta(1e5, p, beta, S.h0)*al^1.5;                   % Eq. (17), ha = 1/alpha
Rv = al*S.h0;
Knc = Ksnc*Rv.^(-1.5*(1 - Snc));                                    % Eq. (20), l = 0
T = 293.15; Mw = 0.018015; g = 9.81; Rg = 8.314;                     % Appendix A, 20 C
rho_sv = 1e-3*exp(31.3716 - 6014.79/T - 7.92495e-3*T)/T;
rho_w  = 998.2;
Da  = 2.12e-5*(T/273.15)^2*1e4*86400;
cH  = Mw*g/(Rg*T)/100;
tau = max(ths - ((ths - thr)*Sc + thr*Snc), 0).^(10/3)/ths^2;
Kv  = rho_sv/rho_w*tau*Da*cH.*exp(-cH*h);
lK = log10(Kc + Knc + Kv);
end

function st = ols_line(x, y)
% ordinary least squares y = a + b x, plus the through-origin fit y = b0 x (Script S1)
x = x(:); y = y(:); N = numel(x);
st = struct('N', N, 'a', NaN, 'b', NaN, 'ciA', NaN, 'ciB', NaN, 'pB', NaN, 'R2', NaN, 'b0', NaN, 'R2_0', NaN);
if N < 4 || std(x) == 0 || std(y) == 0, return, end
X = [ones(N, 1) x]; b = X\y; res = y - X*b; dfree = N - 2;
covb = (sum(res.^2)/dfree)*inv(X'*X); se = sqrt(diag(covb)); %#ok<MINV>
t = b./se; pv = betainc(dfree./(dfree + t.^2), dfree/2, 0.5);
if exist('tinv', 'file'), tq = tinv(0.975, dfree); else, tq = 1.96; end
st.a = b(1); st.b = b(2); st.ciA = tq*se(1); st.ciB = tq*se(2); st.pB = pv(2);
st.R2 = 1 - sum(res.^2)/sum((y - mean(y)).^2);
st.b0 = (x'*y)/(x'*x);
st.R2_0 = 1 - sum((y - st.b0*x).^2)/sum((y - mean(y)).^2);
end

function s = cstr(c)
if isinf(c) && c > 0, s = 'Inf'; elseif isinf(c), s = 'none'; else, s = sprintf('%.2f', c); end
end

function print_report(f, S, SUM, LIMc, crit, LIM, A, isG, kRef, kOpt, bestCount, rangeP2, okP2, SB, DEV, DS, ...
    devReg, devReg2, gReg, chk, ids, npts)
nG = sum(isG);
fprintf(f, '\n========== Closure beta = 1 + c/n: single-start and constrained fits (%d German soils) ==========\n', nG);
fprintf(f, 'single_start = independent bounds, one lsqcurvefit run from x0 = [0.01 0.001 1.5 0.01] (SWRFIT_PDI_EVM.m), n >= 1\n');
fprintf(f, 'constrained  = theta_s > theta_r >= 0 enforced, multi-start + continuation in c\n');
fprintf(f, 'RMSE_K: Eq. (20), l = 0, Ks,c = measured KSAT (counted once), unsaturated K data (Script S6)\n');
fprintf(f, 'C < 0 checked on 1e-12 cm <= h < h0 and as h -> 0+ ("anywhere"), for h >= 1e-6 cm (Figs. 2-3), h >= 0.01 cm, h >= 1 cm\n');
fprintf(f, '\nConsistency: single start at c = 0.1 vs Script S7 (%d soils): max |dRMSE_theta| = %.2e, max |d parameter| = %.2e\n', ...
    chk.nS7, chk.dS7_rmse, chk.dS7_par);
fprintf(f, '             constrained at c <= 5 vs Script S1 (%d fits): max |dRMSE_theta| = %.2e (99th percentile %.2e)\n', ...
    chk.nS1, chk.dS1_rmse, chk.dS1_q);
fprintf(f, '             soil-specific beta vs Script S1 (%d soils): max |dRMSE_free| = %.2e\n', chk.nS1b, chk.dS1_free);
fprintf(f, '             constrained at c = 0.1 vs Script S6 (%d soils): max |dRMSE_theta| = %.2e (%d soils > 1e-6), max |dRMSE_K| = %.2e\n', ...
    chk.nS6, chk.dS6_rmse, chk.nS6diff, chk.dS6_K);
fprintf(f, '             expression of the original code vs overflow-free function (c <= %g): max |d theta| = %.2e\n', ...
    S.cVerbatim, chk.dVerbatim);
for pr = S.procs
    q = SUM(SUM.Procedure == pr, :);
    fprintf(f, '\n--- %s (mean of per-soil RMSE; RMSE_K over %d soils with KSAT) ---\n', plabel(pr), q.NsoilsK(1));
    fprintf(f, '%6s %9s %7s %9s %8s %8s %8s %6s %6s %6s %6s %6s %6s %8s %8s %7s %5s %7s %8s\n', 'c', 'RMSEth', 'd%', ...
        'pooled', 'wet', 'int', 'dry', 'thr>s', 'C<0', 'h>1e-6', 'h>.01', 'h>1', 'meas', 'hNeg.01', 'RMSEK', 'dK%', 'ef0', ...
        'h*med', 'hNegAny');
    for k = 1:height(q)
        fprintf(f, '%6s %9.6f %7.3f %9.6f %8.5f %8.5f %8.5f %6d %6d %6d %6d %6d %6d %8.3g %8.4f %7.2f %5d %7.3g %8.2g\n', ...
            cstr(q.c(k)), q.meanRMSE_theta(k), q.dMeanRMSE_theta_pct(k), q.pooledRMSE_theta(k), q.meanRMSE_theta_wet(k), ...
            q.meanRMSE_theta_int(k), q.meanRMSE_theta_dry(k), q.N_theta_r_gt_theta_s(k), q.N_C_negative_anywhere(k), ...
            q.N_C_negative_h_ge_1e_6cm(k), q.N_C_negative_h_ge_0p01cm(k), q.N_C_negative_h_ge_1cm(k), ...
            q.N_C_negative_measured_range(k), q.max_hNeg_cm(k), q.meanRMSE_log10K(k), q.dMeanRMSE_log10K_pct(k), ...
            q.N_exitflag0(k), q.median_hStar_cm_thr_ge_0p05(k), q.max_hNeg_any_cm(k));
    end
end
fprintf(f, '\nLargest admissible c (C > 0 in every soil for all c'' <= c):\n');
fprintf(f, '   %-12s %12s %12s %12s %12s\n', '', 'anywhere', 'h>=1e-6 cm', 'h>=0.01 cm', 'h>=1 cm');
for ip = 1:numel(S.procs)
    fprintf(f, '   %-12s %12s %12s %12s %12s\n', plabel(S.procs(ip)), cstr(LIMc(ip, 1)), cstr(LIMc(ip, 2)), cstr(LIMc(ip, 3)), cstr(LIMc(ip, 4)));
end
fprintf(f, '\nOnset of C < 0 with independent bounds (soils with C < 0 anywhere; hneg = largest suction with C < 0):\n');
for pr = "OB1"
    a = A.(pr);
    for k = find(S.cGrid >= 0.08 & S.cGrid <= 0.2)
        i = find(isG & a.negAny(:, k) == 1);
        fprintf(f, '   %s c = %s: %d soils', plabel(pr), cstr(S.cGrid(k)), numel(i));
        for j = i(:)'
            fprintf(f, ' | %s n = %.4f thr = %.3f ths = %.3f hneg = %.2g cm', ids(j), a.n(j, k), a.thr(j, k), ...
                a.ths(j, k), a.hNegMax(j, k));
        end
        fprintf(f, '\n');
    end
end
for pr = "OB1"
    cm = LIM.(char("cmax_" + pr)); lim = find(isG & isfinite(cm));
    [~, o] = sort(cm(lim));
    fprintf(f, '\nSoils whose admissible limit (h >= 0.01 cm) lies inside the grid, %s: %d (fit at the first violating c)\n', ...
        plabel(pr), numel(lim));
    for j = o(:)'
        i = lim(j);
        fprintf(f, '   %-9s c_max = %5s | c = %5s: n = %.4f theta_r = %.3f theta_s = %.3f hneg = %.3g cm RMSE = %.5f\n', ids(i), ...
            cstr(cm(i)), cstr(LIM.(char("first_violating_c_" + pr))(i)), LIM.(char("n_" + pr))(i), ...
            LIM.(char("theta_r_" + pr))(i), LIM.(char("theta_s_" + pr))(i), LIM.(char("hNeg_max_cm_" + pr))(i), ...
            LIM.(char("RMSE_theta_" + pr))(i));
    end
end
q = SUM(SUM.Procedure == "P2", :);
fprintf(f, '\nconstrained: pooled-SSE optimum on the grid: c = %s (mean RMSE_theta %.6f, %.2f%% vs c = 0.1; RMSE_K %.4f, %.2f%%)\n', ...
    cstr(S.cGrid(kOpt)), q.meanRMSE_theta(kOpt), q.dMeanRMSE_theta_pct(kOpt), q.meanRMSE_log10K(kOpt), q.dMeanRMSE_log10K_pct(kOpt));
fprintf(f, 'constrained: best grid c of each soil (%d soils): ', sum(okP2));
for k = 1:numel(S.cGrid), if bestCount(k) > 0, fprintf(f, '%s: %d  ', cstr(S.cGrid(k)), bestCount(k)); end, end
fprintf(f, '\nconstrained: soils whose RMSE_theta changes by < 1e-4 over the whole grid: %d of %d (< 1e-5: %d)\n', ...
    sum(rangeP2 < 1e-4), numel(rangeP2), sum(rangeP2 < 1e-5));
dM = A.P2.rmse(okP2, :) - A.P2.rmse(okP2, kRef);
fprintf(f, 'constrained: per-soil RMSE_theta change from c = 0.1 to c = Inf: median %.2e, 90th pct %.2e, max %.2e (decrease); increase in %d soils\n', ...
    median(-dM(:, end)), quantile(-dM(:, end), 0.9), max(-dM(:, end)), sum(dM(:, end) > 1e-9));
% -dSnc/dln h = beta sigma(beta ln(alpha h))/ln(1 + (alpha h0)^beta) rises from 10% to 90% of its plateau over
% 2 ln(9)/beta in ln h, i.e. 2 ln(9)/(beta ln 10) decades around h = 1/alpha (PDI: 2 ln(9) b decades)
fprintf(f, 'Width of the capillary/non-capillary transition (10-90%% of -dSnc/dln h) = 2 ln9/(beta ln10) decades, constrained fits:\n');
for c = [0 0.1 1 5 10 50 100]
    k = find(S.cGrid == c); w = 2*log(9)./((1 + c./A.P2.n(okP2, k))*log(10));
    fprintf(f, '   c = %s: median %.3f decades (range %.3f-%.3f)\n', cstr(c), median(w), min(w), max(w));
end
fprintf(f, '   c = Inf: 0 (kink at h = 1/alpha) | PDI with b = 0.1-0.3: %.2f-%.2f decades\n', 2*log(9)*0.1, 2*log(9)*0.3);
G = isG;
fprintf(f, '\nSoil-specific beta (constrained, 1 <= beta <= %g), %d German soils:\n', S.betaMax, nG);
e = G & ~SB.beta_no_effect;
fprintf(f, '   theta_r < 0.02 (beta has no effect): %d; others %d: at beta = 1: %d, at beta = %g: %d, interior: %d\n', ...
    sum(G & SB.beta_no_effect), sum(e), sum(e & SB.atLowerBound), S.betaMax, sum(e & SB.atUpperBound), ...
    sum(e & ~SB.atLowerBound & ~SB.atUpperBound));
fprintf(f, '   beta* quantiles (0 10 25 50 75 90 100%%, theta_r >= 0.02): %s\n', mat2str(quantile(SB.beta_opt(e), [0 .1 .25 .5 .75 .9 1]), 4));
fprintf(f, '   RMSE gain of the free beta over the closure: median %.2f%%, 90th pct %.2f%% (all soils)\n', ...
    median(SB.gain_pct(G)), quantile(SB.gain_pct(G), 0.9));
fprintf(f, '   BIC (one more parameter): free beta supported (dBIC > 2) in %d of %d soils, strongly (dBIC > 6) in %d, (dBIC > 10) in %d\n', ...
    sum(SB.dBIC(G) > 2), nG, sum(SB.dBIC(G) > 6), sum(SB.dBIC(G) > 10));
fprintf(f, '   beta-1 = %.3f + %.3f/n, R2 = %.3f (N = %d, all German soils)\n', gReg.a, gReg.b, gReg.R2, gReg.N);
fprintf(f, '\n--- Development set (%d Table 1 soils with data): pooled RMSE_theta by c (thr> = theta_r > theta_s, C<0 anywhere) ---\n', ...
    height(DS));
fprintf(f, '%6s', 'c'); for pr = S.procs, fprintf(f, ' %12s %4s %4s', plabel(pr), 'thr>', 'C<0'); end, fprintf(f, '\n');
for k = 1:numel(S.cGrid)
    fprintf(f, '%6s', cstr(S.cGrid(k)));
    for pr = S.procs, fprintf(f, ' %12.6f %4d %4d', DEV.pooled.(pr)(k), DEV.thrGt.(pr)(k), DEV.negAny.(pr)(k)); end
    fprintf(f, '\n');
end
fprintf(f, 'Optimum of the pooled SSE: ');
for pr = S.procs, fprintf(f, '%s c = %s   ', plabel(pr), cstr(S.cGrid(DEV.kOpt.(pr)))); end
fprintf(f, '\nPer soil (c = 0.1): single start n, theta_r, theta_s, RMSE | constrained n, theta_r, RMSE | beta*, c* = n(beta*-1), dBIC, cmax single start\n');
for i = 1:height(DS)
    fprintf(f, '   %-12s %.3f %.3f %.3f %.5f | %.3f %.3f %.5f | beta* %6.3f c* %7.3f dBIC %6.1f | %5s\n', DS.Sample_ID(i), ...
        DS.n_OB1(i), DS.theta_r_OB1(i), DS.theta_s_OB1(i), DS.RMSE_theta_OB1(i), DS.n_P2(i), DS.theta_r_P2(i), ...
        DS.RMSE_theta_P2(i), DS.beta_opt(i), DS.c_free(i), DS.dBIC(i), cstr(DS.cmax_OB1(i)));
end
fprintf(f, 'Development set: beta*-1 = %.3f(+-%.3f) + %.3f(+-%.3f)/n, R2 = %.3f, p = %.2g (N = %d)\n', devReg.a, devReg.ciA, ...
    devReg.b, devReg.ciB, devReg.R2, devReg.pB, devReg.N);
fprintf(f, '                 theta_r >= 0.02: beta*-1 = %.3f + %.3f/n, R2 = %.3f, p = %.2g (N = %d)\n', devReg2.a, devReg2.b, ...
    devReg2.R2, devReg2.pB, devReg2.N);
fprintf(f, 'Retention points: %d German soils, %d in total\n', sum(npts(G)), sum(npts));
end

function make_figures(outDir, S, SUM, DEV, DS)
cx = S.cGrid; cx(cx == 0) = 0.004; cx(isinf(cx)) = 300;
tk = [0.004 0.01 0.1 1 10 100 300]; tl = {'0', '0.01', '0.1', '1', '10', '100', '\infty'};
sty.OB1 = {[0.85 0.33 0.10], 'o-', 'independent bounds, single start'};
sty.P2 = {[0.47 0.67 0.19], 'd-', '\theta_s > \theta_r \geq 0'};
lab = arrayfun(@(p) sty.(p){3}, S.procs, 'UniformOutput', false);
ib = S.procs(S.procs ~= "P2");                     % procedures with independent bounds
f1 = figure('Visible', 'off', 'Position', [60 60 1560 900]);
tiledlayout(f1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% (a) retention error
nexttile; hold on
for pr = S.procs
    q = SUM(SUM.Procedure == pr, :);
    plot(cx, q.meanRMSE_theta, sty.(pr){2}, 'Color', sty.(pr){1}, 'LineWidth', 1.3, 'MarkerSize', 4);
end
set(gca, 'XScale', 'log', 'XTick', tk, 'XTickLabel', tl); xlim([0.003 400]); grid on; box on
xline(S.cRef, 'k--', 'HandleVisibility', 'off'); xlabel('c in \beta = 1 + c/n'); ylabel('mean RMSE_\theta (cm^3 cm^{-3})');
legend(lab, 'Location', 'southwest'); title(sprintf('(a) retention error, %d soils', SUM.Nfits(1)));
% (b) onset of negative capacity with independent bounds
nexttile; hold on
ls = {':', '--', '-'}; crit = {'N_C_negative_anywhere', 'N_C_negative_h_ge_1e_6cm', 'N_C_negative_h_ge_0p01cm'};
for pr = ib
    q = SUM(SUM.Procedure == pr & SUM.c <= 0.3, :);
    for j = 1:3
        plot(q.c, q.(crit{j}), [ls{j} 'o'], 'Color', sty.(pr){1}, 'LineWidth', 1.3, 'MarkerSize', 3.5, ...
            'MarkerFaceColor', sty.(pr){1}, 'HandleVisibility', 'off');
    end
end
q = SUM(SUM.Procedure == "P2" & SUM.c <= 0.3, :);
plot(q.c, q.N_C_negative_anywhere, sty.P2{2}, 'Color', sty.P2{1}, 'LineWidth', 1.3, 'MarkerSize', 3.5, 'HandleVisibility', 'off');
xline(S.cRef, 'k--', 'HandleVisibility', 'off'); xlim([0 0.305]); grid on; box on
xlabel('c in \beta = 1 + c/n'); ylabel('soils with C(h) < 0');
h1 = plot(nan, nan, 'k:', 'LineWidth', 1.3); h2 = plot(nan, nan, 'k--', 'LineWidth', 1.3); h3 = plot(nan, nan, 'k-', 'LineWidth', 1.3);
legend([h1 h2 h3], {'anywhere in 0 < h < h_0', 'h \geq 10^{-6} cm (range of Figs. 2-3)', 'h \geq 0.01 cm (measured range)'}, ...
    'Location', 'northwest');
title('(b) onset of C < 0 (colours as in a)');
% (c) conductivity error
nexttile; hold on
for pr = S.procs
    q = SUM(SUM.Procedure == pr, :);
    plot(cx, q.meanRMSE_log10K, sty.(pr){2}, 'Color', sty.(pr){1}, 'LineWidth', 1.3, 'MarkerSize', 4);
end
set(gca, 'XScale', 'log', 'XTick', tk, 'XTickLabel', tl); xlim([0.003 400]); grid on; box on
xline(S.cRef, 'k--'); xlabel('c in \beta = 1 + c/n'); ylabel('mean RMSE of log_{10} K');
title(sprintf('(c) predicted conductivity, %d soils with K_{sat}', SUM.NsoilsK(1)));
% (d) suction range of negative capacity with independent bounds
nexttile; hold on
for pr = ib
    q = SUM(SUM.Procedure == pr, :);
    y = q.max_hNeg_any_cm; y(~(y > 0)) = NaN;
    plot(cx, y, sty.(pr){2}, 'Color', sty.(pr){1}, 'LineWidth', 1.3, 'MarkerSize', 4);
end
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', tk, 'XTickLabel', tl); xlim([0.003 400]); ylim([1e-13 1e3]);
grid on; box on
yline(S.hplot, 'k--', 'h = 10^{-6} cm (lower end of Figs. 2-3)', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
yline(S.hmin, 'k-', 'h = 0.01 cm (wettest data)', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xline(S.cRef, 'k--', 'HandleVisibility', 'off');
xlabel('c in \beta = 1 + c/n'); ylabel('largest suction with C < 0, all soils (cm)');
legend(arrayfun(@(p) sty.(p){3}, ib, 'UniformOutput', false), 'Location', 'southeast');
title('(d) suction range of C < 0 with independent bounds');
exportgraphics(f1, fullfile(outDir, 'FigS3_closure_sensitivity.png'), 'Resolution', 250); close(f1);

f2 = figure('Visible', 'off', 'Position', [60 60 1400 520]);
tiledlayout(f2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on
for pr = S.procs
    plot(cx, DEV.pooled.(pr), sty.(pr){2}, 'Color', sty.(pr){1}, 'LineWidth', 1.3, 'MarkerSize', 4);
end
set(gca, 'XScale', 'log', 'XTick', tk, 'XTickLabel', tl); xlim([0.003 400]); grid on; box on
xline(S.cRef, 'k--', 'HandleVisibility', 'off'); xlabel('c in \beta = 1 + c/n'); ylabel('pooled RMSE_\theta (cm^3 cm^{-3})');
legend(lab, 'Location', 'northwest');
title(sprintf('(a) the %d Table 1 soils with data', height(DS)));
nexttile; hold on
xi = linspace(0.1, 1.05, 100); plot(xi, 1 + 0.1*xi, 'r-', 'LineWidth', 2);
x = 1./DS.n_free; y = DS.beta_opt; ly = y;
nm = replace(strrep(DS.Sample_ID, "_", "-"), ["Berlin-Sand" "Gilat-Loam" "Pachapa-FSC" "Pachepsky-SL" "Pachepsky-CL"], ...
    ["Berlin sand" "Gilat loam" "Pachapa FSC" "Pachepsky SL" "Pachepsky CL"]);
scatter(x, y, 40, DS.theta_r_free, 'filled');
top = find(y >= 9.9); [~, o] = sort(x(top)); ly(top(o)) = 10*0.88.^(0:numel(o) - 1);   % labels at the bounds
bot = find(y <= 1.05); [~, o] = sort(x(bot)); ly(bot(o)) = 1.12.^(1:numel(o));          % are staggered
for i = 1:height(DS)
    if ly(i) ~= y(i), plot([x(i) x(i) + 0.012], [y(i) ly(i)], '-', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off'); end
    text(x(i) + 0.014, ly(i), char(nm(i)), 'FontSize', 8);
end
set(gca, 'YScale', 'log', 'YTick', [1 1.1 1.2 1.5 2 3 5 10]); ylim([0.95 11.5]); xlim([0.15 0.85]); grid on; box on
cb = colorbar; cb.Label.String = '\theta_r (cm^3 cm^{-3})';
xlabel('1/n'); ylabel('soil-specific \beta^*'); legend('\beta = 1 + 0.1/n', '\beta^*', 'Location', 'northwest');
title('(b) soil-specific \beta^* (\theta_s > \theta_r \geq 0)');
exportgraphics(f2, fullfile(outDir, 'FigS4_Table1_soils.png'), 'Resolution', 250); close(f2);
end
