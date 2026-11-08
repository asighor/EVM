function R = EVM_original_code_recheck(dataDir, nSoilsMax)
%EVM_ORIGINAL_CODE_RECHECK  The previous fitting procedure, reproduced with the fitting code of the previous
% version of the manuscript (SWRFIT_PDI_EVM.m), compared with the fits under theta_s > theta_r >= 0
% (Texts S3, S6 and S7 of the Supporting Information; Table 6 and Figures 2 and 3).
%
%   R = EVM_original_code_recheck              % data files next to this file
%   R = EVM_original_code_recheck(dataDir)
%   R = EVM_original_code_recheck(dataDir, 30) % quick test on the first 30 soils (or a list of IDs)
%   R = EVM_original_code_recheck(dataDir, 'figures')   % redraw the figures from the saved results
%
% Previous procedure: ONE lsqcurvefit run per soil (default options) from
%   x0 = [alpha theta_r n theta_s x5] = [0.01 0.001 1.5 0.01 0.01]   (x5 is not used by the models)
% with 0 <= alpha, theta_r, theta_s <= 1 and 1 <= n <= 10; theta_s > theta_r is not imposed. It reproduces
% the published Tables 2 and 3 and Figure 2 of the previous version. Constrained procedure:
% theta_s > theta_r >= 0 enforced, theta_r = f theta_s with f <= 0.9999, 30-110 starting points (Scripts S2, S6).
% For every fit: theta_r > theta_s, the sign of C on 1e-12 cm <= h < h0, on the plotted range
% 1e-6 <= h <= 1e7 cm, at h >= 0.01 cm and inside the measured range, the sign of C as h -> 0+,
% C(1e-6 cm), and RMSE_theta and RMSE_K computed exactly as in the previous code (Table 3: Ks = the
% first conductivity row of the soil, every conductivity row compared, including KSAT; for EVM also
% with Knc = Ks,nc Snc^0.5 (alpha h0)^(-1.5(1 - Snc)), Snc with beta, i.e. l = 0.5 in Eq. 20, lines
% 402-404 and 448 of the previous code). The capacity is computed with the formulas C_pdi and C_evm of
% the previous code and checked against the implementation of Scripts S2-S6.
%
% Outputs (folder recheck_results next to the data): EVM_recheck_report.txt, EVM_recheck_fits.csv
% (every fit: Procedure "previous" or "constrained"), EVM_recheck_summary.csv, EVM_recheck_table3.csv
% (Table 3 as computed by the previous code), EVM_recheck_table6.csv (Table 6),
% Fig2_previous_vs_constrained_part1.png, Fig2_previous_vs_constrained_part2.png and
% Fig3_TUBS205_previous_vs_constrained.png (Figures 2 and 3), EVM_recheck_results.mat.
% Runtime: about 5 min with 24 parallel workers (R2025b, Optimization and Parallel Computing Toolboxes).

if nargin < 1 || isempty(dataDir)
    dataDir = fileparts(mfilename('fullpath'));
end
outDir = fullfile(dataDir, 'recheck_results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
tStart = tic;
if nargin > 1 && (ischar(nSoilsMax) || isstring(nSoilsMax)) && any(strcmpi(nSoilsMax, 'figures'))
    L = load(fullfile(outDir, 'EVM_recheck_results.mat'), 'R');       % redraw the figures only
    make_figures(outDir, L.R.settings, L.R.fits, L.R.data, L.R.ids, L.R.bestVariant);
    R = L.R; fprintf('Figures redrawn in %s.\n', outDir);
    return
end

%% ---------------- settings ----------------
S.h0 = 1e7; S.c = 0.1; S.fMax = 0.9999; S.nMin = 1.0001; S.hmin = 1e-2;
S.hg = logspace(-12, log10(S.h0) - 1e-9, 20000)';     % 1e-12 cm <= h < h0
S.hp = logspace(-6, 7, 20001)';                        % range plotted in Figures 2 and 3
S.groups = ["TUDD" "UFZ" "FSUJ" "KIT" "TUBS"];
S.procs = ["OB1" "P2"];          % OB1 = previous procedure (published fits); P2 = theta_s > theta_r >= 0
S.models = ["PDI" "EVM"];
xB = [0.01; 0.001; 1.5; 0.01; 0.01];                    % x0 of the previous procedure (line 210 of SWRFIT_PDI_EVM.m)
S.V = struct('tag', {'OB1'}, 'x0', {xB}, 'lbE', {[0; 0; 1; 0; 0]}, 'fitPDI', {true});
S.ms3.PDI = [0.0096 0.0050 0.0085 0.0061 0.0048; 1.40 0.56 0.89 1.00 0.87];   % published Table 3
S.ms3.EVM = [0.0096 0.0051 0.0085 0.0061 0.0048; 1.43 0.62 0.95 1.07 0.94];
S.T1ids = ["FSUJ_023" "FSUJ_035" "TUBS_027" "TUBS_128" "KIT_022" "KIT_009" "TUBS_106"];
% published Table 2 (theta_r theta_s alpha n) of these soils, PDI and EVM
S.ms2.PDI = [0.171 0.460 0.019 1.726; 9e-5 0.589 0.003 1.378; 0.081 0.400 0.003 1.291; 0.245 0.383 0.015 1.464;
    0.301 0.550 0.080 1.662; 0.278 0.658 0.045 1.401; 0.486 0.601 0.020 1.680];
S.ms2.EVM = [0.173 0.461 0.017 1.765; 5e-5 0.589 0.003 1.378; 0.097 0.400 0.003 1.330; 0.249 0.383 0.011 1.581;
    0.303 0.551 0.070 1.725; 0.304 0.658 0.038 1.478; 0.485 0.603 0.014 1.995];
S.ms2K.PDI = [0.35 0.93 0.22 0.28 0.27 0.42 0.38];       % published Table 2 RMSE_K of these soils
S.ms2K.EVM = [0.28 0.93 0.11 0.48 0.33 0.53 0.37];
S.two = ["FSUJ_015" "TUBS_205"];
opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxFunctionEvaluations', 4000, 'MaxIterations', 1000, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-10);
optsLC = optimoptions('lsqcurvefit', 'Algorithm', 'interior-point', 'Display', 'off', ...
    'MaxFunctionEvaluations', 8000, 'MaxIterations', 2000, 'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-10, 'ConstraintTolerance', 1e-9);

%% ---------------- data (read as in the original code) ----------------
ret  = readtable(fullfile(dataDir, 'RetMeas.csv'),  'TextType', 'string');
cond = readtable(fullfile(dataDir, 'CondMeas.csv'), 'TextType', 'string');
ids = unique(ret.Sample_ID, 'stable');
if nargin > 1 && ~isempty(nSoilsMax)
    if isnumeric(nSoilsMax), ids = ids(1:min(nSoilsMax, numel(ids)));
    else,                    ids = ids(ismember(ids, string(nSoilsMax)));
    end
    ids = unique([ids; intersect([S.two S.T1ids]', unique(ret.Sample_ID), 'stable')], 'stable');
end
nS = numel(ids);
[~, ir] = ismember(ret.Sample_ID, ids);
[~, ic] = ismember(cond.Sample_ID, ids);
D = struct('id', cellstr(ids), 'pF', [], 'th', [], 'hm', [], 'Km', [], 'hasKsat', false);
for i = 1:nS
    s = ir == i; k = ic == i;
    D(i).pF = ret.pF(s); D(i).th = ret.theta(s);
    D(i).hm = 10.^cond.pF(k); D(i).Km = cond.k(k);          % every row, file order (Km(1) = Ks)
    D(i).hasKsat = any(k) && cond.MeasType(find(k, 1)) == "KSAT";
end
grp = extractBefore(ids, "_");
fprintf('Loaded %d soils (%d with KSAT as the first conductivity row).\n', nS, sum([D.hasKsat]));

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
T = struct2table(cell2mat(cellfun(@(o) o, out, 'UniformOutput', false)));
T = [table(ids, grp, [D.hasKsat]', arrayfun(@(d) 10^min(d.pF), D), 'VariableNames', ...
    {'Sample_ID', 'Group', 'Ksat_first_row', 'hMinObs_cm'}), T];

% long table: one row per soil, model and procedure
F = table();
vn = {'alpha', 'n', 'theta_r', 'theta_s', 'SSE', 'RMSEth', 'RMSEK_orig', 'RMSEK_orig_l05', 'RMSEK_complex', 'exitflag', 'funcCount', ...
    'thrGtThs', 'nLt1', 'sign0', 'unbounded0', 'strictPos', 'negPlot', 'negMeas', 'neg1cm', 'negObs', 'nObsNeg', ...
    'hNegMax', 'Cmin', 'CminPlot', 'CmaxPlot', 'C1e6', 'thetaExcess', 'thetaExcessObs', 'dCorig'};
for M = S.models
    for pr = S.procs
        px = char(M + "_" + pr + "_");
        tt = table(T.Sample_ID, T.Group, repmat(M, nS, 1), repmat(pr, nS, 1), 'VariableNames', {'Sample_ID', 'Group', 'Model', 'Procedure'});
        for v = vn, tt.(v{1}) = T.([px v{1}]); end
        F = [F; tt]; %#ok<AGROW>
    end
end

% summary by model and procedure
rows = cell(0, 19);
for M = S.models
    for pr = S.procs
        q = F(F.Model == M & F.Procedure == pr, :); ng = q.strictPos == 0;
        rows(end + 1, :) = {M, pr, height(q), mean(q.RMSEth, 'omitnan'), sum(q.thrGtThs == 1), sum(q.nLt1 == 1), ...
            sum(q.unbounded0 == 1), sum(q.exitflag == 0), sum(ng), sum(q.negPlot == 1), sum(q.negMeas == 1), ...
            sum(q.negObs == 1), mx(q.hNegMax(ng)), mn(q.Cmin(ng)), mn(q.CminPlot(q.negPlot == 1)), mx(q.thetaExcess), ...
            sum(q.thrGtThs == 1 & q.strictPos == 1), mean(q.RMSEK_orig, 'omitnan'), mx(q.dCorig)}; %#ok<AGROW>
    end
end
SUM = cell2table(rows, 'VariableNames', {'Model', 'Procedure', 'Nfits', 'mean_RMSE_theta', 'N_theta_r_gt_theta_s', ...
    'N_n_lt_1', 'N_C_unbounded_at_h0', 'N_exitflag0', 'N_C_negative_somewhere', 'N_C_negative_plotted_range', ...
    'N_C_negative_h_ge_0p01cm', 'N_C_negative_measured_range', 'max_hNeg_max_cm', 'min_C_min', 'min_C_plotted', ...
    'max_theta_excess', 'N_theta_r_gt_theta_s_but_C_positive', 'mean_RMSE_K_original_code', ...
    'max_rel_diff_C_original_formula'});
for v = ["Model" "Procedure"], SUM.(v) = string(SUM.(v)); end

% Table 3 as computed by the original code (all soils of a group; RMSE_K with Ks = first K row)
rows = cell(0, 10);
for pr = S.procs
    for g = 1:numel(S.groups)
        for M = S.models
            q = F(F.Model == M & F.Procedure == pr & F.Group == S.groups(g), :);
            if isempty(q), continue, end
            ms = S.ms3.(char(M));
            rows(end + 1, :) = {pr, S.groups(g), M, height(q), mean(q.RMSEth), mean(q.RMSEK_orig), mean(q.RMSEK_orig_l05), ...
                ms(1, g), ms(2, g), sum(q.RMSEK_complex == 1)}; %#ok<AGROW>
        end
    end
end
T3 = cell2table(rows, 'VariableNames', {'Procedure', 'Group', 'Model', 'Nsoils', 'meanRMSE_theta', ...
    'meanRMSE_K_original_code', 'meanRMSE_K_l05', 'published_RMSE_theta', 'published_RMSE_K', 'N_complex_RMSE_K'});
for v = ["Procedure" "Group" "Model"], T3.(v) = string(T3.(v)); end
same = @(v, p, d) abs(round(v, d) - p) < 1e-9 | abs(floor(v*10^d + 1e-9)/10^d - p) < 1e-9;   % rounded or truncated
T3.match_theta = same(T3.meanRMSE_theta, T3.published_RMSE_theta, 4);
T3.match_K = same(T3.meanRMSE_K_original_code, T3.published_RMSE_K, 2);
T3.match_K_l05 = same(T3.meanRMSE_K_l05, T3.published_RMSE_K, 2);

% Table 2 parameters of the German Table 1 soils against the published values (3 decimals; a value
% matches when it rounds or truncates to the published one; values below 0.001 within a factor 2)
rows = cell(0, 6);
for pr = S.procs
    for M = S.models
        ms = S.ms2.(char(M)); ok = nan(numel(S.T1ids), 4); dev = nan(numel(S.T1ids), 4);
        for k = 1:numel(S.T1ids)
            q = F(F.Model == M & F.Procedure == pr & F.Sample_ID == S.T1ids(k), :);
            if isempty(q), continue, end
            v = [q.theta_r q.theta_s q.alpha q.n]; dev(k, :) = abs(v - ms(k, :));
            ok(k, :) = abs(round(v, 3) - ms(k, :)) < 1e-9 | abs(floor(v*1e3 + 1e-9)/1e3 - ms(k, :)) < 1e-9 | ...
                (ms(k, :) < 1e-3 & v < 2*ms(k, :) & v > ms(k, :)/2);
        end
        rows(end + 1, :) = {pr, M, sum(isfinite(ok(:, 1))), sum(ok(:) == 1), sum(all(ok == 1, 2)), max(dev(:, 4))}; %#ok<AGROW>
    end
end
T2 = cell2table(rows, 'VariableNames', {'Procedure', 'Model', 'Nsoils', 'N_parameters_match', ...
    'N_soils_all_four_match', 'max_abs_dev_n'});
for v = ["Procedure" "Model"], T2.(v) = string(T2.(v)); end

% the soils whose previous EVM fit has theta_r > theta_s
two = F(ismember(F.Sample_ID, S.two), :);
act = unique(F.Sample_ID(F.Model == "EVM" & F.Procedure == "OB1" & F.thrGtThs == 1), 'stable');

% match of the previous procedure with the published Tables 2 (56 parameters) and 3 (20 entries)
score = zeros(2, numel(S.V));
for k = 1:numel(S.V)
    q = T3(T3.Procedure == S.V(k).tag, :);
    score(2, k) = sum(q.match_theta) + sum(q.match_K & q.Model == "PDI") + sum(q.match_K_l05 & q.Model == "EVM");
    score(1, k) = sum(T2.N_parameters_match(T2.Procedure == S.V(k).tag));
end
best = string(S.V(1).tag);

% Table 6: negative water capacity of the previous (published) fits and of the constrained fits
rows = cell(0, 23);
for M = S.models
    for pr = [best "P2"]
        q = F(F.Model == M & F.Procedure == pr, :); ng = q.strictPos == 0;
        rows(end + 1, :) = {M, pr, height(q), sum(q.thrGtThs == 1), sum(q.nLt1 == 1), sum(ng), sum(q.negPlot == 1), ...
            sum(q.negMeas == 1), sum(q.neg1cm == 1), sum(q.negObs == 1), sum(q.nObsNeg(ng)), md(q.hNegMax(ng)), ...
            mx(q.hNegMax(ng)), md(q.thetaExcess(ng)), mx(q.thetaExcess(ng)), md(q.thetaExcessObs(ng)), ...
            mx(q.thetaExcessObs(ng)), md(q.Cmin(ng)), mn(q.Cmin(ng)), md(-q.Cmin(ng)./q.CmaxPlot(ng)), ...
            mean(q.RMSEth, 'omitnan'), mean(q.RMSEth(ng)), sum(q.thrGtThs == 1 & q.strictPos == 1)}; %#ok<AGROW>
    end
end
T6 = cell2table(rows, 'VariableNames', {'Model', 'Procedure', 'Nfits', 'N_theta_r_gt_theta_s', 'N_n_lt_1', 'N_C_negative', ...
    'N_C_negative_plotted_range', 'N_C_negative_h_ge_0p01cm', 'N_C_negative_h_ge_1cm', 'N_C_negative_in_measured_range', ...
    'N_retention_obs_with_C_negative', 'median_hNeg_max_cm', 'max_hNeg_max_cm', 'median_theta_excess', ...
    'max_theta_excess', 'median_theta_excess_measured_range', 'max_theta_excess_measured_range', 'median_C_min', ...
    'min_C_min', 'median_neg_to_max_C', 'mean_RMSE_theta', 'mean_RMSE_theta_C_negative_fits', ...
    'N_theta_r_gt_theta_s_but_C_positive'});
for v = ["Model" "Procedure"], T6.(v) = string(T6.(v)); end

% consistency with Script S6 (constrained fits) when its output exists
chk.nS6 = 0; chk.dS6 = NaN;
f6 = fullfile(dataDir, 'manuscript_results', 'EVM_parameter_sets.csv');
if exist(f6, 'file')
    P6 = readtable(f6, 'TextType', 'string'); d = [];
    for M = S.models
        for pr = "P2"
            a = F(F.Model == M & F.Procedure == pr, :); b = P6(P6.Model == M & P6.Procedure == pr, :);
            [tf, j] = ismember(a.Sample_ID, b.Sample_ID);
            d = [d; abs(a.RMSEth(tf) - b.RMSE_theta(j(tf)))]; %#ok<AGROW>
            chk.nS6 = chk.nS6 + sum(tf);
        end
    end
    chk.dS6 = max(d);
end

%% ---------------- outputs ----------------
lab = @(t) replace(replace(t, best, "previous"), "P2", "constrained");
Fw = F;     Fw.Procedure = lab(Fw.Procedure);
SUMw = SUM; SUMw.Procedure = lab(SUMw.Procedure);
T3w = T3(T3.Procedure == best, :); T3w.Procedure = lab(T3w.Procedure);
T6w = T6;   T6w.Procedure = lab(T6w.Procedure);
writetable(Fw,   fullfile(outDir, 'EVM_recheck_fits.csv'));
writetable(SUMw, fullfile(outDir, 'EVM_recheck_summary.csv'));
writetable(T3w,  fullfile(outDir, 'EVM_recheck_table3.csv'));
writetable(T6w,  fullfile(outDir, 'EVM_recheck_table6.csv'));
fid = fopen(fullfile(outDir, 'EVM_recheck_report.txt'), 'w');
if fid < 0, error('Cannot write the report in %s', outDir); end
for f = [1, fid]
    print_report(f, S, F, SUM, T3, T2, T6, two, act, score, best, chk);
end
fclose(fid);
R = struct('settings', S, 'bySoil', T, 'fits', F, 'summary', SUM, 'table3', T3, 'table2', T2, 'table6', T6, ...
    'twoSoils', two, 'bestVariant', best, 'score', score, 'checks', chk, 'data', D, 'ids', ids);
save(fullfile(outDir, 'EVM_recheck_results.mat'), 'R');
make_figures(outDir, S, F, D, ids, best);
fprintf('Outputs written to %s (%.1f min).\n', outDir, toc(tStart)/60);
fprintf('EVM_RECHECK_DONE\n');
end

%% =====================================================================
%% local functions
%% =====================================================================
function v = mx(x), x = x(isfinite(x)); if isempty(x), v = NaN; else, v = max(x); end, end
function v = md(x), x = x(isfinite(x)); if isempty(x), v = NaN; else, v = median(x); end, end
function v = mn(x), x = x(isfinite(x)); if isempty(x), v = NaN; else, v = min(x); end, end

function [fE, fP] = orig_models()
% the two retention models of SWRFIT_PDI_EVM.m, copied verbatim (x = [alpha theta_r n theta_s x5])
fE = @(vGx,h) (vGx(4)-vGx(2))*((1./(1+(vGx(1)*h).^(vGx(3)))).^(1-1/vGx(3))-(1./(1+(vGx(1)*10^7).^(vGx(3)))).^(1-1/vGx(3)))/(1-(1./(1+(vGx(1)*10^7).^(vGx(3)))).^(1-1/vGx(3)))...
    +vGx(2)*(1-log(1+(vGx(1)*h).^(1+0.1/vGx(3)))/log(1+(vGx(1)*10^7).^(1+0.1/vGx(3))));
fP = @(vGx,h) (vGx(4)-vGx(2))*((1./(1+(vGx(1)*h).^vGx(3))).^(1-1/vGx(3))-(1./(1+(vGx(1)*10^7).^vGx(3))).^(1-1/vGx(3)))/(1-(1./(1+(vGx(1)*10^7).^vGx(3))).^(1-1/vGx(3)))...
    +vGx(2)*(log(10^7./h)-0.1*log(10)*(1+2*(1-exp(-(vGx(2)/(vGx(4)-vGx(2)))^2))/vGx(3)^2)*log(1+(1./(vGx(1)*h)).^(1/(0.1*log(10)*(1+2*(1-exp(-(vGx(2)/(vGx(4)-vGx(2)))^2))/vGx(3)^2)))))/log(vGx(1)*10^7);
end

function r = analyse_soil(d, S, opts, optsLC)
h = 10.^d.pF; th = d.th;
[fE, fP] = orig_models();
op = optimoptions('lsqcurvefit', 'Display', 'off');      % all other options at their defaults
ws = warning('off', 'all');
r = struct();
% original code: one lsqcurvefit run
for k = 1:numel(S.V)
    V = S.V(k);
    for M = S.models
        px = char(M + "_" + V.tag + "_");
        if M == "PDI" && ~V.fitPDI                           % same x0 and bounds as the previous variant
            src = char("PDI_" + S.V(k - 1).tag + "_");
            for fn = fieldnames(r)', if startsWith(fn{1}, src), r.([px fn{1}(numel(src) + 1:end)]) = r.(fn{1}); end, end
            continue
        end
        x = nan(5, 1); rn = NaN; ef = NaN; fc = NaN;
        try
            if M == "EVM"
                [x, rn, ~, ef, o] = lsqcurvefit(fE, V.x0, h, th, V.lbE, [1; 1; 10; 1; 1], op);
            else
                [x, rn, ~, ef, o] = lsqcurvefit(fP, V.x0, h, th, [0; 0; 1; 0], [1; 1; 10; 1], op);
            end
            fc = o.funcCount;
        catch
        end
        r = store(r, px, [x(1) x(3) x(2) x(4)], rn, ef, fc, M, d, S, fE, fP);
    end
end
% theta_s > theta_r >= 0 (P2 of Scripts S2 and S6; the independent-bounds search P1 is internal)
for M = S.models
    Fp = fit_procedures(h, th, char(M), S, opts, optsLC);
    for pr = "P2"
        p = Fp.(pr).p;
        rn = sum((shp_theta(h, p, char(M), S) - th).^2);
        r = store(r, char(M + "_" + pr + "_"), p, rn, NaN, NaN, M, d, S, fE, fP);
    end
end
warning(ws);
end

function r = store(r, px, p, rn, ef, fc, M, d, S, fE, fP)
% parameters, errors and capacity diagnostics of one fit (p = [alpha n theta_r theta_s])
h = 10.^d.pF;
x = [p(1); p(3); p(2); p(4); 0.01];
r.([px 'alpha']) = p(1); r.([px 'n']) = p(2); r.([px 'theta_r']) = p(3); r.([px 'theta_s']) = p(4);
r.([px 'SSE']) = rn; r.([px 'exitflag']) = ef; r.([px 'funcCount']) = fc;
if any(~isfinite(p))                                         % the fit failed
    r.([px 'RMSEth']) = NaN; r.([px 'RMSEK_orig']) = NaN; r.([px 'RMSEK_complex']) = NaN; r.([px 'RMSEK_orig_l05']) = NaN;
    c = cap_diag([], char(M), S, h);
else
    if M == "EVM", f = fE; else, f = fP; end
    r.([px 'RMSEth']) = sqrt(mean((f(x, h) - d.th).^2));
    [r.([px 'RMSEK_orig']), r.([px 'RMSEK_complex'])] = orig_rmseK(M, x, d.hm, d.Km, fE, fP, false);
    r.([px 'RMSEK_orig_l05']) = orig_rmseK(M, x, d.hm, d.Km, fE, fP, true);
    c = cap_diag(p, char(M), S, h);
    Co = orig_capacity(M, x, S.hp); Cm = shp_capacity(S.hp, p, char(M), S);
    c.dCorig = max(abs(Co - Cm))/max(abs(Cm));
end
for fn = fieldnames(c)', r.([px fn{1}]) = c.(fn{1}); end
end

function c = cap_diag(p, mdl, S, hObs)
% sign, size and suction range of negative water capacity (flags as 0/1; NaN when the fit failed)
fn = {'thrGtThs', 'nLt1', 'sign0', 'unbounded0', 'strictPos', 'negPlot', 'negMeas', 'neg1cm', 'negObs', 'nObsNeg', ...
    'hNegMax', 'Cmin', 'CminPlot', 'CmaxPlot', 'C1e6', 'thetaExcess', 'thetaExcessObs', 'dCorig'};
c = cell2struct(num2cell(nan(numel(fn), 1)), fn, 1);
if isempty(p), return, end
C = shp_capacity(S.hg, p, mdl, S); Cp = shp_capacity(S.hp, p, mdl, S); th = shp_theta(S.hg, p, mdl, S);
neg = C < 0;
c.thrGtThs = double(p(3) > p(4)); c.nLt1 = double(p(2) < 1);
c.sign0 = lead_sign(p, mdl, S);
c.unbounded0 = double(p(2) < 1 && c.sign0 > 0);               % C -> +Inf as h -> 0+
c.strictPos = double(~any(neg) && ~any(Cp < 0) && c.sign0 > 0);
c.negPlot = double(any(Cp < 0)); c.negMeas = double(any(neg & S.hg >= S.hmin)); c.neg1cm = double(any(neg & S.hg >= 1));
c.negObs = double(any(neg & S.hg >= min(hObs)));               % theta rises inside the measured suction range
c.nObsNeg = sum(shp_capacity(hObs, p, mdl, S) < 0);           % retention observations on a rising branch
if any(neg), c.hNegMax = max(S.hg(neg)); c.Cmin = min(C(neg)); end
if c.sign0 < 0 && ~any(neg), c.hNegMax = S.hg(1); end         % only below 1e-12 cm
c.CminPlot = min(Cp); c.CmaxPlot = max(Cp); c.C1e6 = Cp(1);
c.thetaExcess = max(max(th) - p(4), 0); if c.thetaExcess < 1e-12, c.thetaExcess = 0; end
m = S.hg >= min(hObs);
c.thetaExcessObs = max(max(th(m)) - p(4), 0); if c.thetaExcessObs < 1e-12, c.thetaExcessObs = 0; end
end

function C = orig_capacity(M, x, h)
% -C_pdi and -C_evm of SWRFIT_PDI_EVM.m (lines 350-369 and 392-409), copied verbatim
alpha = x(1); qr = x(2); n = x(3); m = 1 - 1/n; qs = x(4); h0 = 10^7; ha = 1/alpha;
G = @(h) (1./(1+(alpha*h).^n)).^m;
if M == "PDI"
    b0 = 0.1*log(10); b1 = (qr/(qs-qr))^2; b = b0*(1+2*(1-exp(-b1))/(n^2));
    C_pdi = @(h) 1/(1-1*G(h0))*(qs-qr)*(-(1./(1 + (alpha*h).^n)).^m*m.*(alpha*h).^n*n./((1 + (alpha*h).^n).*h))...
        +1*qr*(-1./(h.*log(h0/ha).*(1+(ha./h).^(1/b))));
    C = -C_pdi(h);
else
    b1 = 1+0.1/n;
    C_evm = @(h) 1/(1-G(h0))*(qs-qr)*(-(1./(1 + (alpha*h).^n)).^m*m.*(alpha*h).^n*n./((1 + (alpha*h).^n).*h))...
        +qr*(-(h/ha).^b1*b1./(h.*(1 + (h/ha).^b1).*log(1 + (h0/ha)^b1)));
    C = -C_evm(h);
end
end

function [rk, cplx] = orig_rmseK(M, x, hm, Km, fE, fP, l05)
% RMSE_K of SWRFIT_PDI_EVM.m (lines 371-389 and 418-467): ks = Km(1), every conductivity row compared.
% l05 = true (EVM only): Knc of lines 402-404 and 448, i.e. Snc with beta and the factor Snc^0.5 (l = 0.5
% in Eq. 20), the form that reproduces the RMSE_K of EVM in the published Table 3
rk = NaN; cplx = NaN;
if isempty(Km), return, end
ks = Km(1);
alpha = x(1); qr = x(2); n = x(3); m = 1 - 1/n; qs = x(4); h0 = 10^7; ha = 1/alpha;
T = 293.15; rho_w = 1000; rho_sv = 0.001*exp(31.3716-6014.79/T-7.92495*0.001*T)/T;
Mw = 0.018015; g = 9.81; R = 8.314; Da = 2.414e-5*(T/273.15)^2;
G = @(h) (1./(1+(alpha*h).^n)).^m;
Sc = @(h) (G(h)-G(h0))/(1-G(h0));
if M == "PDI"
    b = 0.1*log(10)*(1+2*(1-exp(-(qr/(qs-qr))^2))/(n^2));
    Sm = fP(x, 10^5);
    Snc = @(h) (log(h0./h)-b*log(1+(ha./h).^(1/b)))/log(h0/ha);
    Kv = @(h) (100*24*60*60) * (rho_sv/rho_w)*Mw*g/(R*T).*((qs-fP(x,h)).^(10/3)./(qs^2)*Da).*exp(-h*10^(-2)*Mw*g./(R*T));
    K = @(h) ks/1*Sc(h).^(1/2).*(1-((1-G(h).^(1/m))/(1-G(h0).^(1/m))).^m).^2 + ...
        86.2*Sm*ha^(-1.5)/1*(h0/ha).^(-1.5*(1-Snc(h))) + 1*Kv(h);
else
    Sm = fE(x, 10^5);
    Snc = @(h) 1-log(G(h))/log(G(h0));                          % line 458 (exponent n, not beta)
    a = 0.75;
    Kv = @(h) (100*24*60*60) * (rho_sv/rho_w)*Mw*g/(R*T).*((qs-fE(x,h)).^(10/3)./(qs^2)*Da).*exp(-h*10^(-2)*Mw*g./(R*T));
    K = @(h) ks*Sc(h).^(1/2).*(1-((1-G(h).^(1/m))/(1-G(h0).^(1/m))).^m).^2 + ...
        86.2*Sm*alpha^(1.5)*Snc(h).^(0.0).*((alpha*h0).^(a*Snc(h)-a)).^2 + 1*Kv(h);
    if l05
        b1 = 1+0.1/n;
        Snc = @(h) (1-log(1+(alpha*h).^b1)/log(1+(alpha*h0).^b1));   % lines 402-404
        K = @(h) ks*Sc(h).^(1/2).*(1-((1-G(h).^(1/m))/(1-G(h0).^(1/m))).^m).^2 + ...
            86.2*Sm*alpha^(1.5)*Snc(h).^(0.5).*((alpha*h0).^(a*Snc(h)-a)).^2 + 1*Kv(h);
    end
end
v = sqrt(mean((log10(K(hm)) - log10(Km)).^2));
cplx = double(~isreal(v) && abs(imag(v)) > 0);
rk = real(v);
end

function s0 = lead_sign(p, mdl, S)
% sign of C as h -> 0+: C ~ A (alpha h)^ec + B (alpha h)^enc, the smaller exponent dominates
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n;
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
if s2 < s1*(1 - 1e-9), p1 = p2; s1 = s2; end
active = p1(3) > S.fMax*p1(4);
if ~active && s1 < s2*(1 - 1e-9), p2 = p1; s2 = s1; end
if active
    [px, sx] = fit_shp(h, th, mdl, 'P2x', start_set('P2x', th, false), S, opts, optsLC);
    if sx < s2*(1 - 1e-9) && px(3) < px(4), p2 = px; end
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
% PDI (Iden and Durner, 2014) and EVM (Eq. 7); valid also for n < 1 (m < 0)
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n;
G = (1 + (al*h).^n).^(-m);
G0 = (1 + (al*S.h0)^n)^(-m);
Sc = (G - G0)/(1 - G0);
if strcmp(mdl, 'PDI')
    if ths == thr, r2 = Inf; else, r2 = (thr/(ths - thr))^2; end
    b = 0.1 + 0.2/n^2*(1 - exp(-r2));
    u = log10(al*h);
    Snc = 1 - (max(u, 0) + b*log1p(exp(-abs(u)/b)))/log10(al*S.h0);
else
    beta = 1 + S.c/n;
    Snc = 1 - log1p((al*h).^beta)/log1p((al*S.h0)^beta);
end
th = (ths - thr)*Sc + thr*Snc;
end

function C = shp_capacity(h, p, mdl, S)
% C = -dtheta/dh (cm-1)
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n; x = al*h;
dG = m*n*al*x.^(n - 1)./(1 + x.^n).^(m + 1);
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

function print_report(f, S, F, SUM, T3, T2, T6, two, act, score, best, chk)
L = @(pr) char(replace(replace(string(pr), best, "previous"), "P2", "constrained"));
fprintf(f, '\n========== Previous fitting procedure vs theta_s > theta_r >= 0 (German database, %d soils) ==========\n', ...
    numel(unique(F.Sample_ID)));
fprintf(f, 'previous   : fitting code of the previous version (SWRFIT_PDI_EVM.m), one lsqcurvefit run per soil (default\n');
fprintf(f, '             options) from x0 = [alpha theta_r n theta_s] = [0.01 0.001 1.5 0.01]; 0 <= alpha, theta_r, theta_s <= 1,\n');
fprintf(f, '             1 <= n <= 10 (independent bounds; theta_s > theta_r not imposed)\n');
fprintf(f, 'constrained: theta_s > theta_r >= 0 enforced (theta_r = f theta_s, f <= %.4f), 30-110 starting points (Scripts S2, S6)\n', S.fMax);
fprintf(f, 'C < 0 checked on 1e-12 cm <= h < h0, on the plotted range 1e-6-1e7 cm, and as h -> 0+ (leading term)\n\n');
fprintf(f, '%-4s %-11s %5s %8s %6s %6s %6s %6s %6s %10s %10s %6s %7s\n', 'Mdl', 'Procedure', 'N', 'RMSEth', 'thr>s', 'C<0', ...
    'plot', 'h>.01', 'meas', 'hNegMax', 'Cmin', 'th>s+', 'RMSEK');
for i = 1:height(SUM)
    s = SUM(i, :);
    fprintf(f, '%-4s %-11s %5d %8.5f %6d %6d %6d %6d %6d %10.3g %10.3g %6d %7.3f\n', s.Model, L(s.Procedure), s.Nfits, ...
        s.mean_RMSE_theta, s.N_theta_r_gt_theta_s, s.N_C_negative_somewhere, s.N_C_negative_plotted_range, ...
        s.N_C_negative_h_ge_0p01cm, s.N_C_negative_measured_range, s.max_hNeg_max_cm, s.min_C_min, ...
        s.N_theta_r_gt_theta_s_but_C_positive, s.mean_RMSE_K_original_code);
end
fprintf(f, '   (thr>s = theta_r > theta_s; plot = C < 0 on 1e-6-1e7 cm; meas = C < 0 inside the measured range;\n');
fprintf(f, '    th>s+ = theta_r > theta_s but C > 0 everywhere; RMSEK = RMSE_K as computed by the previous code)\n');
fprintf(f, '   max relative difference between C of the previous formulas (C_pdi, C_evm) and Scripts S2-S6: %.2e\n', ...
    max(SUM.max_rel_diff_C_original_formula));
for M = S.models
    a = F(F.Model == M & F.Procedure == best, :); b = F(F.Model == M & F.Procedure == "P2", :);
    fprintf(f, '   %s mean RMSE_theta: previous %.6f, constrained %.6f\n', M, mean(a.RMSEth, 'omitnan'), mean(b.RMSEth, 'omitnan'));
end

fprintf(f, '\nPublished values reproduced by the previous procedure: %d of 56 Table 2 parameters (German soils), %d of 20 Table 3 entries\n', ...
    score(1, 1), score(2, 1));
fprintf(f, 'Table 3 as computed by the previous code (all soils of a group; RMSE_K with Ks = first K row, all rows;\n');
fprintf(f, 'l=0.5: EVM Knc with Snc(beta)^0.5 (lines 402-404, 448); = : equal to the published value rounded or truncated):\n');
fprintf(f, '%-5s %-4s %4s %9s %9s %8s %9s %8s\n', 'Group', 'Mdl', 'N', 'RMSEth', 'publ.', 'RMSEK', 'l=0.5', 'publ.');
q3 = T3(T3.Procedure == best, :);
for i = 1:height(q3)
    t = q3(i, :);
    fprintf(f, '%-5s %-4s %4d %9.5f%s %8.4f %8.3f%s %8.3f%s %8.2f\n', t.Group, t.Model, t.Nsoils, t.meanRMSE_theta, ...
        repmat('=', 1, t.match_theta), t.published_RMSE_theta, t.meanRMSE_K_original_code, repmat('=', 1, t.match_K), ...
        t.meanRMSE_K_l05, repmat('=', 1, t.match_K_l05), t.published_RMSE_K);
end
q2 = T2(T2.Procedure == best, :);
for i = 1:height(q2)
    fprintf(f, '   Table 2 %s: %2d of %2d parameters match, all four in %d of %d soils\n', q2.Model(i), q2.N_parameters_match(i), ...
        4*q2.Nsoils(i), q2.N_soils_all_four_match(i), q2.Nsoils(i));
end
for M = S.models
    nm = 0;
    for k = 1:numel(S.T1ids)
        q = F(F.Model == M & F.Procedure == best & F.Sample_ID == S.T1ids(k), :);
        if isempty(q), continue, end
        v = q.RMSEK_orig; if M == "EVM", v = q.RMSEK_orig_l05; end
        nm = nm + (abs(round(v, 2) - S.ms2K.(char(M))(k)) < 1e-9 || abs(floor(v*100 + 1e-9)/100 - S.ms2K.(char(M))(k)) < 1e-9);
    end
    fprintf(f, '   Table 2 RMSE_K %s (EVM with l = 0.5): %d of %d published values reproduced\n', M, nm, numel(S.T1ids));
end

q = F(F.Model == "EVM" & F.Procedure == best & F.thrGtThs == 1, :); nStar = (1 + sqrt(1 + 4*S.c))/2;
fprintf(f, '\nPrevious EVM fits with theta_r > theta_s (%d):\n', height(q));
for i = 1:height(q)
    fprintf(f, '   %-9s thr %.3f ths %.3f n %.4f | C > 0 everywhere: %d | C < 0 for h < %.2g cm (plotted range: %d)\n', ...
        q.Sample_ID(i), q.theta_r(i), q.theta_s(i), q.n(i), q.strictPos(i), q.hNegMax(i), q.negPlot(i));
end
fprintf(f, '   for theta_r > theta_s, C < 0 near saturation iff c/n > n - 1, i.e. n < %.4f: %d of these fits\n', nStar, sum(q.n < nStar));
q = F(F.Model == "PDI" & F.Procedure == best & F.thrGtThs == 1, :);
b = 0.1 + 0.2./q.n.^2.*(1 - exp(-(q.theta_r./(q.theta_s - q.theta_r)).^2)); e = 1./(b*log(10)) - 1; ec = q.n - 1;
fprintf(f, 'Previous PDI fits with theta_r > theta_s: %d (%s), with C < 0: %d\n', height(q), ...
    strjoin(arrayfun(@(g) sprintf('%s %d', g, sum(q.Group == g)), S.groups, 'UniformOutput', false), ', '), sum(q.strictPos == 0));
fprintf(f, '   film exponent 1/(b ln10) - 1 = %.2f-%.2f vs capillary exponent n - 1 = %.2f-%.2f; film exponent larger in %d\n', ...
    min(e), max(e), min(ec), max(ec), sum(e > ec));
for M = S.models
    ids = F.Sample_ID(F.Model == M & F.Procedure == best & F.thrGtThs == 1);
    a = F(F.Model == M & F.Procedure == best & ismember(F.Sample_ID, ids), :);
    c = F(F.Model == M & F.Procedure == "P2" & ismember(F.Sample_ID, ids), :);
    fprintf(f, '   %s: these %d soils, mean RMSE_theta previous %.5f -> constrained %.5f\n', M, numel(ids), mean(a.RMSEth), mean(c.RMSEth));
end
fprintf(f, '\nTUBS_205 (Figure 3):\n');
for M = S.models
    for pr = [best "P2"]
        q = F(F.Sample_ID == "TUBS_205" & F.Model == M & F.Procedure == pr, :);
        if isempty(q), continue, end
        fprintf(f, '   %s %-11s alpha %.4f n %.4f thr %.4f ths %.4f | RMSE %.5f | C < 0: %d (h < %.2g cm), theta - theta_s max %.4f, C_min %.2e, C_max %.2e\n', ...
            M, L(pr), q.alpha, q.n, q.theta_r, q.theta_s, q.RMSEth, q.strictPos == 0, q.hNegMax, q.thetaExcess, q.Cmin, q.CmaxPlot);
    end
end

fprintf(f, '\nTable 6:\n');
fprintf(f, '%-4s %-11s %6s %5s %5s %6s %6s %6s %5s %5s %13s %17s %17s %17s %8s %8s\n', 'Mdl', 'Procedure', 'thr>s', 'C<0', ...
    '>1e-6', '>.01cm', '>1cm', 'meas', 'nObs', 'th>s+', 'hNegMax md/mx', 'excess md/mx', 'excess obs md/mx', ...
    'Cmin md/mn', 'RMSEth', 'RMSEneg');
for i = 1:height(T6)
    t = T6(i, :);
    fprintf(f, '%-4s %-11s %6d %5d %5d %6d %6d %6d %5d %5d %6.3g/%-6.3g %8.2g/%-8.2g %8.2g/%-8.2g %8.2g/%-8.2g %8.5f %8.5f\n', ...
        t.Model, L(t.Procedure), t.N_theta_r_gt_theta_s, t.N_C_negative, t.N_C_negative_plotted_range, ...
        t.N_C_negative_h_ge_0p01cm, t.N_C_negative_h_ge_1cm, t.N_C_negative_in_measured_range, ...
        t.N_retention_obs_with_C_negative, t.N_theta_r_gt_theta_s_but_C_positive, t.median_hNeg_max_cm, t.max_hNeg_max_cm, ...
        t.median_theta_excess, t.max_theta_excess, t.median_theta_excess_measured_range, t.max_theta_excess_measured_range, ...
        t.median_C_min, t.min_C_min, t.mean_RMSE_theta, t.mean_RMSE_theta_C_negative_fits);
    fprintf(f, '     median |C_min|/C_max of the fits with C < 0: %.2f\n', t.median_neg_to_max_C);
end
fprintf(f, '\nConsistency with Script S6: %d constrained fits compared, max |dRMSE_theta| %.2e\n', chk.nS6, chk.dS6);
end

function make_figures(outDir, S, F, D, ids, best)
set(0, 'DefaultFigureVisible', 'off');
% Figures 2 and 3 of the revised manuscript: the previous (published) fits and theta_s > theta_r >= 0.
% Every figure is smaller than the screen: larger off-screen figures lose the top of their axes
% contents in exportgraphics.
blue = [0 0.2 0.9]; red = [0.85 0 0]; P2txt = '\theta_s > \theta_r \geq 0';
hh = logspace(-2, log10(S.h0) - 1e-9, 500)'; hc = logspace(-6, log10(S.h0) - 1e-9, 600)';
grp = extractBefore(ids, "_");
parts = {S.groups(1:3), S.groups(4:5)};
cols = {"PDI", best, 'PDI, independent bounds'; "EVM", best, 'EVM, independent bounds'; ...
        "PDI", "P2", ['PDI, ' P2txt]; "EVM", "P2", ['EVM, ' P2txt]};
let = 'abcdefghijklmnop';
for part = 1:2
    gl = parts{part};
    f2 = figure('Position', [10 10 1640 70 + 300*numel(gl)], 'Color', 'w');
    tl = tiledlayout(f2, numel(gl), 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for r = 1:numel(gl)
        sel = ids(grp == gl(r));
        CC = cell(1, 4); yl = [0 1e-9]; sp = cell(1, 4); spl = cell(1, 4);
        for c = 1:4
            q = F(F.Model == cols{c, 1} & F.Procedure == cols{c, 2}, :);
            [~, j] = ismember(sel, q.Sample_ID); q = q(j, :);
            Cm = zeros(numel(hc), height(q));
            for k = 1:height(q), Cm(:, k) = 1e3*shp_capacity(hc, [q.alpha(k) q.n(k) q.theta_r(k) q.theta_s(k)], char(cols{c, 1}), S); end
            CC{c} = Cm; sp{c} = q.strictPos; spl{c} = q.negPlot;
            if ~isempty(Cm), yl = [min(yl(1), min(Cm(:))) max(yl(2), max(Cm(:)))]; end
        end
        hiNeg = sp{1} == 0;                               % soils whose previous PDI fit has C < 0
        for c = 1:4
            ax = nexttile(tl); hold(ax, 'on');
            if ~isempty(CC{c}), plot(ax, hc, CC{c}, 'LineWidth', 0.5); end
            if any(hiNeg), plot(ax, hc, CC{c}(:, hiNeg), 'k-', 'LineWidth', 1.4); end
            yline(ax, 0, 'k-');
            set(ax, 'XScale', 'log', 'XTick', 10.^(-6:2:6), 'FontSize', 8); xlim(ax, [1e-6 S.h0]);
            ylim(ax, yl + [-0.04 0.06]*diff(yl)); grid(ax, 'on'); box(ax, 'on');
            nneg = sum(sp{c} == 0); npl = sum(spl{c} == 1);
            if nneg == npl, txt = sprintf('n = %d, C < 0 in %d', numel(sel), nneg);
            else,           txt = sprintf('n = %d, C < 0 in %d (only at h < 10^{-9} cm)', numel(sel), nneg);
            end
            title(ax, sprintf('(%s) %s: %s', let((r - 1)*4 + c), gl(r), cols{c, 3}), 'FontSize', 9);
            text(ax, 0.97, 0.06, txt, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 8, 'BackgroundColor', 'w');
            if c == 1, ylabel(ax, 'C (10^{-3} cm^{-1})'); end
            if r == numel(gl), xlabel(ax, 'Suction, h (cm)'); end
        end
    end
    exportgraphics(f2, fullfile(outDir, sprintf('Fig2_previous_vs_constrained_part%d.png', part)), 'Resolution', 300); close(f2);
end
i = find(ids == "TUBS_205", 1);
if ~isempty(i)
    d = D(i); f3 = figure('Position', [10 10 1500 600], 'Color', 'w');
    sets = {"PDI", best, red, ':', 1.8; "PDI", "P2", red, '-', 3.2; "EVM", best, blue, ':', 1.8; "EVM", "P2", blue, '-', 1.3};
    pp = cell(1, 4); lab = cell(1, 4);
    for s = 1:4
        q = F(F.Sample_ID == "TUBS_205" & F.Model == sets{s, 1} & F.Procedure == sets{s, 2}, :);
        pp{s} = [q.alpha q.n q.theta_r q.theta_s];
        if sets{s, 2} == "P2", pr = P2txt; else, pr = 'independent bounds'; end
        lab{s} = sprintf('%s, %s (\\theta_r = %.3f, \\theta_s = %.3f)', sets{s, 1}, pr, pp{s}(3), pp{s}(4));
    end
    a1 = subplot(1, 2, 1); hold(a1, 'on');
    plot(a1, 10.^d.pF, d.th, 'ko', 'MarkerSize', 3);
    for s = 1:4, plot(a1, hh, shp_theta(hh, pp{s}, char(sets{s, 1}), S), sets{s, 4}, 'Color', sets{s, 3}, 'LineWidth', sets{s, 5}); end
    set(a1, 'XScale', 'log', 'XTick', 10.^(-2:7)); xlim(a1, [1e-2 S.h0]); ylim(a1, [0 0.6]); grid(a1, 'on'); box(a1, 'on');
    xlabel(a1, 'Suction, h (cm)'); ylabel(a1, 'Water content, \theta (cm^3 cm^{-3})');
    legend(a1, ['SWC data', lab], 'Location', 'southwest', 'FontSize', 8);
    title(a1, '(a) Silt Loam: TUBS-205, soil water characteristic (fitted)');
    ps = a1.Position;                                     % inset in the empty wet, mid-theta corner
    ai = axes(f3, 'Position', [ps(1) + 0.06*ps(3), ps(2) + 0.30*ps(4), 0.36*ps(3), 0.37*ps(4)]); hold(ai, 'on');
    hz = logspace(-6, 2, 400)'; plot(ai, 10.^d.pF, d.th, 'ko', 'MarkerSize', 3);
    yz = d.th(d.pF <= 2);
    for s = 1:4
        tz = shp_theta(hz, pp{s}, char(sets{s, 1}), S); yz = [yz; tz]; %#ok<AGROW>
        plot(ai, hz, tz, sets{s, 4}, 'Color', sets{s, 3}, 'LineWidth', sets{s, 5});
    end
    set(ai, 'XScale', 'log', 'XTick', 10.^(-6:2:2), 'FontSize', 7); xlim(ai, [1e-6 1e2]);
    ylim(ai, [min(yz) max(yz)] + [-0.1 0.1]*(max(yz) - min(yz))); grid(ai, 'on'); box(ai, 'on');
    title(ai, 'near saturation', 'FontSize', 8, 'FontWeight', 'normal');
    a2 = subplot(1, 2, 2); hold(a2, 'on'); yc = [0 0];
    for s = 1:4
        Cs = 1e3*shp_capacity(hc, pp{s}, char(sets{s, 1}), S);
        plot(a2, hc, Cs, sets{s, 4}, 'Color', sets{s, 3}, 'LineWidth', sets{s, 5}); yc = [min(yc(1), min(Cs)) max(yc(2), max(Cs))];
    end
    yline(a2, 0, 'k-'); set(a2, 'XScale', 'log', 'XTick', 10.^(-6:2:6)); xlim(a2, [1e-6 S.h0]); grid(a2, 'on'); box(a2, 'on');
    ylim(a2, yc + [-0.6 0.08]*diff(yc));                  % room for the legend below the curves
    xlabel(a2, 'Suction, h (cm)'); ylabel(a2, 'Soil water capacity, C (10^{-3} cm^{-1})');
    legend(a2, lab, 'Location', 'southeast', 'FontSize', 8); title(a2, '(b) Silt Loam: TUBS-205, soil water capacity');
    exportgraphics(f3, fullfile(outDir, 'Fig3_TUBS205_previous_vs_constrained.png'), 'Resolution', 300); close(f3);
end
end
