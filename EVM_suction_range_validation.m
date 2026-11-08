function R = EVM_suction_range_validation(dataDir, nSoilsMax)
%EVM_SUCTION_RANGE_VALIDATION  Validation of the van Genuchten (VG), PDI and EVM
% models by suction range, using the German database of Hohenbrink et al. (2023):
% RetMeas.csv and CondMeas.csv (572 soils; 409 with measured KSAT).
%
%   R = EVM_suction_range_validation              % CSV files next to this file
%   R = EVM_suction_range_validation(dataDir)     % CSV files in dataDir
%   R = EVM_suction_range_validation(dataDir, 20) % quick test on the first 20 soils
%
% Answers Reviewer #2 ("Report validation by suction range"):
%   A. Suction ranges (h = 10^pF cm): wet pF < 2, intermediate 2 <= pF < 4.2 (to the
%      permanent wilting point), dry pF >= 4.2 (adsorption-dominated range).
%   B. Retention fits, identical for the three models (Eq. 22, lsqcurvefit, unweighted
%      theta residuals): log10(alpha) in [-5, 0], n in [1.0001, 10] and theta_s > theta_r >= 0
%      enforced as theta_r = f*theta_s (0 <= f <= 0.9999, 0.01 <= theta_s <= 1). 30 starting
%      points (110 where the constraint is active), aligned with the independent-bounds fit and
%      checked with a linear-inequality fit, as in Script S2.
%      Sensitivity W: the same constrained fit with range-balanced weights, so that every
%      suction range with data contributes equally to the objective.
%   C. Conductivity predicted from the retention parameters, with Ks,c = measured KSAT
%      (the KSAT value itself is not scored; errors in log10 units, Eq. 24):
%      VG  : Ks Se^0.5 [1 - (1 - Se^(1/m))^m]^2 (van Genuchten, 1980), capillary only
%      PDI : Kc (Eq. 14) + Knc (Eqs. 16-17) + Kv (Appendix A); Snc of Iden and Durner (2014)
%      EVM : Kc (Eq. 14) + Knc (Eq. 20; l = 0, gamma = 0.75, eta = 2) + Kv (Appendix A)
%   D. Errors per range (RMSE and mean error, per soil and pooled), the coverage of the
%      conductivity data, and the share of the non-capillary terms, (Knc + Kv)/K, in the
%      predicted K at every observation (> 0.5: non-capillary dominated).
%   E. Overlap with the 12-soil development set (Table 1): the 7 German Table 1 soils are
%      flagged (Dev12) and the summary is repeated without them (NoDev12).
%
% Suction h is a positive magnitude (cm), h0 = 1e7 cm, ha = 1/alpha; EVM uses
% beta = 1 + 0.1/n. PDI uses the same h0, which reproduces the PDI parameters of Table 2.
%
% Outputs (folder range_results next to the data):
%   EVM_range_report.txt      all numbers quoted in the manuscript and the reply
%   EVM_range_by_soil.csv     per soil and model: parameters and errors by range (Text S4)
%   EVM_range_summary.csv     errors by model, procedure, subset, variable and range (Table 4)
%   EVM_range_coverage.csv    retention and conductivity observations by group and range
%   FigS8_range_coverage.png, FigS9_range_errors.png, FigS10_range_examples.png,
%   EVM_range_results.mat
% Runtime: about 10-20 min with 24 parallel workers (R2025b, Optimization and Parallel
% Computing Toolboxes; parfor runs serially without a pool).

if nargin < 1 || isempty(dataDir)
    dataDir = fileparts(mfilename('fullpath'));
end
outDir = fullfile(dataDir, 'range_results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
tStart = tic;

%% ---------------- settings ----------------
S.h0     = 1e7;        % oven-dry suction (cm), Schneider and Goss (2012)
S.c      = 0.1;        % EVM closure beta = 1 + c/n
S.fMax   = 0.9999;     % theta_r <= fMax*theta_s, i.e. theta_s > theta_r
S.nMin   = 1.0001;     % lower bound of n (n > 1)
S.pFedge = [2 4.2];    % wet | intermediate | dry
S.rn     = ["wet" "int" "dry" "all"];
S.models = ["VG" "PDI" "EVM"];
S.dev12  = ["FSUJ_023" "FSUJ_035" "TUBS_027" "TUBS_128" "KIT_022" "KIT_009" "TUBS_106"];  % German soils of Table 1
S.hg     = logspace(-12, log10(S.h0) - 1e-9, 20000)';   % monotonicity check on 0 < h < h0
S.hx     = logspace(-2, log10(S.h0) - 1e-9, 4000)';     % capillary/non-capillary crossover of K
S.examples = ["KIT_022" "TUBS_168" "TUBS_119"];          % Figure S10
opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxFunctionEvaluations', 4000, 'MaxIterations', 1000, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-10);
optsLC = optimoptions('lsqcurvefit', 'Algorithm', 'interior-point', 'Display', 'off', ...
    'MaxFunctionEvaluations', 8000, 'MaxIterations', 2000, 'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-10, 'ConstraintTolerance', 1e-9);

%% ---------------- data ----------------
ret  = readtable(fullfile(dataDir, 'RetMeas.csv'),  'TextType', 'string');
cond = readtable(fullfile(dataDir, 'CondMeas.csv'), 'TextType', 'string');
ret  = ret(isfinite(ret.pF) & isfinite(ret.theta), :);
ids  = unique(ret.Sample_ID, 'stable');
if nargin > 1 && ~isempty(nSoilsMax)
    if isnumeric(nSoilsMax), ids = ids(1:min(nSoilsMax, numel(ids)));
    else,                    ids = ids(ismember(ids, string(nSoilsMax)));
    end
end
nS   = numel(ids);
[~, ir] = ismember(ret.Sample_ID, ids);
[~, ic] = ismember(cond.Sample_ID, ids);
D = struct('id', cellstr(ids), 'pF', [], 'th', [], 'mt', [], 'pFK', [], 'K', [], 'mtK', [], 'Ksat', NaN);
for i = 1:nS
    s = ir == i;
    D(i).pF = ret.pF(s); D(i).th = ret.theta(s); D(i).mt = ret.MeasType(s);
    isKs  = ic == i & cond.MeasType == "KSAT";
    unsat = ic == i & cond.MeasType ~= "KSAT" & cond.k > 0;
    D(i).pFK = cond.pF(unsat); D(i).K = cond.k(unsat); D(i).mtK = cond.MeasType(unsat);
    if any(isKs), D(i).Ksat = cond.k(find(isKs, 1)); end
end
grp = extractBefore(ids, "_");
dev = ismember(ids, S.dev12);
fprintf('Loaded %d soils, %d retention points, %d conductivity points, %d soils with measured KSAT.\n', ...
    nS, sum(ir > 0), sum(ic > 0 & cond.MeasType ~= "KSAT"), sum(isfinite([D.Ksat])));

if isempty(gcp('nocreate'))          % parfor runs serially when no pool can start
    try
        pc = parcluster('Processes');             % in-session change only, profile not saved
        pc.NumWorkers = max(pc.NumWorkers, feature('numcores'));
        parpool(pc, feature('numcores'));
    catch
        try, parpool('Processes'); catch, end %#ok<NOCOM>
    end
end

%% ---------------- B-D. fits and errors by range ----------------
out = cell(nS, 1);
parfor i = 1:nS
    out{i} = analyse_soil(D(i), S, opts, optsLC);
end
fprintf('Fits finished (%.1f min).\n', toc(tStart)/60);

T = struct2table(cell2mat(cellfun(@(o) o.row, out, 'UniformOutput', false)));
T = [table(ids, grp, dev, [D.Ksat]', 'VariableNames', {'Sample_ID', 'Group', 'Dev12', 'Ksat_cm_d'}), T];

% pooled observations (retention: all soils; conductivity: soils with measured KSAT)
sid  = cell2mat(cellfun(@(o, i) repmat(i, numel(o.pF), 1), out, num2cell((1:nS)'), 'UniformOutput', false));
sidK = cell2mat(cellfun(@(o, i) repmat(i, numel(o.pFK), 1), out, num2cell((1:nS)'), 'UniformOutput', false));
Q.pF   = cell2mat(cellfun(@(o) o.pF,   out, 'UniformOutput', false));
Q.cls  = cell2mat(cellfun(@(o) o.cls,  out, 'UniformOutput', false));
Q.pFK  = cell2mat(cellfun(@(o) o.pFK,  out, 'UniformOutput', false));
Q.clsK = cell2mat(cellfun(@(o) o.clsK, out, 'UniformOutput', false));
Q.sid = sid; Q.sidK = sidK;
for M = cellstr(S.models)
    for pr = {'P2', 'W'}
        m = M{1}; p = pr{1};
        Q.(m).(p).th = cell2mat(cellfun(@(o) o.res.(m).(p).th, out, 'UniformOutput', false));
        Q.(m).(p).K  = cell2mat(cellfun(@(o) o.res.(m).(p).K,  out, 'UniformOutput', false));
        Q.(m).(p).sh = cell2mat(cellfun(@(o) o.res.(m).(p).sh, out, 'UniformOutput', false));
    end
end

%% ---------------- summaries ----------------
subs = ["All" "TUDD" "UFZ" "FSUJ" "KIT" "TUBS" "NoDev12" "Dev12"];
rows = cell(0, 12);
for M = S.models
    for pr = ["P2" "W"]
        for sb = subs
            ss = find(soil_subset(sb, grp, dev));
            for v = ["theta" "logK"]
                for k = 1:4
                    q = Q.(char(M)).(char(pr));
                    if v == "theta"
                        obs = ismember(Q.sid, ss);
                        if k < 4, obs = obs & Q.cls == k; end
                        e = q.th(obs); nnc = NaN;
                        col = T.(char(M + "_" + pr + "_RMSEth_" + S.rn(k)))(ss);
                    else
                        obs = ismember(Q.sidK, ss);
                        if k < 4, obs = obs & Q.clsK == k; end
                        e = q.K(obs); nnc = sum(q.sh(obs) >= 0.5);
                        col = T.(char(M + "_" + pr + "_RMSEK_" + S.rn(k)))(ss);
                    end
                    col = col(isfinite(col));
                    rows(end + 1, :) = {M, pr, sb, v, S.rn(k), numel(e), numel(col), mean(col), ...
                        median(col), sqrt(mean(e.^2)), mean(e), nnc}; %#ok<AGROW>
                end
            end
        end
    end
end
G = cell2table(rows, 'VariableNames', {'Model', 'Procedure', 'Subset', 'Variable', 'Range', 'Nobs', 'Nsoils', ...
    'meanRMSE', 'medianRMSE', 'pooledRMSE', 'pooledME', 'nNonCapDominated'});
for v = ["Model" "Procedure" "Subset" "Variable" "Range"], G.(v) = string(G.(v)); end

% coverage of the observations by group and range
glc = ["TUDD" "UFZ" "FSUJ" "KIT" "TUBS" "All"];
crow = cell(numel(glc), 18);
for j = 1:numel(glc)
    sel = find(grp == glc(j) | glc(j) == "All");
    pr = vertcat(D(sel).pF); pk = vertcat(D(sel).pFK); mk = string(vertcat(D(sel).mtK));   % string: groups may be empty
    ks = sel(isfinite([D(sel).Ksat]));
    crow(j, :) = {glc(j), numel(sel), numel(pr), sum(pr < S.pFedge(1)), sum(pr >= S.pFedge(1) & pr < S.pFedge(2)), ...
        sum(pr >= S.pFedge(2)), sum(arrayfun(@(d) any(d.pF >= S.pFedge(2)), D(sel))), numel(pk), ...
        sum(pk < S.pFedge(1)), sum(pk >= S.pFedge(1) & pk < S.pFedge(2)), sum(pk >= S.pFedge(2)), sum(pk >= 3), ...
        sum(pk >= 3.5), max([pk; NaN]), sum(mk == "AirEntry"), numel(ks), numel(vertcat(D(ks).pFK)), ...
        sum(arrayfun(@(d) ~isempty(d.pFK), D(sel)))};
end
C = cell2table(crow, 'VariableNames', {'Group', 'Nsoils', 'NobsRet', 'NobsRet_wet', 'NobsRet_int', 'NobsRet_dry', ...
    'NsoilsRetDry', 'NobsK', 'NobsK_wet', 'NobsK_int', 'NobsK_dry', 'NobsK_pF3', 'NobsK_pF35', 'maxpF_K', ...
    'NobsK_AirEntry', 'NsoilsKSAT', 'NobsK_KSATsoils', 'NsoilsK'});
C.Group = string(C.Group);

% consistency with Script S2 (same EVM procedure)
chk = struct('n', 0, 'dSSE', NaN, 'dRMSEK', NaN);
fS2 = fullfile(dataDir, 'positivity_results', 'EVM_positivity_fits.csv');
if exist(fS2, 'file')
    F2 = readtable(fS2, 'TextType', 'string');
    [tf, j] = ismember(ids, F2.Sample_ID);
    if any(tf)
        s1 = T.EVM_P2_RMSEth_all(tf).^2; s2 = F2.P2_RMSE_theta(j(tf)).^2;
        chk.n = sum(tf); chk.dSSE = max(abs(s1 - s2)./s2);
        chk.dRMSEK = max(abs(T.EVM_P2_RMSEK_all(tf) - F2.P2_RMSE_logK(j(tf))), [], 'omitnan');
    end
end

%% ---------------- outputs ----------------
vn = T.Properties.VariableNames;                % constrained fits only (P1 is internal)
Tw = T(:, ~contains(vn, '_P1_') & ~endsWith(vn, '_active') & ~endsWith(vn, '_nFix'));
writetable(Tw, fullfile(outDir, 'EVM_range_by_soil.csv'));
writetable(G, fullfile(outDir, 'EVM_range_summary.csv'));
writetable(C, fullfile(outDir, 'EVM_range_coverage.csv'));
fid = fopen(fullfile(outDir, 'EVM_range_report.txt'), 'w');
if fid < 0, error('Cannot write the report in %s', outDir); end
for f = [1, fid]
    print_report(f, S, T, G, C, Q, nS, chk);
end
fclose(fid);
make_figures(outDir, S, T, Q, D, ids);
R = struct('settings', S, 'bySoil', Tw, 'summary', G, 'coverage', C, 'pooled', Q);
save(fullfile(outDir, 'EVM_range_results.mat'), 'R');
fprintf('Outputs written to %s (%.1f min).\n', outDir, toc(tStart)/60);
fprintf('EVM_RANGE_VALIDATION_DONE\n');
end

%% =====================================================================
%% local functions
%% =====================================================================
function o = analyse_soil(d, S, opts, optsLC)
% Fits VG, PDI and EVM to one soil (P1, P2, W) and computes the errors by range.
t0 = tic;
h = 10.^d.pF; cls = range_class(d.pF, S);
w = zeros(size(h));                    % W: every range with data gets the same total weight
for k = 1:3
    s = cls == k;
    if any(s), w(s) = 1/sum(s); end
end
w = w*numel(h)/sum(w);
hasK = isfinite(d.Ksat) && ~isempty(d.K);
clsK = range_class(d.pFK, S);
r.N_th_wet = sum(cls == 1); r.N_th_int = sum(cls == 2); r.N_th_dry = sum(cls == 3);
r.N_K_wet = sum(clsK == 1); r.N_K_int = sum(clsK == 2); r.N_K_dry = sum(clsK == 3);
r.maxpF_K = NaN;
if ~isempty(d.pFK), r.maxpF_K = max(d.pFK); end
o.pF = d.pF; o.cls = cls; o.pFK = zeros(0, 1); o.clsK = zeros(0, 1);
if hasK, o.pFK = d.pFK; o.clsK = clsK; hK = 10.^d.pFK; lKm = log10(d.K); end
for M = S.models
    mdl = char(M);
    F = fit_procedures(h, d.th, w, mdl, S, opts, optsLC);
    r.([mdl '_active']) = F.active; r.([mdl '_nFix']) = F.nFix;
    p = F.P1.p;                         % independent bounds: Table 2/3 checks only
    r.([mdl '_P1_alpha']) = p(1); r.([mdl '_P1_n']) = p(2); r.([mdl '_P1_theta_r']) = p(3); r.([mdl '_P1_theta_s']) = p(4);
    r.([mdl '_P1_RMSEth']) = sqrt(mean((shp_theta(h, p, mdl, S) - d.th).^2));
    r.([mdl '_P1_RMSEK']) = NaN;
    if hasK, r.([mdl '_P1_RMSEK']) = sqrt(mean((shp_logK(hK, p, mdl, d.Ksat, S) - lKm).^2)); end
    [r.([mdl '_P1_nIncrease']), r.([mdl '_P1_hIncreaseMax'])] = monotonic_check(p, mdl, S);
    for pr = ["P2" "W"]
        p = F.(pr).p; px = [mdl '_' char(pr) '_'];
        r.([px 'alpha']) = p(1); r.([px 'n']) = p(2); r.([px 'theta_r']) = p(3); r.([px 'theta_s']) = p(4);
        res = shp_theta(h, p, mdl, S) - d.th;
        [rm, me] = range_err(res, cls);
        for k = 1:4
            r.([px 'RMSEth_' char(S.rn(k))]) = rm(k); r.([px 'MEth_' char(S.rn(k))]) = me(k);
        end
        [r.([px 'nIncrease']), r.([px 'hIncreaseMax'])] = monotonic_check(p, mdl, S);
        o.res.(mdl).(char(pr)).th = res;
        rm = nan(1, 4); me = nan(1, 4); nnc = nan(1, 4); pfx = NaN; resK = zeros(0, 1); sh = zeros(0, 1);
        if hasK
            [lK, sh] = shp_logK(hK, p, mdl, d.Ksat, S);
            resK = lK - lKm;
            [rm, me] = range_err(resK, clsK);
            nnc = [arrayfun(@(k) sum(sh(clsK == k) >= 0.5), 1:3), sum(sh >= 0.5)];
            if ~strcmp(mdl, 'VG')
                [~, shx] = shp_logK(S.hx, p, mdl, d.Ksat, S);
                k = find(shx >= 0.5, 1);
                if ~isempty(k), pfx = log10(S.hx(k)); end
            end
        end
        for k = 1:4
            r.([px 'RMSEK_' char(S.rn(k))]) = rm(k); r.([px 'MEK_' char(S.rn(k))]) = me(k);
            r.([px 'nNonCapK_' char(S.rn(k))]) = nnc(k);
        end
        r.([px 'pFcross_K']) = pfx;      % suction above which (Knc + Kv) > Kc
        o.res.(mdl).(char(pr)).K = resK; o.res.(mdl).(char(pr)).sh = sh;
    end
end
r.tFit_s = toc(t0);
o.row = r;
end

function F = fit_procedures(h, th, w, mdl, S, opts, optsLC)
% P1: independent bounds; P2: theta_s > theta_r >= 0 (theta_r = f theta_s), aligned with
% P1 as in Script S2 (the P2 set lies inside the P1 set); W: P2 with range-balanced weights.
[p1, ~, s1] = fit_shp(h, th, [], mdl, 'P1', start_set('P1', th, false), S, opts, optsLC);
[p2, ~, s2] = fit_shp(h, th, [], mdl, 'P2', start_set('P2', th, false), S, opts, optsLC);
dense = p1(3) > S.fMax*p1(4);        % constraint active: denser multi-start for both
if dense
    [pb, ~, sb] = fit_shp(h, th, [], mdl, 'P1', start_set('P1', th, true), S, opts, optsLC);
    if sb < s1, p1 = pb; s1 = sb; end
    [pb, ~, sb] = fit_shp(h, th, [], mdl, 'P2', start_set('P2', th, true), S, opts, optsLC);
    if sb < s2, p2 = pb; s2 = sb; end
end
F.nFix = 0;
if s2 < s1*(1 - 1e-9), p1 = p2; s1 = s2; F.nFix = F.nFix + 1; end
F.active = p1(3) > S.fMax*p1(4);
if ~F.active && s1 < s2*(1 - 1e-9), p2 = p1; s2 = s1; F.nFix = F.nFix + 1; end
if F.active
    [px, ~, sx] = fit_shp(h, th, [], mdl, 'P2x', start_set('P2x', th, false), S, opts, optsLC);
    if sx < s2*(1 - 1e-9) && px(3) < px(4), p2 = px; s2 = sx; F.nFix = F.nFix + 1; end
end
F.P1.p = p1; F.P1.sse = s1; F.P2.p = p2; F.P2.sse = s2;
q2 = [log10(p2(1)) p2(2) p2(4) min(p2(3)/p2(4), S.fMax)];
[F.W.p, ~, F.W.sse] = fit_shp(h, th, w, mdl, 'P2', [q2; start_set('P2', th, dense)], S, opts, optsLC);
end

function [p, qBest, sseBest] = fit_shp(h, th, w, mdl, pr, Q0, S, opts, optsLC)
% Least-squares fit of a retention model to theta(h) (Eq. 22) with lsqcurvefit; w = weights.
% q = [log10(alpha) n theta_r theta_s] (P1, P2x) or [log10(alpha) n theta_s f] (P2).
switch pr
    case 'P1',  lb = [-5 S.nMin 0 0];    ub = [0 10 1 1];
    case 'P2',  lb = [-5 S.nMin 0.01 0]; ub = [0 10 1 S.fMax];
    case 'P2x', lb = [-5 S.nMin 0 0.01]; ub = [0 10 1 1];
end
if isempty(w), sw = ones(size(th)); else, sw = sqrt(w); end
fun = @(q, hh) sw.*shp_theta(hh, q2p(q, pr), mdl, S);
y = sw.*th;
ws = warning;                    % interior-point KKT solves warn near the bounds of n
warning('off', 'MATLAB:nearlySingularMatrix'); warning('off', 'MATLAB:singularMatrix');
sseBest = Inf; qBest = Q0(1, :);
for j = 1:size(Q0, 1)
    q0 = min(max(Q0(j, :), lb + 1e-9), ub - 1e-9);
    try
        if strcmp(pr, 'P2x')     % theta_r - fMax*theta_s <= 0
            [q, sse] = lsqcurvefit(fun, q0, h, y, lb, ub, [0 0 1 -S.fMax], 0, [], [], [], optsLC);
            if q(3) >= q(4), continue, end
        else
            [q, sse] = lsqcurvefit(fun, q0, h, y, lb, ub, opts);
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
% the 30 starting points of Script S2; dense adds an 80-point grid
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
    case 'P1',  Q0 = [la nn rr t0];          % theta_r = rr
    case 'P2',  Q0 = [la nn t0 rr];          % f = rr
    case 'P2x', Q0 = [la nn rr.*t0 t0];      % theta_r = rr*theta_s, feasible start
end
end

function p = q2p(q, pr)
if strcmp(pr, 'P2'), p = [10^q(1) q(2) q(4)*q(3) q(3)];
else,                p = [10^q(1) q(2) q(3) q(4)];
end
end

function th = shp_theta(h, p, mdl, S)
% VG (van Genuchten, 1980), PDI (Iden and Durner, 2014) and EVM (Eq. 7) retention functions
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
% PDI non-capillary saturation: Snc = 1 + [x - xa + b ln(1 + exp((xa - x)/b))]/(xa - x0),
% x = log10 h, xa = log10 ha (ha = 1/alpha), x0 = log10 h0, written without overflow
al = p(1); n = p(2); thr = p(3); ths = p(4);
if ths == thr, r2 = Inf; else, r2 = (thr/(ths - thr))^2; end
b = 0.1 + 0.2/n^2*(1 - exp(-r2));
u = log10(al*h);                                  % x - xa
Snc = 1 - (max(u, 0) + b*log1p(exp(-abs(u)/b)))/log10(al*h0);
end

function [lK, share] = shp_logK(h, p, mdl, Ksc, S)
% log10 of the predicted K (cm/d) and the share of the non-capillary terms, (Knc + Kv)/K.
% The capillary term is evaluated in logarithms, so that it does not underflow as n -> 1.
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n;
lG = -m*log1p((al*h).^n);                                    % ln Gamma(h)
if strcmp(mdl, 'VG')
    lG0 = -Inf; lSc = lG;                                    % Se = Gamma, no h0 scaling
else
    lG0 = -m*log1p((al*S.h0)^n);
    lSc = lG + log1p(-exp(lG0 - lG)) - log1p(-exp(lG0));     % ln[(Gamma - Gamma0)/(1 - Gamma0)]
end
le = lG/m; ld = lG0/m;                                       % ln Gamma^(1/m), ln Gamma0^(1/m)
u  = (exp(le) - exp(ld))/(1 - exp(ld));
lB = log(-expm1(m*log1p(-u)));                               % ln{1 - [(1 - G^(1/m))/(1 - G0^(1/m))]^m}
t  = u < 1e-10;                                              % B = m u for small u
lB(t) = log(m) + le(t) + log1p(-exp(ld - le(t))) - log1p(-exp(ld));
lKc = log(Ksc) + 0.5*lSc + 2*lB;                             % Eq. (14); VG-Mualem for VG
if strcmp(mdl, 'VG')
    lK = lKc/log(10); share = zeros(size(h));
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
thm = shp_theta(1e5, p, mdl, S);                             % theta at h = 1e5 cm, Eq. (17)
Knc = 86.2*thm*al^1.5*(al*S.h0).^(-1.5*(1 - Snc));           % Eqs. (16)-(17); Eq. (20) with l = 0
Kv  = vapor_K(h, (ths - thr)*Sc + thr*Snc, ths);             % Appendix A
K   = exp(lKc) + Knc + Kv;
lK  = log10(K); share = (Knc + Kv)./K;
end

function Kv = vapor_K(h, th, ths)
% isothermal vapour conductivity (cm/d), Appendix A, 20 C
T = 293.15; Mw = 0.018015; g = 9.81; Rg = 8.314;
rho_sv = 1e-3*exp(31.3716 - 6014.79/T - 7.92495e-3*T)/T;     % kg m-3
rho_w  = 998.2;                                              % kg m-3
Da  = 2.12e-5*(T/273.15)^2*1e4*86400;                        % cm2 d-1
cH  = Mw*g/(Rg*T)/100;                                       % cm-1
tau = max(ths - th, 0).^(10/3)/ths^2;                        % Eq. (A4)
Kv  = rho_sv/rho_w*tau*Da*cH.*exp(-cH*h);                    % Eqs. (A1)-(A3)
end

function [nInc, hMax] = monotonic_check(p, mdl, S)
% grid intervals on 1e-12 cm <= h < h0 where theta increases with suction (C < 0)
inc = diff(shp_theta(S.hg, p, mdl, S)) > 1e-13;
nInc = sum(inc); hMax = NaN;
if nInc > 0, hMax = S.hg(find(inc, 1, 'last') + 1); end
end

function c = range_class(pF, S)
c = 1 + (pF >= S.pFedge(1)) + (pF >= S.pFedge(2));      % 1 wet, 2 intermediate, 3 dry
end

function [rm, me] = range_err(res, cls)
% RMSE and mean error per range (1-3) and over all observations (4)
rm = nan(1, 4); me = nan(1, 4);
for k = 1:3
    s = cls == k;
    if any(s), rm(k) = sqrt(mean(res(s).^2)); me(k) = mean(res(s)); end
end
if ~isempty(res), rm(4) = sqrt(mean(res.^2)); me(4) = mean(res); end
end

function s = soil_subset(sb, grp, dev)
switch sb
    case "All",     s = true(size(grp));
    case "NoDev12", s = ~dev;
    case "Dev12",   s = dev;
    otherwise,      s = grp == sb;
end
end

function g = getG(G, M, pr, sb, v, rg)
g = G(G.Model == M & G.Procedure == pr & G.Subset == sb & G.Variable == v & G.Range == rg, :);
end

function print_report(f, S, T, G, C, Q, nS, chk)
mn = S.models; rl = ["wet (pF < 2)" "intermediate (2-4.2)" "dry (pF >= 4.2)" "all"];
fprintf(f, '\n========== Validation by suction range: VG, PDI and EVM, German database (%d soils) ==========\n', nS);
fprintf(f, 'Ranges      : wet pF < %.1f | intermediate %.1f <= pF < %.1f | dry pF >= %.1f (h = 10^pF cm)\n', ...
    S.pFedge(1), S.pFedge(1), S.pFedge(2), S.pFedge(2));
fprintf(f, 'Fits (P2)   : Eq. (22) unweighted; log10(alpha) in [-5, 0], n in [%.4f, 10], theta_r = f*theta_s, f in [0, %.4f], theta_s in [0.01, 1]\n', ...
    S.nMin, S.fMax);
fprintf(f, 'Sensitivity : W = P2 with range-balanced weights (each range with data has the same total weight)\n');
fprintf(f, 'K           : Ks,c = measured KSAT (not scored); VG capillary only; PDI and EVM Kc + Knc + Kv; log10 units\n');
fprintf(f, 'Statistics  : mean (median) of per-soil RMSE over soils with data in the range [N soils]; pooled RMSE and mean error (bias = predicted - measured) over all observations\n');

fprintf(f, '\n--- Data coverage (observations) ---\n');
for i = 1:height(C)
    fprintf(f, ['   %-5s %3d soils | retention %6d: wet %6d, intermediate %6d, dry %5d (%.2f%%; %3d soils with dry data)\n' ...
        '         conductivity %6d (%3d soils): wet %5d, intermediate %6d, dry %d | pF >= 3: %4d, pF >= 3.5: %3d, max pF %.2f, AirEntry %3d | ' ...
        'soils with KSAT %3d (%5d K points)\n'], C.Group(i), C.Nsoils(i), C.NobsRet(i), C.NobsRet_wet(i), C.NobsRet_int(i), ...
        C.NobsRet_dry(i), 100*C.NobsRet_dry(i)/C.NobsRet(i), C.NsoilsRetDry(i), C.NobsK(i), C.NsoilsK(i), C.NobsK_wet(i), ...
        C.NobsK_int(i), C.NobsK_dry(i), C.NobsK_pF3(i), C.NobsK_pF35(i), C.maxpF_K(i), C.NobsK_AirEntry(i), C.NsoilsKSAT(i), ...
        C.NobsK_KSATsoils(i));
end

for pr = ["P2" "W"]
    if pr == "P2", ttl = 'Table 4 (P2: unweighted, theta_s > theta_r >= 0)'; else, ttl = 'Sensitivity W (range-balanced weights)'; end
    fprintf(f, '\n--- %s, all soils ---\n', ttl);
    for v = ["theta" "logK"]
        if v == "theta", fprintf(f, '   Retention, RMSE_theta (cm3 cm-3)\n'); else, fprintf(f, '   Conductivity, RMSE of log10 K (soils with KSAT)\n'); end
        for k = 1:4
            if v == "logK" && k == 3 && getG(G, mn(1), pr, "All", v, S.rn(k)).Nobs == 0
                fprintf(f, '      %-21s no observations\n', rl(k)); continue
            end
            g1 = getG(G, mn(1), pr, "All", v, S.rn(k));
            fprintf(f, '      %-21s N = %6d obs |', rl(k), g1.Nobs);
            for M = mn
                g = getG(G, M, pr, "All", v, S.rn(k));
                if v == "theta"
                    fprintf(f, ' %s: %.4f (%.4f) [%3d], pooled %.4f, bias %+.4f |', M, g.meanRMSE, g.medianRMSE, g.Nsoils, g.pooledRMSE, g.pooledME);
                else
                    fprintf(f, ' %s: %.3f (%.3f) [%3d], pooled %.3f, bias %+.3f |', M, g.meanRMSE, g.medianRMSE, g.Nsoils, g.pooledRMSE, g.pooledME);
                end
            end
            fprintf(f, '\n');
        end
    end
end

fprintf(f, '\n--- Retention by soil group (P2): mean per-soil RMSE_theta, wet / intermediate / dry, and dry-range bias ---\n');
for sb = ["TUDD" "UFZ" "FSUJ" "KIT" "TUBS"]
    fprintf(f, '   %-5s', sb);
    for M = mn
        a = arrayfun(@(k) getG(G, M, "P2", sb, "theta", S.rn(k)).meanRMSE, 1:3);
        b = getG(G, M, "P2", sb, "theta", "dry");
        fprintf(f, ' | %s %.4f / %.4f / %.4f, bias %+.4f', M, a, b.pooledME);
    end
    fprintf(f, '\n');
end
fprintf(f, '\n--- Conductivity by soil group (P2): mean per-soil RMSE of log10 K, wet / intermediate, and bias ---\n');
for sb = ["TUDD" "UFZ" "FSUJ" "KIT" "TUBS"]
    fprintf(f, '   %-5s', sb);
    for M = mn
        a = arrayfun(@(k) getG(G, M, "P2", sb, "logK", S.rn(k)).meanRMSE, [1 2 4]);
        b = getG(G, M, "P2", sb, "logK", "all");
        fprintf(f, ' | %s %.3f / %.3f (all %.3f), bias %+.3f', M, a, b.pooledME);
    end
    fprintf(f, '\n');
end

fprintf(f, '\n--- Dry-range retention bias (P2) by fitted theta_r ---\n');
for M = mn
    thr0 = T.(char(M + "_P2_theta_r")) <= 1e-6; me = T.(char(M + "_P2_MEth_dry")); rm = T.(char(M + "_P2_RMSEth_dry"));
    fprintf(f, '   %-3s theta_r = 0: %3d soils, dry RMSE %.4f, mean error %+.4f | theta_r > 0: %3d soils, dry RMSE %.4f, mean error %+.4f\n', M, ...
        sum(thr0 & isfinite(me)), mean(rm(thr0), 'omitnan'), mean(me(thr0), 'omitnan'), sum(~thr0 & isfinite(me)), ...
        mean(rm(~thr0), 'omitnan'), mean(me(~thr0), 'omitnan'));
end

fprintf(f, '\n--- Conductivity observations in the non-capillary (film + vapour) dominated range, (Knc + Kv)/K >= 0.5 ---\n');
nK = sum(isfinite(T.Ksat_cm_d) & (T.N_K_wet + T.N_K_int + T.N_K_dry) > 0);
for M = ["PDI" "EVM"]
    for pr = ["P2" "W"]
        nc = arrayfun(@(k) getG(G, M, pr, "All", "logK", S.rn(k)).nNonCapDominated, 1:4);
        no = arrayfun(@(k) getG(G, M, pr, "All", "logK", S.rn(k)).Nobs, 1:4);
        ns = sum(T.(char(M + "_" + pr + "_nNonCapK_all")) > 0, 'omitnan');
        fprintf(f, '   %-3s %-2s: %5d of %5d K observations (%.1f%%) in %3d of %3d soils: wet %4d of %4d, intermediate %5d of %5d, dry %d\n', ...
            M, pr, nc(4), no(4), 100*nc(4)/no(4), ns, nK, nc(1), no(1), nc(2), no(2), nc(3));
    end
    px = T.(char(M + "_P2_pFcross_K")); mx = T.maxpF_K; k = isfinite(T.Ksat_cm_d) & isfinite(mx);
    nn = T.(char(M + "_P2_n")); thr = T.(char(M + "_P2_theta_r")); hit = T.(char(M + "_P2_nNonCapK_all")) > 0;
    fprintf(f, '   %-3s P2: crossover pF (Knc + Kv = Kc) quantiles 0/5/25/50/75/95/100%%: %s; K data beyond the crossover in %d of %d soils\n', ...
        M, mat2str(quantile(px(k), [0 .05 .25 .5 .75 .95 1]), 3), sum(k & mx >= px), sum(k));
    fprintf(f, '        soils with non-capillary dominated K data: n < 1.15 (capillary term vanishes as m -> 0): %d (theta_r = 0: %d) | n > 3 (steep, coarse): %d | other: %d\n', ...
        sum(hit & nn < 1.15), sum(hit & nn < 1.15 & thr <= 1e-6), sum(hit & nn > 3), sum(hit & nn >= 1.15 & nn <= 3));
    s = Q.(char(M)).P2.sh >= 0.5; e = Q.(char(M)).P2.K;
    fprintf(f, '        P2 log10 K error where non-capillary dominated: RMSE %.3f, bias %+.3f (N %d) | capillary dominated: RMSE %.3f, bias %+.3f (N %d)\n', ...
        sqrt(mean(e(s).^2)), mean(e(s)), sum(s), sqrt(mean(e(~s).^2)), mean(e(~s)), sum(~s));
end
fprintf(f, '   Max pF of the K data per soil (soils with K data): quantiles 0/5/25/50/75/95/100%%: %s\n', ...
    mat2str(quantile(T.maxpF_K(isfinite(T.maxpF_K)), [0 .05 .25 .5 .75 .95 1]), 3));

fprintf(f, '\n--- Overlap with the 12-soil development set: the 7 German Table 1 soils (Dev12) vs the others (NoDev12), P2 ---\n');
for sb = ["NoDev12" "Dev12"]
    fprintf(f, '   %-7s', sb);
    for M = mn
        a = getG(G, M, "P2", sb, "theta", "all"); b = getG(G, M, "P2", sb, "theta", "dry"); c = getG(G, M, "P2", sb, "logK", "all");
        fprintf(f, ' | %s: RMSE_theta %.5f [%d], dry %.4f, RMSE_K %.3f [%d]', M, a.meanRMSE, a.Nsoils, b.meanRMSE, c.meanRMSE, c.Nsoils);
    end
    fprintf(f, '\n');
end
for M = mn
    a = getG(G, M, "P2", "All", "theta", "all"); b = getG(G, M, "P2", "NoDev12", "theta", "all");
    c = getG(G, M, "P2", "All", "logK", "all");  d = getG(G, M, "P2", "NoDev12", "logK", "all");
    fprintf(f, '   %-3s: excluding Dev12 changes mean RMSE_theta by %+.2e (%+.2f%%) and mean RMSE_K by %+.4f (%+.2f%%)\n', M, ...
        b.meanRMSE - a.meanRMSE, 100*(b.meanRMSE/a.meanRMSE - 1), d.meanRMSE - c.meanRMSE, 100*(d.meanRMSE/c.meanRMSE - 1));
end
ii = find(T.Dev12)';
for i = ii
    fprintf(f, '   %-9s N theta %3d (dry %d), N K %3d (max pF %.2f) | EVM RMSE_theta wet/int/dry %.4f/%.4f/%.4f, RMSE_K wet/int %.3f/%.3f\n', ...
        T.Sample_ID(i), T.N_th_wet(i) + T.N_th_int(i) + T.N_th_dry(i), T.N_th_dry(i), T.N_K_wet(i) + T.N_K_int(i), T.maxpF_K(i), ...
        T.EVM_P2_RMSEth_wet(i), T.EVM_P2_RMSEth_int(i), T.EVM_P2_RMSEth_dry(i), T.EVM_P2_RMSEK_wet(i), T.EVM_P2_RMSEK_int(i));
end

msG = ["TUDD" "UFZ" "FSUJ" "KIT" "TUBS" "All"];
fprintf(f, '\n--- Table 3: mean RMSE_theta (cm3 cm-3) / RMSE_K (log10 units) by group, P2 ---\n');
for j = 1:numel(msG)
    sel = T.Group == msG(j) | msG(j) == "All";
    fprintf(f, '   %-5s N = %3d (NK = %3d)', msG(j), sum(sel), sum(sel & isfinite(T.EVM_P2_RMSEK_all)));
    for M = ["VG" "PDI" "EVM"]
        fprintf(f, ' | %s %.4f / %.3f', M, mean(T.(char(M + "_P2_RMSEth_all"))(sel)), ...
            mean(T.(char(M + "_P2_RMSEK_all"))(sel), 'omitnan'));
    end
    fprintf(f, '\n');
end

fprintf(f, '\n--- Fit diagnostics ---\n');
for M = mn
    nn = T.(char(M + "_P2_n")); thr = T.(char(M + "_P2_theta_r")); rk = T.(char(M + "_P2_RMSEK_all"));
    fprintf(f, '   %-3s P2: n = %.4f in %d, n = 10 in %d, theta_r = 0 in %d | RMSE_K > 2 in %d soils, max %.1f\n', ...
        M, S.nMin, sum(nn <= S.nMin + 1e-6), sum(nn >= 10 - 1e-6), sum(thr <= 1e-6), sum(rk > 2), max(rk));
    fprintf(f, '        theta increasing somewhere on 1e-12 cm <= h < h0 (C < 0): P2 %d fits, W %d\n', ...
        sum(T.(char(M + "_P2_nIncrease")) > 0), sum(T.(char(M + "_W_nIncrease")) > 0));
end
if chk.n > 0
    fprintf(f, '   EVM P2 vs Script S2 (%d soils): max relative SSE difference %.1e, max |RMSE_K difference| %.1e\n', chk.n, chk.dSSE, chk.dRMSEK);
end
fprintf(f, '   fitting time per soil (three models, all fits): median %.1f s, max %.1f s\n', median(T.tFit_s), max(T.tFit_s));
end

function make_figures(outDir, S, T, Q, D, ids)
col = [0.20 0.45 0.80; 0.85 0.45 0.10; 0.15 0.60 0.30];          % VG, PDI, EVM
shade = [0.93 0.96 1.00; 0.97 0.97 0.97; 1.00 0.94 0.90];        % wet, intermediate, dry
mn = S.models;

% Figure S8: coverage of the data by suction range
f9 = figure('Visible', 'off', 'Position', [100 100 1400 900]);
ed = -2.05:0.1:6.55; ctr = ed(1:end-1) + 0.05;
mt = vertcat(D.mt); pf = vertcat(D.pF); mtK = vertcat(D.mtK); pfK = vertcat(D.pFK);
subplot(2, 2, 1); hold on
range_patches(S, shade, [0.5 2e4]);
for t = ["HyProp" "AirEntry" "WP4"]
    c = histcounts(pf(mt == t), ed); c(c == 0) = NaN;
    stairs(ctr - 0.05, c, 'LineWidth', 1.4, 'DisplayName', sprintf('%s (%d)', t, sum(mt == t)));
end
set(gca, 'YScale', 'log'); xlim([-2.1 6.6]); ylim([0.5 2e4]); grid on; box on
xlabel('pF = log_{10}(h/cm)'); ylabel('observations per 0.1 pF');
legend('Location', 'northeast'); title(sprintf('(a) retention: %d observations', numel(pf)));
subplot(2, 2, 2); hold on
range_patches(S, shade, [0.5 2e4]);
for t = ["HyProp" "AirEntry"]
    c = histcounts(pfK(mtK == t), ed); c(c == 0) = NaN;
    stairs(ctr - 0.05, c, 'LineWidth', 1.4, 'DisplayName', sprintf('%s (%d)', t, sum(mtK == t)));
end
nKs = sum(isfinite([D.Ksat]));
plot(-2, nKs, 'ks', 'MarkerFaceColor', 'k', 'DisplayName', sprintf('KSAT (%d, input K_{s,c})', nKs));
xline(max(pfK), 'r--', sprintf('max pF %.2f', max(pfK)), 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
set(gca, 'YScale', 'log'); xlim([-2.1 6.6]); ylim([0.5 2e4]); grid on; box on
xlabel('pF = log_{10}(h/cm)'); ylabel('observations per 0.1 pF');
legend('Location', 'northeast'); title(sprintf('(b) conductivity: %d observations, none at pF \\geq %.1f', numel(pfK), S.pFedge(2)));
subplot(2, 2, 3); hold on
range_patches(S, shade, [0 1]);
k = isfinite(T.Ksat_cm_d) & isfinite(T.maxpF_K);
x = sort(T.maxpF_K(k)); stairs(x, (1:numel(x))'/numel(x), 'k-', 'LineWidth', 1.6, 'DisplayName', 'driest K observation');
for M = ["PDI" "EVM"]
    x = T.(char(M + "_P2_pFcross_K"))(k); x(~isfinite(x)) = 7.5; x = sort(x);
    stairs(x, (1:numel(x))'/numel(x), '-', 'Color', col(mn == M, :), 'LineWidth', 1.4, ...
        'DisplayName', sprintf('%s crossover, K_{nc} + K_v = K_c', M));
end
xlim([-2.1 7]); ylim([0 1]); grid on; box on
xlabel('pF'); ylabel(sprintf('fraction of the %d soils with K_{s,c}', sum(k)));
legend('Location', 'southeast'); title('(c) driest K observation and capillary/non-capillary crossover');
subplot(2, 2, 4); hold on
range_patches(S, shade, [0 1]);
eb = -0.1:0.2:4.1; cb = eb(1:end-1) + 0.1;
for M = ["PDI" "EVM"]
    s = Q.(char(M)).P2.sh; bin = discretize(Q.pFK, eb);
    fr = accumarray(bin(isfinite(bin)), double(s(isfinite(bin)) >= 0.5), [numel(cb) 1], @mean, NaN);
    nb = accumarray(bin(isfinite(bin)), 1, [numel(cb) 1]); fr(nb < 20) = NaN;
    plot(cb, fr, '-o', 'Color', col(mn == M, :), 'MarkerSize', 4, 'LineWidth', 1.4, ...
        'DisplayName', sprintf('%s: %d of %d observations', M, sum(s >= 0.5), numel(s)));
end
xlim([-2.1 7]); ylim([0 1]); grid on; box on
xlabel('pF'); ylabel('fraction with (K_{nc} + K_v)/K \geq 0.5');
legend('Location', 'northwest'); title('(d) K observations dominated by the non-capillary terms');
exportgraphics(f9, fullfile(outDir, 'FigS8_range_coverage.png'), 'Resolution', 300); close(f9);

% Figure S9: errors by suction range
f10 = figure('Visible', 'off', 'Position', [100 100 1500 900]);
eb = -2.25:0.25:6.5; cb = eb(1:end-1) + 0.125;
subplot(2, 3, 1); hold on
range_patches(S, shade, [-0.08 0.08]);
for M = mn
    [md, lo, hi] = binned(Q.pF, Q.(char(M)).P2.th, eb, 20);
    plot(cb, md, '-', 'Color', col(mn == M, :), 'LineWidth', 1.6, 'DisplayName', char(M));
    plot(cb, lo, ':', cb, hi, ':', 'Color', col(mn == M, :), 'HandleVisibility', 'off');
end
yline(0, 'k-', 'HandleVisibility', 'off'); xlim([-2.1 6.6]); ylim([-0.08 0.08]); grid on; box on
xlabel('pF'); ylabel('\theta_{fit} - \theta_{meas} (cm^3 cm^{-3})');
legend('Location', 'northwest'); title('(a) retention residuals: median, 10-90%');
subplot(2, 3, 2);
box_by_range(T, S, "P2", "RMSEth", col, 1:3);
set(gca, 'YScale', 'log'); ylabel('per-soil RMSE_\theta (cm^3 cm^{-3})'); title('(b) retention RMSE by range');
subplot(2, 3, 3);
v = zeros(3, 6); lab = strings(1, 6); j = 0;
for M = mn
    for pr = ["P2" "W"]
        j = j + 1; lab(j) = M + " " + pr;
        v(:, j) = arrayfun(@(k) sqrt(mean(Q.(char(M)).(char(pr)).th(Q.cls == k).^2)), 1:3);
    end
end
b = bar(categorical(["wet" "intermediate" "dry"], ["wet" "intermediate" "dry"]), v);
for j = 1:6
    b(j).FaceColor = col(ceil(j/2), :);
    if mod(j, 2) == 0, b(j).FaceAlpha = 0.45; end
end
legend(b, lab, 'Location', 'northwest'); grid on; ylabel('pooled RMSE_\theta (cm^3 cm^{-3})');
title('(c) unweighted (P2) vs range-balanced (W) fits');
subplot(2, 3, 4); hold on
range_patches(S, shade, [-3 3]);
for M = mn
    [md, lo, hi] = binned(Q.pFK, Q.(char(M)).P2.K, eb, 20);
    plot(cb, md, '-', 'Color', col(mn == M, :), 'LineWidth', 1.6, 'DisplayName', char(M));
    plot(cb, lo, ':', cb, hi, ':', 'Color', col(mn == M, :), 'HandleVisibility', 'off');
end
yline(0, 'k-', 'HandleVisibility', 'off'); xlim([-2.1 6.6]); ylim([-3 3]); grid on; box on
xlabel('pF'); ylabel('log_{10} K_{pred} - log_{10} K_{meas}');
legend('Location', 'northwest'); title('(d) conductivity residuals (no data at pF \geq 4.2)');
subplot(2, 3, 5);
box_by_range(T, S, "P2", "RMSEK", col, 1:2);
set(gca, 'YScale', 'log'); ylabel('per-soil RMSE of log_{10} K'); title('(e) conductivity RMSE by range');
subplot(2, 3, 6); hold on
range_patches(S, shade, [-3 3]);
s = Q.EVM.P2.sh >= 0.5; e = Q.EVM.P2.K;
[md1, lo1, hi1] = binned(Q.pFK(~s), e(~s), eb, 20); [md2, lo2, hi2] = binned(Q.pFK(s), e(s), eb, 20);
plot(cb, md1, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.6, 'DisplayName', sprintf('capillary dominated (%d)', sum(~s)));
plot(cb, lo1, ':', cb, hi1, ':', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
plot(cb, md2, '-', 'Color', col(3, :), 'LineWidth', 1.6, 'DisplayName', sprintf('non-capillary dominated (%d)', sum(s)));
plot(cb, lo2, ':', cb, hi2, ':', 'Color', col(3, :), 'HandleVisibility', 'off');
yline(0, 'k-', 'HandleVisibility', 'off'); xlim([-2.1 6.6]); ylim([-3 3]); grid on; box on
xlabel('pF'); ylabel('log_{10} K_{pred} - log_{10} K_{meas}');
legend('Location', 'northwest'); title('(f) EVM residuals by dominant K term');
exportgraphics(f10, fullfile(outDir, 'FigS9_range_errors.png'), 'Resolution', 300); close(f10);

% Figure S10: example soils
ex = S.examples(ismember(S.examples, ids));
if ~isempty(ex)
    f11 = figure('Visible', 'off', 'Position', [100 100 1500 860]);
    hh = logspace(-2, log10(S.h0) - 1e-9, 600)';
    for j = 1:numel(ex)
        i = find(ids == ex(j), 1); d = D(i);
        subplot(2, numel(ex), j); hold on
        range_patches(S, shade, [0 0.75]);
        dr = d.mt == "WP4";
        plot(d.pF(~dr), d.th(~dr), '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 6, 'DisplayName', 'HYPROP and air entry');
        plot(d.pF(dr), d.th(dr), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'WP4');
        for M = mn
            p = pget(T, M, "P2", i);
            plot(log10(hh), shp_theta(hh, p, char(M), S), '-', 'Color', col(mn == M, :), 'LineWidth', 1.4, ...
                'DisplayName', sprintf('%s: dry RMSE %.3f', M, T.(char(M + "_P2_RMSEth_dry"))(i)));
        end
        xlim([-2 7]); ylim([0 max(d.th)*1.15]); grid on; box on
        xlabel('pF'); ylabel('\theta (cm^3 cm^{-3})'); legend('Location', 'southwest');
        title(sprintf('(%s) %s: retention (fitted)', char('a' + j - 1), strrep(ex(j), "_", "-")));
        subplot(2, numel(ex), j + numel(ex)); hold on
        range_patches(S, shade, [1e-12 1e5]);
        if isfinite(d.Ksat)
            plot(-2, d.Ksat, 'k^', 'MarkerFaceColor', 'k', 'DisplayName', 'KSAT = K_{s,c}');
            plot(d.pFK, d.K, 'ko', 'MarkerSize', 3, 'DisplayName', 'measured K');
            for M = mn
                p = pget(T, M, "P2", i);
                plot(log10(hh), 10.^shp_logK(hh, p, char(M), d.Ksat, S), '-', 'Color', col(mn == M, :), 'LineWidth', 1.4, ...
                    'DisplayName', sprintf('%s: RMSE %.2f', M, T.(char(M + "_P2_RMSEK_all"))(i)));
            end
            p = pget(T, "EVM", "P2", i);
            [lK, sh] = shp_logK(hh, p, 'EVM', d.Ksat, S);
            plot(log10(hh), 10.^lK.*(1 - sh), '--', 'Color', [0.4 0.4 0.4], 'DisplayName', 'EVM K_c');
            plot(log10(hh), 10.^lK.*sh, ':', 'Color', col(3, :), 'LineWidth', 1.6, 'DisplayName', 'EVM K_{nc} + K_v');
        end
        set(gca, 'YScale', 'log'); xlim([-2 7]); ylim([1e-12 1e5]); grid on; box on
        xlabel('pF'); ylabel('K (cm d^{-1})'); legend('Location', 'southwest');
        title(sprintf('(%s) %s: conductivity (predicted)', char('a' + numel(ex) + j - 1), strrep(ex(j), "_", "-")));
    end
    exportgraphics(f11, fullfile(outDir, 'FigS10_range_examples.png'), 'Resolution', 300); close(f11);
end
end

function range_patches(S, shade, yl)
x = [-2.1 S.pFedge 7.1];
for k = 1:3
    patch([x(k) x(k + 1) x(k + 1) x(k)], [yl(1) yl(1) yl(2) yl(2)], shade(k, :), 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
xline(S.pFedge, 'k:', 'HandleVisibility', 'off');
end

function [md, lo, hi] = binned(x, y, ed, nmin)
bin = discretize(x, ed); ok = isfinite(bin) & isfinite(y); nb = numel(ed) - 1;
md = accumarray(bin(ok), y(ok), [nb 1], @median, NaN);
lo = accumarray(bin(ok), y(ok), [nb 1], @(v) quantile(v, 0.1), NaN);
hi = accumarray(bin(ok), y(ok), [nb 1], @(v) quantile(v, 0.9), NaN);
n = accumarray(bin(ok), 1, [nb 1]);
md(n < nmin) = NaN; lo(n < nmin) = NaN; hi(n < nmin) = NaN;
end

function box_by_range(T, S, pr, var, col, ks)
names = ["wet" "intermediate" "dry"];
x = []; y = []; c = [];
for M = S.models
    for k = ks
        v = T.(char(M + "_" + pr + "_" + var + "_" + S.rn(k))); v = v(isfinite(v));
        x = [x; repmat(names(k), numel(v), 1)]; y = [y; v]; c = [c; repmat(M, numel(v), 1)]; %#ok<AGROW>
    end
end
b = boxchart(categorical(x, names(ks)), y, 'GroupByColor', categorical(c, S.models), 'MarkerStyle', '.');
for j = 1:numel(b), b(j).BoxFaceColor = col(j, :); b(j).MarkerColor = col(j, :); end
legend('Location', 'northwest'); grid on; box on
end

function p = pget(T, M, pr, i)
p = [T.(char(M + "_" + pr + "_alpha"))(i) T.(char(M + "_" + pr + "_n"))(i) ...
     T.(char(M + "_" + pr + "_theta_r"))(i) T.(char(M + "_" + pr + "_theta_s"))(i)];
end
