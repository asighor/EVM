function R = EVM_VG_comparison(dataDir, nSoilsMax)
%EVM_VG_COMPARISON  Documented comparison of the conventional van Genuchten (VG)
% model with EVM (and PDI), using the German database of Hohenbrink et al. (2023):
% RetMeas.csv and CondMeas.csv (572 soils; 409 with measured KSAT).
%
%   R = EVM_VG_comparison              % CSV files next to this file
%   R = EVM_VG_comparison(dataDir)     % CSV files in dataDir
%   R = EVM_VG_comparison(dataDir, 20) % quick test on the first 20 soils
%
% Answers Reviewer #2 ("Include a documented conventional VG comparison"):
%   A. Models, each with four fitted parameters (alpha, n, theta_r, theta_s):
%      VG  : theta = theta_r + (theta_s - theta_r) Gamma, Gamma = [1 + (alpha h)^n]^(-m),
%            m = 1 - 1/n (van Genuchten, 1980); K = Ks,c Gamma^0.5 [1 - (1 - Gamma^(1/m))^m]^2
%            (Mualem, 1976a), capillary only
%      PDI : Iden and Durner (2014); K = Kc (Eq. 14) + Knc (Eqs. 16-17) + Kv (Appendix A)
%      EVM : Eq. (7) with beta = 1 + 0.1/n; K = Kc (Eq. 14) + Knc (Eq. 20) + Kv (Appendix A)
%   B. Constraints and fitting, identical for the three models (Eq. 22, lsqcurvefit,
%      unweighted theta residuals): log10(alpha) in [-5, 0] (alpha in cm-1), n in
%      [1.0001, 10], theta_s > theta_r >= 0 enforced as theta_r = f*theta_s (0 <= f <= 0.9999,
%      0.01 <= theta_s <= 1); 30 starting points (110 where the constraint is active),
%      aligned with the independent-bounds fit (P1) and checked with a linear-inequality
%      fit, as in Scripts S2 and S3. Ks,c = measured KSAT (not scored). Errors: RMSE_theta
%      (Eq. 23) and RMSE of log10 K (Eq. 24), per soil and by suction range (wet pF < 2,
%      intermediate 2 <= pF < 4.2, dry pF >= 4.2).
%   C. Paired comparison per soil (EVM - VG, PDI - VG, EVM - PDI): mean difference with a
%      95% bootstrap interval, wins/losses and a two-sided sign test.
%   D. Dry-range prediction test X: every model is fitted to the data at pF < 4.2 only and
%      predicts the dry-range water contents. W: fits with range-balanced weights.
%   E. Attribution of the conductivity differences to the dominant term of the EVM
%      prediction (capillary Kc vs non-capillary Knc + Kv).
%   F. Endpoint behaviour and admissibility: theta(0), theta(h0), -dtheta/dpF at h0,
%      monotonicity (C > 0) with independent bounds and with the constraint, theta_r
%      relative to the measured dry-range water contents, and log10 K at pF 5, 6 and h0.
%
% Suction h is a positive magnitude (cm), h0 = 1e7 cm (oven dryness), ha = 1/alpha.
%
% Outputs (folder vg_results next to the data):
%   EVM_VG_report.txt     all numbers quoted in the manuscript and the reply
%   EVM_VG_by_soil.csv    per soil and model: parameters, errors and endpoints (Text S5)
%   EVM_VG_paired.csv     paired comparisons by procedure, variable and range (Table S7)
%   FigS11_VG_comparison.png, FigS12_VG_endpoints.png, EVM_VG_results.mat
% Runtime: about 10 min with 24 parallel workers (R2025b, Optimization and Parallel
% Computing Toolboxes; parfor runs serially without a pool).

if nargin < 1 || isempty(dataDir)
    dataDir = fileparts(mfilename('fullpath'));
end
outDir = fullfile(dataDir, 'vg_results');
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
S.procs  = ["P2" "W" "X"];
S.pFK    = [5 6 log10(1e7) - 1e-9];                       % log10 K reported at these pF
S.dev12  = ["FSUJ_023" "FSUJ_035" "TUBS_027" "TUBS_128" "KIT_022" "KIT_009" "TUBS_106"];
S.hg     = logspace(-12, log10(S.h0) - 1e-9, 20000)';     % monotonicity check on 0 < h < h0
S.hx     = logspace(-2, log10(S.h0) - 1e-9, 4000)';       % capillary/non-capillary crossover
S.nBoot  = 2000; S.seed = 1;                               % bootstrap of the paired means
S.tieTol = 1e-9;                                           % |difference| below this is a tie
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
D = struct('id', cellstr(ids), 'pF', [], 'th', [], 'mt', [], 'pFK', [], 'K', [], 'Ksat', NaN);
for i = 1:nS
    s = ir == i;
    D(i).pF = ret.pF(s); D(i).th = ret.theta(s); D(i).mt = ret.MeasType(s);
    isKs  = ic == i & cond.MeasType == "KSAT";
    unsat = ic == i & cond.MeasType ~= "KSAT" & cond.k > 0;
    D(i).pFK = cond.pF(unsat); D(i).K = cond.k(unsat);
    if any(isKs), D(i).Ksat = cond.k(find(isKs, 1)); end
end
grp = extractBefore(ids, "_");
dev = ismember(ids, S.dev12);
fprintf('Loaded %d soils, %d retention points, %d soils with measured KSAT.\n', ...
    nS, sum(ir > 0), sum(isfinite([D.Ksat])));

if isempty(gcp('nocreate'))          % parfor runs serially when no pool can start
    try
        pc = parcluster('Processes');             % in-session change only, profile not saved
        pc.NumWorkers = max(pc.NumWorkers, feature('numcores'));
        parpool(pc, feature('numcores'));
    catch
        try, parpool('Processes'); catch, end %#ok<NOCOM>
    end
end

%% ---------------- B, D-F. fits, errors and endpoints ----------------
out = cell(nS, 1);
parfor i = 1:nS
    out{i} = analyse_soil(D(i), S, opts, optsLC);
end
fprintf('Fits finished (%.1f min).\n', toc(tStart)/60);

T = struct2table(cell2mat(cellfun(@(o) o.row, out, 'UniformOutput', false)));
T = [table(ids, grp, dev, [D.Ksat]', 'VariableNames', {'Sample_ID', 'Group', 'Dev12', 'Ksat_cm_d'}), T];

% pooled observations for the figures and the attribution of the K differences
Q.pF   = cell2mat(cellfun(@(o) o.pF, out, 'UniformOutput', false));
Q.cls  = cell2mat(cellfun(@(o) o.cls, out, 'UniformOutput', false));
Q.pFK  = cell2mat(cellfun(@(o) o.pFK, out, 'UniformOutput', false));
for M = cellstr(S.models)
    m = M{1};
    Q.(m).P2th = cell2mat(cellfun(@(o) o.res.(m).P2th, out, 'UniformOutput', false));
    Q.(m).Xth  = cell2mat(cellfun(@(o) o.res.(m).Xth,  out, 'UniformOutput', false));
    Q.(m).K    = cell2mat(cellfun(@(o) o.res.(m).K,    out, 'UniformOutput', false));
    Q.(m).sh   = cell2mat(cellfun(@(o) o.res.(m).sh,   out, 'UniformOutput', false));
end

%% ---------------- C. paired comparisons ----------------
pairs = ["EVM" "VG"; "PDI" "VG"; "EVM" "PDI"];
rows = cell(0, 14);
for j = 1:size(pairs, 1)
    for pr = S.procs
        for v = ["theta" "logK"]
            for k = 1:4
                if v == "logK" && k == 3, continue, end          % no conductivity data in the dry range
                if v == "theta", nm = "_RMSEth_"; else, nm = "_RMSEK_"; end
                a = T.(char(pairs(j, 1) + "_" + pr + nm + S.rn(k)));
                b = T.(char(pairs(j, 2) + "_" + pr + nm + S.rn(k)));
                st = paired_stats(a, b, S);
                rows(end + 1, :) = {pairs(j, 1) + " - " + pairs(j, 2), pr, v, S.rn(k), st.n, st.meanA, st.meanB, ...
                    st.meanD, st.lo, st.hi, st.medD, st.win, st.loss, st.p}; %#ok<AGROW>
            end
        end
    end
end
PT = cell2table(rows, 'VariableNames', {'Comparison', 'Procedure', 'Variable', 'Range', 'Nsoils', 'meanRMSE_first', ...
    'meanRMSE_second', 'meanDiff', 'CI95_lo', 'CI95_hi', 'medianDiff', 'firstLower', 'firstHigher', 'signTest_p'});
for v = ["Comparison" "Procedure" "Variable" "Range"], PT.(v) = string(PT.(v)); end

% consistency with Script S3 (the same P2 fits)
chk = struct('n', 0, 'dth', NaN, 'dK', NaN);
f3 = fullfile(dataDir, 'range_results', 'EVM_range_by_soil.csv');
if exist(f3, 'file')
    T3 = readtable(f3, 'TextType', 'string');
    [tf, j] = ismember(ids, T3.Sample_ID);
    if any(tf)
        chk.n = sum(tf); chk.dth = 0; chk.dK = 0;
        for M = S.models
            a = T.(char(M + "_P2_RMSEth_all"))(tf); b = T3.(char(M + "_P2_RMSEth_all"))(j(tf));
            chk.dth = max(chk.dth, max(abs(a - b)));
            a = T.(char(M + "_P2_RMSEK_all"))(tf); b = T3.(char(M + "_P2_RMSEK_all"))(j(tf));
            chk.dK = max(chk.dK, max(abs(a - b), [], 'omitnan'));
        end
    end
end

%% ---------------- outputs ----------------
vn = T.Properties.VariableNames;                % constrained fits only (P1 is internal)
Tw = T(:, ~contains(vn, '_P1_') & ~endsWith(vn, '_active') & ~endsWith(vn, '_nFix'));
writetable(Tw, fullfile(outDir, 'EVM_VG_by_soil.csv'));
writetable(PT, fullfile(outDir, 'EVM_VG_paired.csv'));
fid = fopen(fullfile(outDir, 'EVM_VG_report.txt'), 'w');
if fid < 0, error('Cannot write the report in %s', outDir); end
for f = [1, fid]
    print_report(f, S, T, PT, Q, nS, chk);
end
fclose(fid);
make_figures(outDir, S, T, Q, PT);
R = struct('settings', S, 'bySoil', Tw, 'paired', PT, 'pooled', Q);
save(fullfile(outDir, 'EVM_VG_results.mat'), 'R');
fprintf('Outputs written to %s (%.1f min).\n', outDir, toc(tStart)/60);
fprintf('EVM_VG_COMPARISON_DONE\n');
end

%% =====================================================================
%% local functions
%% =====================================================================
function o = analyse_soil(d, S, opts, optsLC)
% Fits VG, PDI and EVM to one soil (P1, P2, W, X); errors by range and endpoint diagnostics.
t0 = tic;
h = 10.^d.pF; cls = range_class(d.pF, S);
w = zeros(size(h));                    % W: every range with data gets the same total weight
for k = 1:3
    s = cls == k;
    if any(s), w(s) = 1/sum(s); end
end
w = w*numel(h)/sum(w);
inX = cls < 3;                         % X: fit to pF < 4.2, predict the dry range
hasK = isfinite(d.Ksat) && ~isempty(d.K);
clsK = range_class(d.pFK, S);
r.N_th_wet = sum(cls == 1); r.N_th_int = sum(cls == 2); r.N_th_dry = sum(cls == 3);
r.N_K = numel(d.K);
r.minTheta_dry = NaN; r.maxpF_theta = max(d.pF);
if any(cls == 3), r.minTheta_dry = min(d.th(cls == 3)); end
o.pF = d.pF; o.cls = cls; o.pFK = zeros(0, 1);
if hasK, o.pFK = d.pFK; hK = 10.^d.pFK; lKm = log10(d.K); end
for M = S.models
    mdl = char(M);
    F = fit_procedures(h, d.th, w, inX, mdl, S, opts, optsLC);
    r.([mdl '_active']) = F.active; r.([mdl '_nFix']) = F.nFix;
    p = F.P1.p;                         % independent bounds: admissibility only
    r.([mdl '_P1_theta_r']) = p(3); r.([mdl '_P1_theta_s']) = p(4);
    r.([mdl '_P1_RMSEth_all']) = sqrt(mean((shp_theta(h, p, mdl, S) - d.th).^2));
    r.([mdl '_P1_nIncrease']) = monotonic_check(p, mdl, S);
    for pr = S.procs
        p = F.(pr).p; px = [mdl '_' char(pr) '_'];
        r.([px 'alpha']) = p(1); r.([px 'n']) = p(2); r.([px 'theta_r']) = p(3); r.([px 'theta_s']) = p(4);
        res = shp_theta(h, p, mdl, S) - d.th;
        [rm, me] = range_err(res, cls);
        for k = 1:4
            r.([px 'RMSEth_' char(S.rn(k))]) = rm(k); r.([px 'MEth_' char(S.rn(k))]) = me(k);
        end
        rm = nan(1, 4);
        if hasK
            [lK, sh] = shp_logK(hK, p, mdl, d.Ksat, S);
            rm = range_err(lK - lKm, clsK);
            if pr == "P2", o.res.(mdl).K = lK - lKm; o.res.(mdl).sh = sh; end
        elseif pr == "P2"
            o.res.(mdl).K = zeros(0, 1); o.res.(mdl).sh = zeros(0, 1);
        end
        for k = [1 2 4]
            r.([px 'RMSEK_' char(S.rn(k))]) = rm(k);
        end
        if pr == "P2", o.res.(mdl).P2th = res; end
        if pr == "X",  o.res.(mdl).Xth = res(cls == 3); end
    end
    % endpoints and admissibility of the constrained fit (P2)
    p = F.P2.p; px = [mdl '_P2_'];
    r.([px 'nIncrease']) = monotonic_check(p, mdl, S);
    r.([px 'theta0_minus_ths']) = shp_theta(0, p, mdl, S) - p(4);
    r.([px 'theta_h0']) = shp_theta(S.h0*(1 - 1e-12), p, mdl, S);
    r.([px 'slope_h0']) = (shp_theta(S.h0*10^(-1e-4), p, mdl, S) - shp_theta(S.h0*(1 - 1e-12), p, mdl, S))/1e-4;
    lK = nan(1, 3);
    if hasK, lK = shp_logK(10.^S.pFK(:), p, mdl, d.Ksat, S)'; end
    r.([px 'log10K_pF5']) = lK(1); r.([px 'log10K_pF6']) = lK(2); r.([px 'log10K_h0']) = lK(3);
end
% conductivity errors of every model at the points where the EVM prediction is capillary
% (Kc) or non-capillary (Knc + Kv) dominated
sh = o.res.EVM.sh; r.EVM_fracNonCapK = NaN;
if hasK, r.EVM_fracNonCapK = mean(sh >= 0.5); end
for M = S.models
    e = o.res.(char(M)).K; rc = NaN; rnc = NaN;
    if any(sh < 0.5),  rc  = sqrt(mean(e(sh < 0.5).^2)); end
    if any(sh >= 0.5), rnc = sqrt(mean(e(sh >= 0.5).^2)); end
    r.([char(M) '_P2_RMSEK_EVMcap']) = rc; r.([char(M) '_P2_RMSEK_EVMnoncap']) = rnc;
end
r.tFit_s = toc(t0);
o.row = r;
end

function F = fit_procedures(h, th, w, inX, mdl, S, opts, optsLC)
% P1: independent bounds; P2: theta_s > theta_r >= 0 (theta_r = f theta_s), aligned with
% P1 as in Script S2; W: P2 with range-balanced weights; X: P2 fitted to pF < 4.2 only.
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
if all(inX)
    F.X = F.P2;                      % no dry-range data: nothing to predict
else
    [F.X.p, ~, F.X.sse] = fit_shp(h(inX), th(inX), [], mdl, 'P2', [q2; start_set('P2', th(inX), dense)], ...
        S, opts, optsLC);
end
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

function nInc = monotonic_check(p, mdl, S)
% grid intervals on 1e-12 cm <= h < h0 where theta increases with suction (C < 0)
nInc = sum(diff(shp_theta(S.hg, p, mdl, S)) > 1e-13);
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

function st = paired_stats(a, b, S)
% paired difference a - b over the soils where both are finite: mean with a 95% bootstrap
% interval (soils resampled), median, wins/losses and a two-sided sign test (exact binomial)
ok = isfinite(a) & isfinite(b); a = a(ok); b = b(ok); d = a - b; n = numel(d);
st = struct('n', n, 'meanA', mean(a), 'meanB', mean(b), 'meanD', mean(d), 'lo', NaN, 'hi', NaN, ...
    'medD', median(d), 'win', sum(d < -S.tieTol), 'loss', sum(d > S.tieTol), 'p', NaN);
if n == 0, return, end
rs = RandStream('mt19937ar', 'Seed', S.seed);
mb = mean(d(randi(rs, n, n, S.nBoot)), 1);
q = quantile(mb, [0.025 0.975]); st.lo = q(1); st.hi = q(2);
nn = st.win + st.loss; kmin = min(st.win, st.loss);
st.p = 1;
if nn > 0, st.p = min(1, 2*betainc(0.5, nn - kmin, kmin + 1)); end
end

function g = getP(PT, cmp, pr, v, rg)
g = PT(PT.Comparison == cmp & PT.Procedure == pr & PT.Variable == v & PT.Range == rg, :);
end

function print_report(f, S, T, PT, Q, nS, chk)
mn = S.models; rl = ["wet (pF < 2)" "intermediate (2-4.2)" "dry (pF >= 4.2)" "all"];
fprintf(f, '\n========== Conventional van Genuchten (VG) vs PDI and EVM, German database (%d soils) ==========\n', nS);
fprintf(f, 'Models     : VG theta = thr + (ths - thr) Gamma, m = 1 - 1/n, Mualem K (capillary only); PDI and EVM Kc + Knc + Kv\n');
fprintf(f, 'Fits       : Eq. (22) unweighted; log10(alpha) in [-5, 0], n in [%.4f, 10], theta_r = f*theta_s, f in [0, %.4f], theta_s in [0.01, 1]\n', ...
    S.nMin, S.fMax);
fprintf(f, 'Procedures : P2 all data (main) | W range-balanced weights | X fitted to pF < %.1f, dry range predicted\n', S.pFedge(2));
fprintf(f, 'Paired     : first - second per soil; mean [95%% bootstrap CI, %d resamples]; first lower / higher; sign-test p\n', S.nBoot);

fprintf(f, '\n--- Mean per-soil errors (P2) ---\n');
for k = 1:4
    fprintf(f, '   RMSE_theta %-21s', rl(k));
    for M = mn, fprintf(f, ' | %s %.4f', M, mean(T.(char(M + "_P2_RMSEth_" + S.rn(k))), 'omitnan')); end
    fprintf(f, '\n');
end
for k = [1 2 4]
    fprintf(f, '   RMSE_K     %-21s', rl(k));
    for M = mn, fprintf(f, ' | %s %.3f (median %.3f)', M, mean(T.(char(M + "_P2_RMSEK_" + S.rn(k))), 'omitnan'), ...
            median(T.(char(M + "_P2_RMSEK_" + S.rn(k))), 'omitnan')); end
    fprintf(f, '\n');
end

for cmp = ["EVM - VG" "PDI - VG" "EVM - PDI"]
    fprintf(f, '\n--- Paired comparison %s ---\n', cmp);
    for pr = S.procs
        for v = ["theta" "logK"]
            for k = 1:4
                g = getP(PT, cmp, pr, v, S.rn(k));
                if isempty(g), continue, end
                if v == "theta", fm = '%+.5f [%+.5f, %+.5f], median %+.5f'; else, fm = '%+.3f [%+.3f, %+.3f], median %+.3f'; end
                fprintf(f, ['   %-2s %-5s %-21s N %3d: mean %.5g vs %.5g, diff ' fm ' | lower %3d, higher %3d, p = %.1e\n'], ...
                    pr, v, rl(k), g.Nsoils, g.meanRMSE_first, g.meanRMSE_second, g.meanDiff, g.CI95_lo, g.CI95_hi, ...
                    g.medianDiff, g.firstLower, g.firstHigher, g.signTest_p);
            end
        end
    end
end

fprintf(f, '\n--- Dry-range retention by procedure: mean per-soil RMSE_theta and mean error (soils with dry data) ---\n');
for pr = S.procs
    fprintf(f, '   %-2s', pr);
    for M = mn
        fprintf(f, ' | %s %.4f (%+.4f)', M, mean(T.(char(M + "_" + pr + "_RMSEth_dry")), 'omitnan'), ...
            mean(T.(char(M + "_" + pr + "_MEth_dry")), 'omitnan'));
    end
    fprintf(f, ' | wet RMSE: VG %.4f PDI %.4f EVM %.4f\n', mean(T.(char("VG_" + pr + "_RMSEth_wet")), 'omitnan'), ...
        mean(T.(char("PDI_" + pr + "_RMSEth_wet")), 'omitnan'), mean(T.(char("EVM_" + pr + "_RMSEth_wet")), 'omitnan'));
end
e = cell2mat(cellfun(@(m) Q.(m).Xth, cellstr(mn), 'UniformOutput', false));
fprintf(f, '   X pooled dry-range prediction RMSE: VG %.4f, PDI %.4f, EVM %.4f (%d observations)\n', ...
    sqrt(mean(e(:, 1).^2)), sqrt(mean(e(:, 2).^2)), sqrt(mean(e(:, 3).^2)), size(e, 1));

fprintf(f, '\n--- Conductivity by the dominant term of the EVM prediction (P2, pooled observations) ---\n');
s = Q.EVM.sh >= 0.5;
for cl = [false true]
    if cl, lab = 'Knc + Kv > Kc'; sel = s; else, lab = 'Kc dominant  '; sel = ~s; end
    fprintf(f, '   %s (%5d obs):', lab, sum(sel));
    for M = mn
        ee = Q.(char(M)).K(sel);
        fprintf(f, ' | %s RMSE %.3f, mean %+.3f, median %+.3f', M, sqrt(mean(ee.^2)), mean(ee), median(ee));
    end
    fprintf(f, '\n');
end
k = isfinite(T.EVM_P2_RMSEK_all);
fprintf(f, '   per soil: RMSE_K of VG minus EVM, median %+.3f for soils with >= 50%% non-capillary dominated points (N %d), %+.3f for the others (N %d)\n', ...
    median(T.VG_P2_RMSEK_all(k & T.EVM_fracNonCapK >= 0.5) - T.EVM_P2_RMSEK_all(k & T.EVM_fracNonCapK >= 0.5)), ...
    sum(k & T.EVM_fracNonCapK >= 0.5), median(T.VG_P2_RMSEK_all(k & T.EVM_fracNonCapK < 0.5) - ...
    T.EVM_P2_RMSEK_all(k & T.EVM_fracNonCapK < 0.5)), sum(k & T.EVM_fracNonCapK < 0.5));

fprintf(f, '\n--- Endpoint behaviour and admissibility (P2 unless stated) ---\n');
for M = mn
    t0 = T.(char(M + "_P2_theta0_minus_ths")); th0 = T.(char(M + "_P2_theta_h0")); sl = T.(char(M + "_P2_slope_h0"));
    fprintf(f, '   %-3s theta(0) - theta_s max |%.1e| | theta(h0): median %.4f, max %.4f, > 0.01 in %d, > 0.05 in %d | -dtheta/dpF at h0: median %.2e\n', ...
        M, max(abs(t0)), median(th0), max(th0), sum(th0 > 0.01), sum(th0 > 0.05), median(sl));
    fprintf(f, '        theta increasing somewhere (C < 0): %d fits | n > 1 in all fits: %d\n', ...
        sum(T.(char(M + "_P2_nIncrease")) > 0), all(T.(char(M + "_P2_n")) > 1));
    fprintf(f, '        log10 K (cm/d) median at pF 5 / 6 / h0: %.2f / %.2f / %.2f | theta_r = 0 in %d fits, f at upper bound in %d\n', ...
        median(T.(char(M + "_P2_log10K_pF5")), 'omitnan'), median(T.(char(M + "_P2_log10K_pF6")), 'omitnan'), ...
        median(T.(char(M + "_P2_log10K_h0")), 'omitnan'), sum(T.(char(M + "_P2_theta_r")) <= 1e-6), ...
        sum(T.(char(M + "_P2_theta_r")) >= (S.fMax - 1e-6)*T.(char(M + "_P2_theta_s"))));
end
v = T.VG_P2_theta_h0; q = quantile(v, [0 .05 .25 .5 .75 .95 1]);
fprintf(f, '   VG theta(h0) quantiles 0/5/25/50/75/95/100%%: %s; theta(h0)/theta_s median %.3f\n', mat2str(q, 3), median(v./T.VG_P2_theta_s));
dd = isfinite(T.minTheta_dry) & T.VG_P2_theta_r > 1e-6;
fprintf(f, '   VG theta_r > lowest measured dry-range theta: %d of %d soils with dry data and theta_r > 0 (of %d with dry data); median excess %.4f\n', ...
    sum(dd & T.VG_P2_theta_r > T.minTheta_dry), sum(dd), sum(isfinite(T.minTheta_dry)), ...
    median(T.VG_P2_theta_r(dd & T.VG_P2_theta_r > T.minTheta_dry) - T.minTheta_dry(dd & T.VG_P2_theta_r > T.minTheta_dry)));
for pf = ["pF5" "pF6" "h0"]
    a = T.(char("VG_P2_log10K_" + pf)); b = T.(char("EVM_P2_log10K_" + pf)); k = isfinite(a) & isfinite(b);
    fprintf(f, '   log10 K at %-3s: EVM - VG median %+.1f orders of magnitude (5-95%%: %+.1f to %+.1f; N %d)\n', pf, ...
        median(b(k) - a(k)), quantile(b(k) - a(k), 0.05), quantile(b(k) - a(k), 0.95), sum(k));
end

fprintf(f, '\n--- Table 3 by group (P2): mean RMSE_theta / RMSE_K ---\n');
for g = ["TUDD" "UFZ" "FSUJ" "KIT" "TUBS" "All"]
    sel = T.Group == g | g == "All";
    fprintf(f, '   %-5s N %3d (%3d)', g, sum(sel), sum(sel & isfinite(T.EVM_P2_RMSEK_all)));
    for M = mn
        fprintf(f, ' | %s %.4f / %.3f', M, mean(T.(char(M + "_P2_RMSEth_all"))(sel)), mean(T.(char(M + "_P2_RMSEK_all"))(sel), 'omitnan'));
    end
    fprintf(f, '\n');
end

fprintf(f, '\n--- VG fits of the German Table 1 soils (for Table 2) ---\n');
for i = find(T.Dev12)'
    fprintf(f, '   %-9s VG: theta_r %.3f theta_s %.3f alpha %.4f n %.3f RMSE_theta %.4f RMSE_K %.2f | EVM %.4f / %.2f | PDI %.4f / %.2f\n', ...
        T.Sample_ID(i), T.VG_P2_theta_r(i), T.VG_P2_theta_s(i), T.VG_P2_alpha(i), T.VG_P2_n(i), T.VG_P2_RMSEth_all(i), ...
        T.VG_P2_RMSEK_all(i), T.EVM_P2_RMSEth_all(i), T.EVM_P2_RMSEK_all(i), T.PDI_P2_RMSEth_all(i), T.PDI_P2_RMSEK_all(i));
end

fprintf(f, '\n--- Checks ---\n');
if chk.n > 0
    fprintf(f, '   P2 fits vs Script S3 (%d soils): max |RMSE_theta difference| %.1e, max |RMSE_K difference| %.1e\n', chk.n, chk.dth, chk.dK);
end
fprintf(f, '   fitting time per soil (three models, all fits): median %.1f s, max %.1f s\n', median(T.tFit_s), max(T.tFit_s));
end

function make_figures(outDir, S, T, Q, PT)
col = [0.20 0.45 0.80; 0.85 0.45 0.10; 0.15 0.60 0.30];          % VG, PDI, EVM
mn = S.models;

% Figure S11: errors of VG, PDI and EVM
f12 = figure('Visible', 'off', 'Position', [100 100 1500 900]);
subplot(2, 3, 1);
x = T.VG_P2_RMSEth_all; y = T.EVM_P2_RMSEth_all; z = T.VG_P2_theta_r <= 1e-6;
lim = [min([x; y])/1.3 max([x; y])*1.3];
loglog(x(~z), y(~z), '.', 'Color', col(1, :), 'MarkerSize', 9); hold on
loglog(x(z), y(z), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 9); loglog(lim, lim, 'k--');
axis([lim lim]); axis square; grid on
xlabel('RMSE_\theta, VG (cm^3 cm^{-3})'); ylabel('RMSE_\theta, EVM (cm^3 cm^{-3})');
g = getP(PT, "EVM - VG", "P2", "theta", "all");
legend(sprintf('VG \\theta_r > 0 (%d)', sum(~z)), sprintf('VG \\theta_r = 0 (%d)', sum(z)), '1:1', 'Location', 'northwest');
title(sprintf('(a) retention: EVM lower in %d of %d soils', g.firstLower, g.Nsoils));
subplot(2, 3, 2);
x = T.VG_P2_RMSEK_all; y = T.EVM_P2_RMSEK_all; k = isfinite(x) & isfinite(y);
lim = [min([x(k); y(k)])/1.3 max([x(k); y(k)])*1.3];
scatter(x(k), y(k), 10, T.EVM_fracNonCapK(k), 'filled'); hold on
plot(lim, lim, 'k--'); set(gca, 'XScale', 'log', 'YScale', 'log'); axis([lim lim]); grid on; box on
colormap(gca, parula); cb = colorbar; cb.Label.String = 'fraction of points with K_{nc} + K_v > K_c (EVM)'; clim([0 1]);
xlabel('RMSE of log_{10} K, VG'); ylabel('RMSE of log_{10} K, EVM');
g = getP(PT, "EVM - VG", "P2", "logK", "all");
title(sprintf('(b) conductivity: EVM lower in %d of %d soils', g.firstLower, g.Nsoils));
subplot(2, 3, 3);
names = ["wet" "intermediate" "dry"]; xg = []; yg = []; cg = [];
for c = ["PDI" "EVM"]
    for k = 1:3
        v = T.(char(c + "_P2_RMSEth_" + S.rn(k))) - T.(char("VG_P2_RMSEth_" + S.rn(k))); v = v(isfinite(v));
        xg = [xg; repmat(names(k), numel(v), 1)]; yg = [yg; v]; cg = [cg; repmat(c + " - VG", numel(v), 1)]; %#ok<AGROW>
    end
end
b = boxchart(categorical(xg, names), yg, 'GroupByColor', categorical(cg, ["PDI - VG" "EVM - VG"]), 'MarkerStyle', '.');
b(1).BoxFaceColor = col(2, :); b(2).BoxFaceColor = col(3, :); b(1).MarkerColor = col(2, :); b(2).MarkerColor = col(3, :);
yline(0, 'k-', 'HandleVisibility', 'off'); ylim([-0.08 0.03]); grid on; legend('Location', 'southwest');
ylabel('\Delta RMSE_\theta (cm^3 cm^{-3})'); title('(c) per-soil difference to VG by range');
subplot(2, 3, 4);
s = Q.EVM.sh >= 0.5; cn = ["K_c dominant" "K_{nc} + K_v dominant"];
xg = []; yg = []; cg = [];
for M = mn
    e = abs(Q.(char(M)).K); lab = strings(size(e)); lab(~s) = "K_c dominant"; lab(s) = "K_{nc} + K_v dominant";
    xg = [xg; lab]; yg = [yg; e]; cg = [cg; repmat(M, numel(e), 1)]; %#ok<AGROW>
end
b = boxchart(categorical(xg, cn), max(yg, 1e-3), 'GroupByColor', categorical(cg, mn), 'MarkerStyle', 'none');
for j = 1:3, b(j).BoxFaceColor = col(j, :); end
set(gca, 'YScale', 'log'); grid on; legend('Location', 'northwest');
ylabel('|log_{10} K_{pred} - log_{10} K_{meas}|');
title(sprintf('(d) conductivity errors by dominant EVM term (%d / %d points)', sum(~s), sum(s)));
subplot(2, 3, 5);
labs = ["P2: all data" "W: range-balanced" "X: fit to pF < 4.2"]; xg = []; yg = []; cg = [];
for j = 1:3
    for M = mn
        v = T.(char(M + "_" + S.procs(j) + "_RMSEth_dry")); v = v(isfinite(v));
        xg = [xg; repmat(labs(j), numel(v), 1)]; yg = [yg; v]; cg = [cg; repmat(M, numel(v), 1)]; %#ok<AGROW>
    end
end
b = boxchart(categorical(xg, labs), yg, 'GroupByColor', categorical(cg, mn), 'MarkerStyle', '.');
for j = 1:3, b(j).BoxFaceColor = col(j, :); b(j).MarkerColor = col(j, :); end
set(gca, 'YScale', 'log'); grid on; legend('Location', 'southwest');
ylabel('dry-range RMSE_\theta (cm^3 cm^{-3})'); title('(e) dry range (pF \geq 4.2): fitted and predicted');
subplot(2, 3, 6); hold on
ed = 4.2:0.2:6.4; cb2 = ed(1:end-1) + 0.1; pfd = Q.pF(Q.cls == 3);
for j = 1:3
    e = Q.(char(mn(j))).Xth;
    [md, lo, hi] = binned(pfd, e, ed, 10);
    plot(cb2, md, '-o', 'Color', col(j, :), 'LineWidth', 1.6, 'MarkerSize', 4, 'DisplayName', ...
        sprintf('%s: RMSE %.4f', mn(j), sqrt(mean(e.^2))));
    plot(cb2, lo, ':', cb2, hi, ':', 'Color', col(j, :), 'HandleVisibility', 'off');
end
yline(0, 'k-', 'HandleVisibility', 'off'); grid on; box on; xlim([4.1 6.5]); ylim([-0.1 0.3]);
xlabel('pF'); ylabel('\theta_{pred} - \theta_{meas} (cm^3 cm^{-3})'); legend('Location', 'northwest');
title('(f) X: dry range predicted from pF < 4.2');
exportgraphics(f12, fullfile(outDir, 'FigS11_VG_comparison.png'), 'Resolution', 300); close(f12);

% Figure S12: endpoint behaviour
f13 = figure('Visible', 'off', 'Position', [100 100 1500 460]);
subplot(1, 3, 1); hold on
v = sort(T.VG_P2_theta_h0); stairs(v, (1:numel(v))'/numel(v), '-', 'Color', col(1, :), 'LineWidth', 1.6);
v = sort(T.VG_P2_theta_r); stairs(v, (1:numel(v))'/numel(v), '--', 'Color', col(1, :), 'LineWidth', 1.2);
xline(max(T.EVM_P2_theta_h0), '-', 'Color', col(3, :), 'LineWidth', 2);
grid on; box on; xlim([-0.005 0.3]); ylim([0 1]);
xlabel('\theta (cm^3 cm^{-3})'); ylabel('fraction of soils');
legend('VG: \theta(h_0)', 'VG: \theta_r', sprintf('EVM: \\theta(h_0) = 0 (PDI max %.1e)', max(T.PDI_P2_theta_h0)), 'Location', 'southeast');
title(sprintf('(a) water content at oven dryness, median VG %.3f', median(T.VG_P2_theta_h0)));
subplot(1, 3, 2);
k = isfinite(T.minTheta_dry) & T.VG_P2_theta_r > 1e-6; x = T.VG_P2_theta_r(k); y = T.minTheta_dry(k);
plot(x, y, '.', 'Color', col(1, :), 'MarkerSize', 9); hold on; plot([0 0.3], [0 0.3], 'k--');
axis([0 0.25 0 0.25]); axis square; grid on
xlabel('fitted VG \theta_r (cm^3 cm^{-3})'); ylabel('lowest measured \theta at pF \geq 4.2');
title(sprintf('(b) VG asymptote above the data in %d of %d soils', sum(y < x), numel(x)));
subplot(1, 3, 3);
labs = ["pF 5" "pF 6" "pF 7 (h_0)"]; fl = ["pF5" "pF6" "h0"]; xg = []; yg = []; cg = [];
for j = 1:3
    for M = mn
        v = T.(char(M + "_P2_log10K_" + fl(j))); v = v(isfinite(v));
        xg = [xg; repmat(labs(j), numel(v), 1)]; yg = [yg; v]; cg = [cg; repmat(M, numel(v), 1)]; %#ok<AGROW>
    end
end
nlow = sum(yg < -30 & cg == "VG");
b = boxchart(categorical(xg, labs), max(yg, -30), 'GroupByColor', categorical(cg, mn), 'MarkerStyle', '.');
for j = 1:3, b(j).BoxFaceColor = col(j, :); b(j).MarkerColor = col(j, :); end
ylim([-31 0]); grid on; legend('Location', 'southwest'); ylabel('predicted log_{10} K (cm d^{-1})');
title({'(c) predicted K where no K data exist', sprintf('(%d VG values below -30 plotted at -30)', nlow)});
exportgraphics(f13, fullfile(outDir, 'FigS12_VG_endpoints.png'), 'Resolution', 300); close(f13);
end

function [md, lo, hi] = binned(x, y, ed, nmin)
bin = discretize(x, ed); ok = isfinite(bin) & isfinite(y); nb = numel(ed) - 1;
md = accumarray(bin(ok), y(ok), [nb 1], @median, NaN);
lo = accumarray(bin(ok), y(ok), [nb 1], @(v) quantile(v, 0.1), NaN);
hi = accumarray(bin(ok), y(ok), [nb 1], @(v) quantile(v, 0.9), NaN);
n = accumarray(bin(ok), 1, [nb 1]);
md(n < nmin) = NaN; lo(n < nmin) = NaN; hi(n < nmin) = NaN;
end
