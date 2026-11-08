function R = EVM_beta_closure_analysis(dataDir, nSoilsMax)
%EVM_BETA_CLOSURE_ANALYSIS  Documents the empirical closure beta = 1 + 0.1/n
% of the extended van Genuchten-Mualem (EVM) model with the German database
% of Hohenbrink et al. (2023): RetMeas.csv and CondMeas.csv (572 soils).
%
%   R = EVM_beta_closure_analysis              % CSV files next to this file
%   R = EVM_beta_closure_analysis(dataDir)     % CSV files in dataDir
%   R = EVM_beta_closure_analysis(dataDir, 20) % quick test on the first 20 soils
%   R = EVM_beta_closure_analysis(dataDir, 'figures')   % redraw Figs. S1 and S2 from the saved results
%
% Answers Reviewer #2 (comment on beta) and Reviewer #4 (comment 2):
%   A. Eq. (6) is monotonic for every beta > 0. beta = 1 and beta > 1 differ
%      in near-saturation smoothness: C(0+) > 0 for beta = 1, C(0+) = 0 for
%      beta > 1.
%   B. Closure scan: every soil refitted for each c in beta = 1 + c/n with
%      theta_s > theta_r >= 0 (RMSE_theta, RMSE_K and capacity positivity).
%      The scan with the previous single-start procedure (independent
%      bounds) is Script S8 (EVM_beta_closure_recheck.m).
%   C. Soil-specific optimized beta (5-parameter fits) and a test of the
%      linear relationship between beta and 1/n.
%
% Fitting (lsqcurvefit, Eq. 22, unweighted theta residuals): theta_s > theta_r >= 0
%   enforced through theta_r = f*theta_s, 0 <= f <= 0.9999 (P2). Internally, a
%   search with independent bounds aligns the constrained optimum (Script S2);
%   its results are not reported.
%
% Suction h is a positive magnitude (cm): h = 10^pF, h0 = 1e7 cm, ha = 1/alpha.
% K: Kc from Eq. (14) (already contains Ks,c; Ks,c = measured KSAT),
%    Knc from Eq. (20) with l = 0, gamma = 0.75, eta = 2 and
%    Ks,nc = 86.2*theta(1e5 cm)*ha^-1.5 (Eq. 17), Kv from Appendix A (20 C).
%    RMSE_K is in log10 units (Eq. 24), computed on the unsaturated K data.
%
% Outputs (folder beta_closure_results next to the data):
%   EVM_beta_closure_report.txt  all numbers quoted in the manuscript/reply
%   EVM_closure_summary.csv      per c: RMSE_theta, RMSE_K, positivity counts
%   EVM_closure_scan.csv         per soil and c: parameters and diagnostics
%   EVM_soil_specific_beta.csv   per soil: free beta*, closure beta (Table S1 for all soils)
%   FigS1_Snc_monotonic_smoothness.png, FigS2_soil_specific_beta.png, EVM_beta_closure_results.mat
% Runtime: about 6 min with 24 parallel workers (R2025b, Optimization and
% Parallel Computing Toolboxes; parfor runs serially without a pool).

if nargin < 1 || isempty(dataDir)
    dataDir = fileparts(mfilename('fullpath'));
end
outDir = fullfile(dataDir, 'beta_closure_results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
tStart = tic;
if nargin > 1 && (ischar(nSoilsMax) || isstring(nSoilsMax)) && any(strcmpi(nSoilsMax, 'figures'))
    L = load(fullfile(outDir, 'EVM_beta_closure_results.mat'), 'R');   % redraw the figures only
    S = L.R.settings; TL = L.R.scanTable; SB = L.R.soilSpecific; ids = string(SB.Sample_ID);
    pr = string(TL.Procedure); sid = string(TL.Sample_ID); scan = cell(numel(ids), 1);
    for i = 1:numel(ids)                  % fitted parameters by procedure and c, from the scan table
        for P = "P2"
            q = TL(sid == ids(i) & pr == P, :); [~, o] = sort(q.c); q = q(o, :);
            scan{i}.(P).p = [q.alpha q.n q.theta_r q.theta_s]; scan{i}.(P).beta = q.beta;
        end
    end
    make_figures(outDir, S, L.R.summary, SB, scan, [], ids, find(abs(S.cGrid - S.cRef) < 1e-12), find(S.cGrid == 0));
    R = L.R; fprintf('Figures redrawn in %s.\n', outDir);
    return
end

%% ---------------- settings ----------------
S.h0      = 1e7;      % oven-dry suction (cm), Schneider and Goss (2012)
S.hmin    = 1e-2;     % wettest suction in the database (pF = -2)
S.nGrid   = 6000;     % log-spaced points of the capacity check
S.cGrid   = [0 0.01 0.02 0.05 0.08 0.09 0.10 0.11 0.12 0.15 0.20 0.30 0.50 1 2 5];
S.cRef    = 0.10;     % coefficient used in the manuscript
S.betaMax = 10;       % upper bound of the soil-specific beta
S.hg      = logspace(log10(S.hmin), log10(S.h0) - 1e-9, S.nGrid)';
S.hgWet   = logspace(-8, log10(S.hmin), 1000)';   % extra check of strict positivity below hmin
S.hg1cm   = S.hg(S.hg >= 1);                       % alternative range h >= 1 cm (pF >= 0)
opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxFunctionEvaluations', 4000, 'MaxIterations', 1000, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-10);
nc   = numel(S.cGrid);
kRef = find(abs(S.cGrid - S.cRef) < 1e-12);
k0   = find(S.cGrid == 0);

%% ---------------- data ----------------
ret  = readtable(fullfile(dataDir, 'RetMeas.csv'),  'TextType', 'string');
cond = readtable(fullfile(dataDir, 'CondMeas.csv'), 'TextType', 'string');
ret  = ret(isfinite(ret.pF) & isfinite(ret.theta), :);
ids  = unique(ret.Sample_ID, 'stable');
if nargin > 1 && ~isempty(nSoilsMax), ids = ids(1:min(nSoilsMax, numel(ids))); end
nS   = numel(ids);
[~, ir] = ismember(ret.Sample_ID, ids);
[~, ic] = ismember(cond.Sample_ID, ids);
D = struct('id', cellstr(ids), 'h', [], 'th', [], 'hK', [], 'K', [], 'Ksat', NaN);
for i = 1:nS
    D(i).h  = 10.^ret.pF(ir == i);
    D(i).th = ret.theta(ir == i);
    isKs  = ic == i & cond.MeasType == "KSAT";
    unsat = ic == i & cond.MeasType ~= "KSAT" & cond.k > 0;
    D(i).hK = 10.^cond.pF(unsat);
    D(i).K  = cond.k(unsat);
    if any(isKs), D(i).Ksat = cond.k(find(isKs, 1)); end
end
grp = extractBefore(ids, "_");
fprintf('Loaded %d soils, %d retention points, %d soils with measured KSAT.\n', ...
    nS, height(ret), sum(isfinite([D.Ksat])));

if isempty(gcp('nocreate'))          % parfor runs serially when no pool can start
    try
        pc = parcluster('Processes');             % in-session change only, profile not saved
        pc.NumWorkers = max(pc.NumWorkers, feature('numcores'));
        parpool(pc, feature('numcores'));
    catch
        try, parpool('Processes'); catch, end %#ok<NOCOM>
    end
end

%% ---------------- B. closure scan: beta = 1 + c/n ----------------
scan = cell(nS, 1);
parfor i = 1:nS
    scan{i} = scan_one_soil(D(i), S, opts);
end
fprintf('Closure scan finished (%.1f min).\n', toc(tStart)/60);

procs = {'P2'};                      % reported procedure (the independent-bounds search is internal)
SUM = struct();
for ip = 1:numel(procs)
    pr = procs{ip};
    RM  = cell2mat(cellfun(@(s) s.(pr).rmse', scan, 'UniformOutput', false));
    SSE = cell2mat(cellfun(@(s) s.(pr).sse',  scan, 'UniformOutput', false));
    NEG = cell2mat(cellfun(@(s) s.(pr).neg',  scan, 'UniformOutput', false));
    N1  = cell2mat(cellfun(@(s) s.(pr).neg1cm', scan, 'UniformOutput', false));
    STR = cell2mat(cellfun(@(s) s.(pr).strictNeg', scan, 'UniformOutput', false));
    RGT = cell2mat(cellfun(@(s) s.(pr).thrGtThs',  scan, 'UniformOutput', false));
    RK  = cell2mat(cellfun(@(s) s.(pr).rmseK', scan, 'UniformOutput', false));
    np  = cellfun(@(s) s.npts, scan);
    hasK = all(isfinite(RK), 2);
    T = table(S.cGrid', mean(RM)', median(RM)', sqrt(sum(SSE)'/sum(np)), ...
        sum(NEG)', sum(STR)', sum(RGT)', sum(NEG & RGT)', mean(RK(hasK, :))', median(RK(hasK, :))', ...
        'VariableNames', {'c', 'meanRMSE_theta', 'medianRMSE_theta', 'pooledRMSE_theta', ...
        'nNegC_measRange', 'nNegC_strict_h0plus', 'nThrGtThs', 'nNegC_with_ThrGtThs', ...
        'meanRMSE_logK', 'medianRMSE_logK'});
    T.nNegC_h1cm = sum(N1)';
    T.dMeanRMSE_theta_pct = 100*(T.meanRMSE_theta/T.meanRMSE_theta(kRef) - 1);
    T.dMeanRMSE_logK_pct  = 100*(T.meanRMSE_logK/T.meanRMSE_logK(kRef) - 1);
    T.Procedure = repmat(string(pr), nc, 1);
    SUM.(pr) = T; SUM.([pr '_nK']) = sum(hasK);
    SUM.([pr '_RM']) = RM; SUM.([pr '_NEG']) = NEG; SUM.([pr '_RK']) = RK;
end

%% ---------------- C. soil-specific beta (P2, beta free in [1, betaMax]) ----------------
ss = cell(nS, 1);
parfor i = 1:nS
    ss{i} = soil_specific(D(i), scan{i}, S, opts, kRef);
end
SB = struct2table(cell2mat(ss));
SB.Sample_ID = ids; SB.Group = grp;
SB.RMSE_closure = cellfun(@(s) s.P2.rmse(kRef), scan);
SB.RMSE_beta1   = cellfun(@(s) s.P2.rmse(k0), scan);
SB.beta_closure = 1 + S.cRef./cellfun(@(s) s.P2.p(kRef, 2), scan);
SB.gain_vs_closure_pct = 100*(SB.RMSE_closure - SB.RMSE_free)./SB.RMSE_closure;
SB.profile_range = max(SB.RMSE_beta1, SB.RMSE_closure) - SB.RMSE_free;
SB.beta_no_effect = SB.theta_r < 0.02;          % theta_r*Snc negligible: beta cannot affect theta(h)
SB.atLowerBound = SB.beta_opt <= 1.01;
SB.atUpperBound = SB.beta_opt >= 0.99*S.betaMax;
SB.identifiable = SB.theta_r >= 0.05 & ~SB.atLowerBound & ~SB.atUpperBound;
SB = movevars(SB, {'Sample_ID', 'Group'}, 'Before', 1);

x = 1./SB.n; y = SB.beta_opt - 1;
subs = {true(nS, 1), ~SB.beta_no_effect, SB.theta_r >= 0.05, SB.identifiable};
subName = {'all soils', 'theta_r >= 0.02 (beta affects theta)', 'theta_r >= 0.05', ...
    'identifiable (theta_r >= 0.05, 1.01 < beta* < 0.99*betaMax)'};
REG = struct([]);
for k = 1:numel(subs)
    sel = subs{k};
    REG(k).subset = subName{k}; REG(k).N = sum(sel);
    REG(k).invn = ols_line(x(sel), y(sel));          % beta - 1 = a + b/n
    REG(k).n    = ols_line(SB.n(sel), SB.beta_opt(sel)); % beta = a + b*n
end
fprintf('Soil-specific fits finished (%.1f min).\n', toc(tStart)/60);

%% ---------------- A. near-saturation smoothness numbers ----------------
pRef = cell2mat(cellfun(@(s) s.P2.p(kRef, :), scan, 'UniformOutput', false));
bRef = 1 + S.cRef./pRef(:, 2);
selA = pRef(:, 3) >= 0.05;
C0beta1 = pRef(:, 3).*pRef(:, 1)./log1p(pRef(:, 1)*S.h0);   % C(0+) of Eq. (9) for beta = 1 (cm^-1)
hStar   = (bRef - 1).^(1./bRef)./pRef(:, 1);                % suction of max non-capillary capacity (cm)

%% ---------------- outputs ----------------
% long table of the closure scan
rows = cell(nS*nc*numel(procs), 17); r = 0;
for i = 1:nS
    for ip = 1:numel(procs)
        s = scan{i}.(procs{ip});
        for k = 1:nc
            r = r + 1;
            rows(r, :) = {ids(i), grp(i), string(procs{ip}), S.cGrid(k), s.p(k, 1), s.p(k, 2), ...
                s.p(k, 3), s.p(k, 4), s.beta(k), s.rmse(k), s.sse(k), scan{i}.npts, ...
                s.neg(k), s.hNegMax(k), ~s.strictNeg(k), s.rmseK(k), numel(D(i).K)};
        end
    end
end
TL = cell2table(rows, 'VariableNames', {'Sample_ID', 'Group', 'Procedure', 'c', 'alpha', 'n', ...
    'theta_r', 'theta_s', 'beta', 'RMSE_theta', 'SSE', 'Npts', 'NegC_measRange', 'hNegMax_cm', ...
    'StrictPositive_h0plus', 'RMSE_logK', 'NK'});
writetable(TL, fullfile(outDir, 'EVM_closure_scan.csv'));
writetable(SUM.P2, fullfile(outDir, 'EVM_closure_summary.csv'));
writetable(SB, fullfile(outDir, 'EVM_soil_specific_beta.csv'));

fid = fopen(fullfile(outDir, 'EVM_beta_closure_report.txt'), 'w');
if fid < 0, error('Cannot write the report in %s', outDir); end
for f = [1, fid]
    print_report(f, S, SUM, REG, SB, ids, C0beta1, hStar, selA, nS, kRef);
end
fclose(fid);

make_figures(outDir, S, SUM, SB, scan, D, ids, kRef, k0);
R = struct('settings', S, 'summary', SUM, 'soilSpecific', SB, 'regression', REG, 'scanTable', TL);
save(fullfile(outDir, 'EVM_beta_closure_results.mat'), 'R');
fprintf('Outputs written to %s (%.1f min).\n', outDir, toc(tStart)/60);
fprintf('EVM_BETA_ANALYSIS_DONE\n');
end

%% =====================================================================
%% local functions
%% =====================================================================
function out = scan_one_soil(d, S, opts)
% Refits the four EVM parameters for every c of S.cGrid (beta = 1 + c/n) with
% theta_s > theta_r >= 0 (P2), multi-start plus continuation in c; the search with
% independent bounds (P1) only aligns the P2 optimum and is not reported.
nc = numel(S.cGrid);
kRef = find(abs(S.cGrid - S.cRef) < 1e-12);
out.npts = numel(d.th);
Q = struct();
for ip = 1:2
    pr = sprintf('P%d', ip);
    % forward sweep in c (multi-start + previous solution, 12 extra starts at
    % c = 0 and c = cRef), then a backward sweep, so that every c keeps the
    % lowest SSE found from its neighbours
    Q0 = start_set(pr, d.th); Qx = start_extra(pr, d.th);
    qAll = NaN(nc, 4); sAll = Inf(nc, 1); qPrev = [];
    for k = 1:nc
        Qk = [qPrev; Q0];
        if k == 1 || k == kRef, Qk = [Qk; Qx]; end %#ok<AGROW>
        [~, qAll(k, :), sAll(k)] = fit_evm(d.h, d.th, 'closure', S.cGrid(k), pr, Qk, S, opts);
        qPrev = qAll(k, :);
    end
    for k = nc-1:-1:1
        [~, q, sse] = fit_evm(d.h, d.th, 'closure', S.cGrid(k), pr, qAll(k+1, :), S, opts);
        if sse < sAll(k), qAll(k, :) = q; sAll(k) = sse; end
    end
    Q.(pr).q = qAll; Q.(pr).sse = sAll;
end
% the P2 set lies inside the P1 set: a lower P2 SSE is also a P1 solution, and
% a P1 optimum with theta_r < theta_s is also the P2 optimum
for k = 1:nc
    p1 = q2p(Q.P1.q(k, :), 'P1'); p2 = q2p(Q.P2.q(k, :), 'P2');
    if Q.P2.sse(k) < Q.P1.sse(k)*(1 - 1e-12)
        Q.P1.q(k, :) = [log10(p2(1)) p2(2) p2(3) p2(4)]; Q.P1.sse(k) = Q.P2.sse(k);
    elseif p1(3) <= 0.9999*p1(4) && p1(4) >= 0.01 && Q.P1.sse(k) < Q.P2.sse(k)*(1 - 1e-12)
        Q.P2.q(k, :) = [log10(p1(1)) p1(2) p1(4) p1(3)/p1(4)]; Q.P2.sse(k) = Q.P1.sse(k);
    end
end
for ip = 1:2
    pr = sprintf('P%d', ip);
    s.p = NaN(nc, 4); s.q = NaN(nc, 4); s.beta = NaN(nc, 1); s.rmse = NaN(nc, 1);
    s.sse = NaN(nc, 1); s.neg = false(nc, 1); s.hNegMax = NaN(nc, 1); s.neg1cm = false(nc, 1);
    s.strictNeg = false(nc, 1); s.thrGtThs = false(nc, 1); s.rmseK = NaN(nc, 1);
    for k = 1:nc
        c = S.cGrid(k); q = Q.(pr).q(k, :); sse = Q.(pr).sse(k); p = q2p(q, pr);
        beta = 1 + c/p(2);
        C = evm_capacity(S.hg, p, beta, S.h0);
        s.p(k, :) = p; s.q(k, :) = q; s.beta(k) = beta;
        s.sse(k) = sse; s.rmse(k) = sqrt(sse/numel(d.th));
        s.neg(k) = any(C < 0);
        if s.neg(k), s.hNegMax(k) = max(S.hg(C < 0)); end
        s.neg1cm(k) = any(C(S.hg >= 1) < 0);
        s.thrGtThs(k) = p(3) > p(4);
        % sign of C as h -> 0+: if theta_r > theta_s the capillary term is
        % negative and the non-capillary term dominates only when beta < n
        s.strictNeg(k) = s.neg(k) || any(evm_capacity(S.hgWet, p, beta, S.h0) < 0) || ...
            (p(3) > p(4) && beta >= p(2));
        if isfinite(d.Ksat) && ~isempty(d.K)
            Kp = evm_K(d.hK, p, beta, d.Ksat, S.h0);
            s.rmseK(k) = sqrt(mean((log10(Kp) - log10(d.K)).^2));
        end
    end
    out.(pr) = s;
end
end

function o = soil_specific(d, sc, S, opts, kRef)
% Soil-specific beta: alpha, n, theta_r, theta_s and beta fitted together
% (procedure P2, so C > 0 for every beta > 0 by the proof of Section 2.2).
qa = sc.P2.q(kRef, :); qb = sc.P2.q(end, :);
Q0 = [repmat(qa, 5, 1), [1.05; 1.5; 3; 6; 9]; qb, min(sc.P2.beta(end), S.betaMax - 0.01)];
[p, q, sse] = fit_evm(d.h, d.th, 'free', NaN, 'P2', Q0, S, opts);
o.alpha = p(1); o.n = p(2); o.theta_r = p(3); o.theta_s = p(4);
o.beta_opt = q(5); o.RMSE_free = sqrt(sse/numel(d.th));
end

function [p, qBest, sseBest] = fit_evm(h, th, rule, c, pr, Q0, S, opts)
% Least-squares fit of Eq. (7) to theta(h) (Eq. 22) with lsqcurvefit.
% q = [log10(alpha) n theta_r theta_s] (P1) or [log10(alpha) n theta_s f] (P2),
% plus beta as q(5) when rule = 'free'.
if strcmp(pr, 'P1'), lb = [-5 1.0001 0 0];    ub = [0 10 1 1];
else,                lb = [-5 1.0001 0.01 0]; ub = [0 10 1 0.9999];
end
if strcmp(rule, 'free'), lb = [lb 1]; ub = [ub S.betaMax]; end
fun = @(q, hh) evm_theta(hh, q2p(q, pr), beta_of(q, rule, c), S.h0);
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

function b = beta_of(q, rule, c)
switch rule
    case 'closure', b = 1 + c/q(2);
    case 'fixed',   b = c;
    otherwise,      b = q(5);
end
end

function th = evm_theta(h, p, beta, h0)
% Eq. (7): theta = (ths - thr) Sc + thr Snc, with Eqs. (3), (4) and (6)
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n;
G   = (1 + (al*h).^n).^(-m);
G0  = (1 + (al*h0)^n)^(-m);
Sc  = (G - G0)/(1 - G0);
Snc = 1 - log1p((al*h).^beta)/log1p((al*h0)^beta);
th  = (ths - thr)*Sc + thr*Snc;
end

function C = evm_capacity(h, p, beta, h0)
% Eq. (9) with h > 0: C = (ths - thr)(-dSc/dh) + thr(-dSnc/dh), Eqs. (10)-(11)
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n; x = al*h;
G0   = (1 + (al*h0)^n)^(-m);
dSc  = m*n*al*x.^(n - 1)./(1 + x.^n).^(m + 1)/(1 - G0);
dSnc = beta*al*x.^(beta - 1)./((1 + x.^beta)*log1p((al*h0)^beta));
C    = (ths - thr)*dSc + thr*dSnc;
end

function K = evm_K(h, p, beta, Ksc, h0)
% K = Kc + Knc + Kv (Eq. 21 with Ks,c counted once, because Eq. 14 is dimensional)
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n;
G   = (1 + (al*h).^n).^(-m);
G0  = (1 + (al*h0)^n)^(-m);
Sc  = max((G - G0)/(1 - G0), 0);
Snc = max(1 - log1p((al*h).^beta)/log1p((al*h0)^beta), 0);
Kc  = Ksc*sqrt(Sc).*(1 - ((1 - G.^(1/m))/(1 - G0^(1/m))).^m).^2;   % Eq. (14)
thm  = evm_theta(1e5, p, beta, h0);                                 % theta at h = 1e5 cm
Ksnc = 86.2*thm*al^1.5;                                             % Eq. (17), ha = 1/alpha
Knc  = Ksnc*(al*h0).^(-1.5*(1 - Snc));                              % Eq. (20), l = 0
T = 293.15; Mw = 0.018015; g = 9.81; Rg = 8.314;                    % Appendix A, 20 C
rho_sv = 1e-3*exp(31.3716 - 6014.79/T - 7.92495e-3*T)/T;            % kg m-3
rho_w  = 998.2;                                                     % kg m-3
Da  = 2.12e-5*(T/273.15)^2*1e4*86400;                               % cm2 d-1
cH  = Mw*g/(Rg*T)/100;                                              % cm-1
th  = (ths - thr)*Sc + thr*Snc;
tau = max(ths - th, 0).^(10/3)/ths^2;                               % Eq. (A4)
Kv  = rho_sv/rho_w*tau*Da*cH.*exp(-cH*h);                           % Eqs. (A1)-(A3)
K   = max(Kc + Knc + Kv, realmin);
end

function st = ols_line(x, y)
% Ordinary least squares y = a + b x, plus the through-origin fit y = b0 x.
x = x(:); y = y(:); N = numel(x);
st = struct('a', NaN, 'b', NaN, 'ciA', NaN, 'ciB', NaN, 'pB', NaN, 'R2', NaN, 'b0', NaN, 'R2_0', NaN);
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

function print_report(f, S, SUM, REG, SB, ids, C0beta1, hStar, selA, nS, kRef)
fprintf(f, '\n================ EVM closure beta = 1 + c/n : German database (%d soils) ================\n', nS);
fprintf(f, 'Closure scan: alpha, n, theta_r, theta_s refitted for each c (Eq. 22), theta_s > theta_r >= 0 enforced\n');
fprintf(f, 'Capacity checked for C(h) > 0 on %.0e <= h < %.0e cm (and as h -> 0+)\n', S.hmin, S.h0);
for pr = {'P2'}
    T = SUM.(pr{1});
    lab = 'theta_s > theta_r >= 0 enforced (P2)';
    fprintf(f, '\n--- %s; RMSE_K over %d soils with measured KSAT ---\n', lab, SUM.([pr{1} '_nK']));
    fprintf(f, '%6s %10s %10s %10s %8s %8s %8s %8s %8s %9s %8s %9s\n', 'c', 'meanRMSEt', 'dRMSEt%', 'poolRMSEt', ...
        'nNegC', 'nStrict', 'nNeg1cm', 'thr>ths', 'neg&>', 'meanRMSEK', 'dRMSEK%', 'medRMSEK');
    for k = 1:height(T)
        fprintf(f, '%6.2f %10.6f %10.3f %10.6f %8d %8d %8d %8d %8d %9.4f %8.3f %9.4f\n', T.c(k), T.meanRMSE_theta(k), ...
            T.dMeanRMSE_theta_pct(k), T.pooledRMSE_theta(k), T.nNegC_measRange(k), T.nNegC_strict_h0plus(k), ...
            T.nNegC_h1cm(k), T.nThrGtThs(k), T.nNegC_with_ThrGtThs(k), T.meanRMSE_logK(k), T.dMeanRMSE_logK_pct(k), ...
            T.medianRMSE_logK(k));
    end
end
fprintf(f, '\nSoil-specific beta (P2, beta free in [1, %g]):\n', S.betaMax);
e = ~SB.beta_no_effect;
fprintf(f, '   theta_r < 0.02 (beta has no effect on theta): %d soils; remaining %d soils:\n', sum(~e), sum(e));
fprintf(f, '   beta* quantiles (0, 10, 25, 50, 75, 90, 100%%): %s\n', mat2str(quantile(SB.beta_opt(e), [0 .1 .25 .5 .75 .9 1]), 4));
fprintf(f, '   at lower bound (<= 1.01): %d | at upper bound: %d | interior: %d | identifiable (theta_r >= 0.05, interior): %d\n', ...
    sum(SB.atLowerBound & e), sum(SB.atUpperBound & e), sum(e & ~SB.atLowerBound & ~SB.atUpperBound), sum(SB.identifiable));
fprintf(f, '   RMSE gain of free beta over the closure: median %.2e (%.2f%%), 90th pct %.2f%%, max %.2f%%\n', ...
    median(SB.RMSE_closure - SB.RMSE_free), median(SB.gain_vs_closure_pct), ...
    quantile(SB.gain_vs_closure_pct, 0.9), max(SB.gain_vs_closure_pct));
fprintf(f, '   soils whose RMSE changes < 1e-4 between beta = 1, closure and beta*: %d of %d\n', sum(SB.profile_range < 1e-4), nS);
g = SB.gain_vs_closure_pct(e); ga = SB.RMSE_closure(e) - SB.RMSE_free(e);
fprintf(f, '   theta_r >= 0.02: RMSE gain median %.2f%% (%.2e cm3 cm-3), IQR %.2f-%.2f%%, > 5%%: %d, > 10%%: %d, RMSE range < 1e-4: %d\n', ...
    median(g), median(ga), quantile(g, 0.25), quantile(g, 0.75), sum(g > 5), sum(g > 10), sum(SB.profile_range(e) < 1e-4));
dn = SB.beta_opt(SB.identifiable) - SB.n(SB.identifiable);
fprintf(f, '   identifiable soils: beta* - n median %.2f (IQR %.2f-%.2f)\n', median(dn), quantile(dn, 0.25), quantile(dn, 0.75));
for k = 1:numel(REG)
    a = REG(k).invn; b = REG(k).n;
    fprintf(f, '   [%s] N = %d: beta-1 = %.3f(+-%.3f) + %.3f(+-%.3f)/n, R2 = %.3f, p = %.2g; through origin b0 = %.3f, R2 = %.3f | beta = %.3f + %.3f n, R2 = %.3f\n', ...
        REG(k).subset, REG(k).N, a.a, a.ciA, a.b, a.ciB, a.R2, a.pB, a.b0, a.R2_0, b.a, b.b, b.R2);
end
t1 = ["FSUJ_023" "FSUJ_035" "TUBS_027" "TUBS_128" "KIT_022" "KIT_009" "TUBS_106"];
fprintf(f, '\nTable-1 German soils (P2): closure beta, soil-specific beta*, RMSE closure / beta*\n');
for s = t1
    i = find(ids == s, 1);
    if isempty(i), continue, end
    fprintf(f, '   %-9s n = %.3f theta_r = %.3f  beta_closure = %.4f  beta* = %.3f  RMSE %.5f / %.5f\n', s, ...
        SB.n(i), SB.theta_r(i), SB.beta_closure(i), SB.beta_opt(i), SB.RMSE_closure(i), SB.RMSE_free(i));
end
fprintf(f, '\nNear-saturation smoothness (P2 fits, c = %.2f, soils with theta_r >= 0.05, N = %d):\n', S.cRef, sum(selA));
fprintf(f, '   beta = 1      : C(0+) = theta_r*alpha/ln(1+alpha*h0), median %.2e cm^-1 (range %.2e - %.2e)\n', ...
    median(C0beta1(selA)), min(C0beta1(selA)), max(C0beta1(selA)));
fprintf(f, '   beta = 1+0.1/n: C(0+) = 0; non-capillary capacity peaks at h* = (beta-1)^(1/beta)/alpha, median %.2f cm (range %.3f - %.1f)\n', ...
    median(hStar(selA)), min(hStar(selA)), max(hStar(selA)));
fprintf(f, '   c = %.2f summary: mean RMSE_theta = %.6f (theta_s > theta_r >= 0)\n', S.cGrid(kRef), SUM.P2.meanRMSE_theta(kRef));
end

function make_figures(outDir, S, SUM, SB, scan, D, ids, kRef, k0)
% Figure S1: Eq. (6) is monotonic for every beta; beta controls C near saturation
f1 = figure('Visible', 'off', 'Position', [100 100 1100 420]);
x = logspace(-6, 5, 2000)'; x0 = 1e5; bet = [0.8 1 1+0.1/2 1+0.1/1.2 1.5 2];
lab = {'\beta = 0.8', '\beta = 1', '\beta = 1+0.1/n, n = 2', '\beta = 1+0.1/n, n = 1.2', '\beta = 1.5', '\beta = 2'};
subplot(1, 2, 1); hold on
for b = bet, semilogx(x, 1 - log1p(x.^b)/log1p(x0^b), 'LineWidth', 1.4); end
set(gca, 'XScale', 'log'); xlabel('h/h_a'); ylabel('S_{nc} (Eq. 6)'); box on; grid on
title('(a) S_{nc} decreases monotonically for every \beta > 0'); legend(lab, 'Location', 'southwest');
subplot(1, 2, 2); hold on
for b = bet, loglog(x, b*x.^(b - 1)./((1 + x.^b)*log1p(x0^b)), 'LineWidth', 1.4); end
set(gca, 'XScale', 'log', 'YScale', 'log'); xlim([1e-6 10]); box on; grid on
xlabel('h/h_a'); ylabel('-dS_{nc}/d(h/h_a)');
title('(b) h\rightarrow0: \infty (\beta<1), finite (\beta=1), 0 (\beta>1)');
exportgraphics(f1, fullfile(outDir, 'FigS1_Snc_monotonic_smoothness.png'), 'Resolution', 300); close(f1);

% Figure S2: soil-specific beta
f3 = figure('Visible', 'off', 'Position', [100 100 1100 420]);
subplot(1, 2, 1);
e = ~SB.beta_no_effect;
scatter(1./SB.n(~e), SB.beta_opt(~e), 10, [0.75 0.75 0.75], 'filled'); hold on
scatter(1./SB.n(e), SB.beta_opt(e), 16, SB.theta_r(e), 'filled');
xi = linspace(0.1, 1, 100); plot(xi, 1 + 0.1*xi, 'r-', 'LineWidth', 2);
set(gca, 'YScale', 'log', 'YTick', [1 1.1 1.2 1.5 2 3 5 10]); ylim([0.98 10.5]);
cb = colorbar; cb.Label.String = '\theta_r (cm^3 cm^{-3})'; grid on; box on
xlabel('1/n'); ylabel('soil-specific \beta^*');
legend('\theta_r < 0.02 (\beta has no effect)', '\beta^*, \theta_r \geq 0.02', '\beta = 1 + 0.1/n', 'Location', 'northwest');
title('(a) \beta^* from retention data vs 1/n');
subplot(1, 2, 2);
edges = [-Inf 0.01 0.1 1 5 10 Inf];
cnt = histcounts(SB.gain_vs_closure_pct(e), edges);
bar(cnt); grid on
set(gca, 'XTickLabel', {'<0.01', '0.01-0.1', '0.1-1', '1-5', '5-10', '>10'});
xlabel('RMSE_\theta reduction by a free \beta relative to the closure (%)'); ylabel('soils (\theta_r \geq 0.02)');
title(sprintf('(b) median reduction %.2f%% (all soils)', median(SB.gain_vs_closure_pct)));
exportgraphics(f3, fullfile(outDir, 'FigS2_soil_specific_beta.png'), 'Resolution', 300); close(f3);

end
