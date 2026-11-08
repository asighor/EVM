function R = EVM_positivity_constraint_analysis(dataDir, nSoilsMax)
%EVM_POSITIVITY_CONSTRAINT_ANALYSIS  Aligns the positivity proof of the EVM soil
% water capacity (Section 2.2) with the fitting procedure (Eq. 22), using the
% German database of Hohenbrink et al. (2023): RetMeas.csv and CondMeas.csv
% (572 soils). The closure beta = 1 + 0.1/n is used throughout.
%
%   R = EVM_positivity_constraint_analysis              % CSV files next to this file
%   R = EVM_positivity_constraint_analysis(dataDir)     % CSV files in dataDir
%   R = EVM_positivity_constraint_analysis(dataDir, 20) % quick test on the first 20 soils
%   R = EVM_positivity_constraint_analysis(dataDir, ["TUBS_205" "FSUJ_015"])  % selected soils
%
% Text S3 of the Supporting Information (positivity, endpoints and order of continuity):
%   A. Proposition: for alpha > 0, n > 1, beta > 0, 0 < h0 < inf and
%      theta_s > theta_r >= 0, C(h) > 0 for every 0 < h < h0 (Eq. 9).
%   B. Fitting (lsqcurvefit, Eq. 22, unweighted theta residuals), theta_s > theta_r >= 0 enforced
%      through theta_r = f*theta_s, 0 <= f <= fMax < 1 (P2), from 30 starting points (110 where
%      the constraint is active). Internally, an independent-bounds search (P1) aligns the
%      constrained optimum, and where the constraint is active the optimum is checked with the
%      constraint written as theta_r - fMax*theta_s <= 0 (linear inequality, interior-point; P2x).
%      Only the constrained fits are reported.
%   C. Strict positivity on 0 < h < h0 (20,000-point grid from 1e-12 cm to h0, plus
%      the sign of the leading term as h -> 0+), separated from the endpoints:
%      theta(0) = theta_s and C(0+) = 0 (n > 1, beta > 1); theta(h0) = 0 and C(h0) > 0.
%   D. Order of continuity: C is continuous at h = 0, so theta(h) is C1 on [0, h0);
%      dC/dh ~ h^(min(n,beta)-2) is unbounded as h -> 0+, so theta is not C2 there.
%   E. RMSE_theta and RMSE_K of the constrained fits (Table 3 by soil group).
%
% Suction h is a positive magnitude (cm): h = 10^pF, h0 = 1e7 cm, ha = 1/alpha.
% K = Kc + Knc + Kv as in EVM_beta_closure_analysis.m: Kc from Eq. (14) with
% Ks,c = measured KSAT (counted once), Knc from Eq. (20) with l = 0, gamma = 0.75,
% eta = 2, Kv from Appendix A (20 C). RMSE_K is in log10 units (Eq. 24).
%
% Outputs (folder positivity_results next to the data):
%   EVM_positivity_report.txt   all numbers quoted in the manuscript and the reply
%   EVM_positivity_fits.csv     per soil: parameters, errors and positivity, endpoint and
%                               continuity diagnostics of the constrained fit (P2_ columns)
%   EVM_positivity_groups.csv   mean errors by soil group
%   EVM_TUBS205_curves.csv      theta(h) and C(h) of the Figure 3 soil
%   FigS6_endpoints_continuity.png, FigS7_capacity_all_soils.png, EVM_positivity_results.mat
% Runtime: a few minutes with 24 parallel workers (R2025b, Optimization and
% Parallel Computing Toolboxes; parfor runs serially without a pool).

if nargin < 1 || isempty(dataDir)
    dataDir = fileparts(mfilename('fullpath'));
end
outDir = fullfile(dataDir, 'positivity_results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
tStart = tic;

%% ---------------- settings ----------------
S.h0    = 1e7;        % oven-dry suction (cm), Schneider and Goss (2012)
S.c     = 0.1;        % closure beta = 1 + c/n
S.hmin  = 1e-2;       % wettest suction in the database (pF = -2)
S.hLow  = 1e-12;      % wet end of the strict-positivity grid (cm)
S.nGrid = 20000;      % log-spaced points of the positivity check on 0 < h < h0
S.fMax  = 0.9999;     % theta_r <= fMax*theta_s, i.e. theta_s > theta_r
S.nMin  = 1.0001;     % lower bound of n (n > 1)
S.hg    = logspace(log10(S.hLow), log10(S.h0) - 1e-9, S.nGrid)';
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

%% ---------------- B. fits with theta_s > theta_r >= 0 (P2; P1 internal) ----------------
fits = cell(nS, 1);
parfor i = 1:nS
    fits{i} = fit_one_soil(D(i), S, opts, optsLC);
end
fprintf('Fits finished (%.1f min).\n', toc(tStart)/60);

T1 = struct2table(cell2mat(cellfun(@(o) o.P1, fits, 'UniformOutput', false)));
T2 = struct2table(cell2mat(cellfun(@(o) o.P2, fits, 'UniformOutput', false)));
T1.Properties.VariableNames = strcat('P1_', T1.Properties.VariableNames);
T2.Properties.VariableNames = strcat('P2_', T2.Properties.VariableNames);
FT = table(ids, grp, arrayfun(@(d) numel(d.th), D), arrayfun(@(d) numel(d.K), D), [D.Ksat]', ...
    'VariableNames', {'Sample_ID', 'Group', 'Npts', 'NK', 'Ksat_cm_d'});
FT = [FT, T1, T2];
FT.active      = cellfun(@(o) o.active, fits);        % P1 optimum violates theta_s > theta_r
FT.fixP1       = cellfun(@(o) o.fixP1, fits);
FT.fixP2       = cellfun(@(o) o.fixP2, fits);
FT.fixP2x      = cellfun(@(o) o.fixP2x, fits);
FT.SSE_P2_reparam = cellfun(@(o) o.sseP2_raw, fits);
FT.SSE_P2x     = cellfun(@(o) o.sseP2x, fits);
FT.relDiff_P2x = (FT.SSE_P2x - FT.SSE_P2_reparam)./FT.SSE_P2_reparam;
FT.relDiff_P2_P1_feasible = (FT.SSE_P2_reparam - cellfun(@(o) o.sseP1_raw, fits))./FT.SSE_P2_reparam;
FT.relDiff_P2_P1_feasible(FT.active) = NaN;
FT.tFit_s      = cellfun(@(o) o.tFit, fits);
FT.dRMSE_theta = FT.P2_RMSE_theta - FT.P1_RMSE_theta;
FT.dRMSE_logK  = FT.P2_RMSE_logK - FT.P1_RMSE_logK;

% Table 3 by soil group
gl = ["TUDD" "UFZ" "FSUJ" "KIT" "TUBS" "All"];
rows = cell(2*numel(gl), 10); r = 0;
for pr = ["P1" "P2"]
    for g = gl
        sel = grp == g | g == "All";
        k = sel & isfinite(FT.(pr + "_RMSE_logK"));
        r = r + 1;
        rows(r, :) = {pr, g, sum(sel), sum(k), mean(FT.(pr + "_RMSE_theta")(sel)), ...
            mean(FT.(pr + "_RMSE_logK")(k)), median(FT.(pr + "_RMSE_logK")(k)), ...
            sum(FT.(pr + "_thrGtThs")(sel)), sum(~FT.(pr + "_strictPositive")(sel)), sum(FT.(pr + "_negMeasRange")(sel))};
    end
end
G = cell2table(rows, 'VariableNames', {'Procedure', 'Group', 'N', 'NK', 'meanRMSE_theta', ...
    'meanRMSE_logK', 'medianRMSE_logK', 'nThrGtThs', 'nNotStrictlyPositive', 'nNegC_measRange'});
G.Procedure = string(G.Procedure); G.Group = string(G.Group);

%% ---------------- outputs (constrained fits only) ----------------
vn = FT.Properties.VariableNames;
FTw = FT(:, ~startsWith(vn, 'P1_') & ~ismember(vn, {'active', 'fixP1', 'fixP2', 'relDiff_P2_P1_feasible', ...
    'dRMSE_theta', 'dRMSE_logK'}));
Gw = G(G.Procedure == "P2", :);
writetable(FTw, fullfile(outDir, 'EVM_positivity_fits.csv'));
writetable(Gw,  fullfile(outDir, 'EVM_positivity_groups.csv'));
i = find(ids == "TUBS_205", 1);
if ~isempty(i)
    hh = logspace(-4, log10(S.h0) - 1e-9, 500)';
    p2 = pvec(FT, 'P2', i);
    TC = table(hh, log10(hh), evm_theta(hh, p2, 1 + S.c/p2(2), S.h0), evm_capacity(hh, p2, 1 + S.c/p2(2), S.h0), ...
        'VariableNames', {'h_cm', 'pF', 'theta_P2', 'C_P2_per_cm'});
    writetable(TC, fullfile(outDir, 'EVM_TUBS205_curves.csv'));
end

fid = fopen(fullfile(outDir, 'EVM_positivity_report.txt'), 'w');
if fid < 0, error('Cannot write the report in %s', outDir); end
for f = [1, fid]
    print_report(f, S, FT, G, nS);
end
fclose(fid);

make_figures(outDir, S, FT, D, ids, grp);
R = struct('settings', S, 'fits', FTw, 'groups', Gw);
save(fullfile(outDir, 'EVM_positivity_results.mat'), 'R');
fprintf('Outputs written to %s (%.1f min).\n', outDir, toc(tStart)/60);
fprintf('EVM_POSITIVITY_ANALYSIS_DONE\n');
end

%% =====================================================================
%% local functions
%% =====================================================================
function o = fit_one_soil(d, S, opts, optsLC)
% P1 and P2 fits of one soil. The P2 set lies inside the P1 set, so a lower
% P2 SSE replaces P1, and a feasible P1 optimum is also the P2 optimum. Where
% the P1 optimum violates the constraint, P2x checks P2 from independent starts.
t0 = tic;
[p1, ~, s1] = fit_evm(d.h, d.th, 'P1', start_set('P1', d.th, false), S, opts, optsLC);
[p2, ~, s2] = fit_evm(d.h, d.th, 'P2', start_set('P2', d.th, false), S, opts, optsLC);
if p1(3) > S.fMax*p1(4)          % constraint active: denser multi-start for both procedures
    [pb, ~, sb] = fit_evm(d.h, d.th, 'P1', start_set('P1', d.th, true), S, opts, optsLC);
    if sb < s1, p1 = pb; s1 = sb; end
    [pb, ~, sb] = fit_evm(d.h, d.th, 'P2', start_set('P2', d.th, true), S, opts, optsLC);
    if sb < s2, p2 = pb; s2 = sb; end
end
o.sseP1_raw = s1; o.sseP2_raw = s2; o.sseP2x = NaN; o.fixP2x = false;
o.fixP1 = s2 < s1*(1 - 1e-9);
if o.fixP1, p1 = p2; s1 = s2; end
o.active = p1(3) > S.fMax*p1(4);
o.fixP2 = ~o.active && s1 < s2*(1 - 1e-9);
if o.fixP2, p2 = p1; s2 = s1; end
if o.active
    [px, ~, sx] = fit_evm(d.h, d.th, 'P2x', start_set('P2x', d.th, false), S, opts, optsLC);
    o.sseP2x = sx;
    o.fixP2x = sx < s2*(1 - 1e-9) && px(3) < px(4);
    if o.fixP2x, p2 = px; s2 = sx; end
end
o.P1 = diagnose(d, p1, s1, S);
o.P2 = diagnose(d, p2, s2, S);
o.tFit = toc(t0);
end

function d = diagnose(dd, p, sse, S)
% Positivity, endpoint and continuity diagnostics of one fitted parameter set.
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n; beta = 1 + S.c/n;
d.alpha = al; d.n = n; d.theta_r = thr; d.theta_s = ths; d.beta = beta; d.f = thr/ths;
d.SSE = sse; d.RMSE_theta = sqrt(sse/numel(dd.th));
d.thrGtThs = thr > ths;
C  = evm_capacity(S.hg, p, beta, S.h0);
th = evm_theta(S.hg, p, beta, S.h0);
neg = C <= 0;
d.nGridNeg = sum(neg);
d.hNegMin = NaN; d.hNegMax = NaN;
if any(neg), d.hNegMin = min(S.hg(neg)); d.hNegMax = max(S.hg(neg)); end
d.negMeasRange = any(neg & S.hg >= S.hmin);
% sign of C as h -> 0+: C ~ A (alpha h)^(n-1) + B (alpha h)^(beta-1); with
% theta_r > theta_s (A < 0) the term with the smaller exponent decides
G0 = (1 + (al*S.h0)^n)^(-m); L = log1p((al*S.h0)^beta);
A = (ths - thr)*m*n/(1 - G0); B = thr*beta/L;
if A >= 0,         s0 = sign(A + B);
elseif beta < n,   s0 = sign(B);
elseif beta > n,   s0 = sign(A);
else,              s0 = sign(A + B);
end
d.sign0 = s0;
d.strictPositive = ~any(neg) && s0 > 0;
d.thetaExceed = max(max(th) - ths, 0);            % water content above theta_s
if d.thetaExceed < 1e-12, d.thetaExceed = 0; end  % round-off
if thr > 0, d.e0 = min(n, beta) - 1; else, d.e0 = n - 1; end   % C ~ (alpha h)^e0 as h -> 0+
d.C_hLow_rel = C(1)/max(C);                       % C(1e-12 cm)/max C: how close C is to 0
d.Ch0 = evm_capacity(S.h0, p, beta, S.h0);        % C at the dry endpoint (cm^-1)
d.h0Ch0 = S.h0*d.Ch0;                             % -dtheta/d ln h at h0
d.C0beta1 = thr*al/log1p(al*S.h0);                % C(0+) if beta were 1 (cm^-1)
d.theta0_minus_ths = evm_theta(0, p, beta, S.h0) - ths;
d.theta_h0 = evm_theta(S.h0, p, beta, S.h0);
d.RMSE_logK = NaN;
if isfinite(dd.Ksat) && ~isempty(dd.K)
    Kp = evm_K(dd.hK, p, beta, dd.Ksat, S.h0);
    d.RMSE_logK = sqrt(mean((log10(Kp) - log10(dd.K)).^2));
end
d.nAtLower = n <= S.nMin + 1e-6;   d.nAtUpper = n >= 10 - 1e-6;
d.thrAtZero = thr <= 1e-6;         d.fAtUpper = d.f >= S.fMax - 1e-6;
d.thrAtOne = thr >= 1 - 1e-6;      d.thsAtOne = ths >= 1 - 1e-6;
d.alphaAtBound = al <= 1e-5*(1 + 1e-6) || al >= 1 - 1e-6;
end

function [p, qBest, sseBest] = fit_evm(h, th, pr, Q0, S, opts, optsLC)
% Least-squares fit of Eq. (7) to theta(h) (Eq. 22) with lsqcurvefit.
% q = [log10(alpha) n theta_r theta_s] (P1, P2x) or [log10(alpha) n theta_s f] (P2).
switch pr
    case 'P1',  lb = [-5 S.nMin 0 0];    ub = [0 10 1 1];
    case 'P2',  lb = [-5 S.nMin 0.01 0]; ub = [0 10 1 S.fMax];
    case 'P2x', lb = [-5 S.nMin 0 0.01]; ub = [0 10 1 1];
end
fun = @(q, hh) evm_theta(hh, q2p(q, pr), 1 + S.c/q(2), S.h0);
ws = warning;                    % interior-point KKT solves warn near the bounds of n
warning('off', 'MATLAB:nearlySingularMatrix'); warning('off', 'MATLAB:singularMatrix');
sseBest = Inf; qBest = Q0(1, :);
for j = 1:size(Q0, 1)
    q0 = min(max(Q0(j, :), lb + 1e-9), ub - 1e-9);
    try
        if strcmp(pr, 'P2x')     % theta_r - fMax*theta_s <= 0
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
% the 18 starts of EVM_beta_closure_analysis.m plus 12 starts that include
% n -> 1, large n and theta_r (or f) -> 1; dense adds an 80-point grid
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

function p = pvec(FT, pr, i)
p = [FT.([pr '_alpha'])(i) FT.([pr '_n'])(i) FT.([pr '_theta_r'])(i) FT.([pr '_theta_s'])(i)];
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

function dC = evm_dCdh(h, p, beta, h0)
% dC/dh = -d2theta/dh2 of Eq. (7), h > 0
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n; x = al*h;
G0 = (1 + (al*h0)^n)^(-m);
A  = (ths - thr)*m*n/(1 - G0); B = thr*beta/log1p((al*h0)^beta);
dC = al^2*(A*x.^(n - 2).*(1 + x.^n).^(-(m + 2)).*((n - 1)*(1 + x.^n) - (m + 1)*n*x.^n) + ...
           B*x.^(beta - 2)./(1 + x.^beta).^2.*((beta - 1)*(1 + x.^beta) - beta*x.^beta));
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

function print_report(f, S, FT, G, nS)
fprintf(f, '\n========== EVM positivity, endpoints and order of continuity: German database (%d soils), beta = 1 + %.1f/n ==========\n', nS, S.c);
fprintf(f, 'Proposition (Section 2.2): alpha > 0, n > 1, beta > 0, h0 < inf, theta_s > theta_r >= 0  =>  C(h) > 0 for 0 < h < h0\n');
fprintf(f, 'Fits (P2)  : log10(alpha) in [-5, 0], n in [%.4f, 10], theta_s in [0.01, 1], theta_r = f*theta_s, f in [0, %.4f]\n', S.nMin, S.fMax);
fprintf(f, 'Check (P2x): theta_r - %.4f*theta_s <= 0 as a linear inequality (interior-point) where the constraint is active\n', S.fMax);
fprintf(f, 'Positivity : C on %d log-spaced suctions %.0e cm <= h < h0 = %.0e cm, plus the sign of the leading term as h -> 0+\n', ...
    S.nGrid, S.hLow, S.h0);
rm = FT.P2_RMSE_theta; rk = FT.P2_RMSE_logK; hk = isfinite(rk);
fprintf(f, '\n--- Constrained fits ---\n');
fprintf(f, '   RMSE_theta: mean %.6f, median %.6f, pooled %.6f cm3 cm-3 | RMSE_K (%d soils with KSAT): mean %.4f, median %.4f log10 units\n', ...
    mean(rm), median(rm), sqrt(sum(FT.P2_SSE)/sum(FT.Npts)), sum(hk), mean(rk(hk)), median(rk(hk)));
fprintf(f, '   theta_r > theta_s: %d | not strictly positive on 0 < h < h0: %d (C < 0 within h >= %.0e cm: %d) | theta > theta_s somewhere: %d\n', ...
    sum(FT.P2_thrGtThs), sum(~FT.P2_strictPositive), S.hmin, sum(FT.P2_negMeasRange), sum(FT.P2_thetaExceed > 0));
fprintf(f, '   bounds hit: n = %.4f: %d | n = 10: %d | theta_r = 0: %d | f = fMax: %d | theta_s = 1: %d | alpha bound: %d\n', ...
    S.nMin, sum(FT.P2_nAtLower), sum(FT.P2_nAtUpper), sum(FT.P2_thrAtZero), sum(FT.P2_fAtUpper), sum(FT.P2_thsAtOne), ...
    sum(FT.P2_alphaAtBound));

fprintf(f, '\nTable 3 (EVM) by group (RMSE_theta cm3 cm-3 / RMSE_K log10 units):\n');
for g = ["TUDD" "UFZ" "FSUJ" "KIT" "TUBS" "All"]
    b = G(G.Procedure == "P2" & G.Group == g, :);
    fprintf(f, '   %-5s N = %3d (NK = %3d): %.4f / %.3f\n', g, b.N, b.NK, b.meanRMSE_theta, b.meanRMSE_logK);
end
t1 = ["FSUJ_023" "FSUJ_035" "TUBS_027" "TUBS_128" "KIT_022" "KIT_009" "TUBS_106"];
fprintf(f, '\nTable-2 German soils:\n');
for s = t1
    i = find(FT.Sample_ID == s, 1);
    if isempty(i), continue, end
    fprintf(f, '   %-9s thr %.3f ths %.3f alpha %.3f n %.3f RMSE %.4f\n', s, FT.P2_theta_r(i), FT.P2_theta_s(i), ...
        FT.P2_alpha(i), FT.P2_n(i), FT.P2_RMSE_theta(i));
end

e0 = FT.P2_e0; th0 = abs(FT.P2_theta0_minus_ths); thh0 = abs(FT.P2_theta_h0);
fprintf(f, '\nEndpoints and order of continuity:\n');
fprintf(f, '   theta(0) = theta_s: max |deviation| %.1e | theta(h0) = 0: max |theta(h0)| %.1e\n', max(th0), max(thh0));
fprintf(f, '   n > 1 and beta > 1 in all %d fits -> C(0+) = 0 analytically; C ~ (alpha h)^e0, e0 = min(n,beta) - 1 (theta_r > 0) or n - 1\n', ...
    sum(FT.P2_n > 1 & FT.P2_beta > 1));
fprintf(f, '   e0 quantiles (0, 5, 25, 50, 75, 95, 100%%): %s\n', mat2str(quantile(e0, [0 .05 .25 .5 .75 .95 1]), 3));
z = FT.P2_theta_r <= 1e-6;
fprintf(f, '   C not differentiable at h = 0 (e0 < 1, theta C1 but not C2): %d of %d | theta_r <= 1e-6: %d fits, largest n %.3f, n >= 2: %d\n', ...
    sum(e0 < 1), nS, sum(z), max(FT.P2_n(z)), sum(z & FT.P2_n >= 2));
fprintf(f, '   C(1e-12 cm)/max C > 0.5 (C reaches 0 only far below 1e-12 cm): %d soils (n < 1.01: %d)\n', ...
    sum(FT.P2_C_hLow_rel > 0.5), sum(FT.P2_n < 1.01));
fprintf(f, '   dry end: C(h0) > 0 in %d of %d fits; C(h0) median %.2e cm^-1 (range %.2e - %.2e); -dtheta/dpF at h0 = ln(10) h0 C(h0): median %.4f (range %.4f - %.4f)\n', ...
    sum(FT.P2_Ch0 > 0), nS, median(FT.P2_Ch0), min(FT.P2_Ch0), max(FT.P2_Ch0), median(log(10)*FT.P2_h0Ch0), ...
    min(log(10)*FT.P2_h0Ch0), max(log(10)*FT.P2_h0Ch0));
s5 = FT.P2_theta_r >= 0.05;
fprintf(f, '   for reference, beta = 1 would give C(0+) = theta_r alpha/ln(1+alpha h0): median %.2e cm^-1 (theta_r >= 0.05, N = %d)\n', ...
    median(FT.P2_C0beta1(s5)), sum(s5));

fprintf(f, '\nConsistency checks:\n');
k = isfinite(FT.relDiff_P2x);
if any(k)
    fprintf(f, '   linear-inequality check (%d soils with an active constraint): relative SSE difference min %+.1e, max %+.1e; optimum replaced in %d\n', ...
        sum(k), min(FT.relDiff_P2x(k)), max(FT.relDiff_P2x(k)), sum(FT.fixP2x));
end
fprintf(f, '   fitting time per soil: median %.1f s, max %.1f s\n', median(FT.tFit_s), max(FT.tFit_s));
end

function make_figures(outDir, S, FT, D, ids, grp)
% Figure S6: endpoints and order of continuity, KIT-022 (P2)
i = find(ids == "KIT_022", 1);
if ~isempty(i)
    p = pvec(FT, 'P2', i); b = 1 + S.c/p(2);
    f7 = figure('Visible', 'off', 'Position', [100 100 1450 430]);
    hh = logspace(log10(S.hLow), log10(S.h0) - 1e-9, 4000)';
    subplot(1, 3, 1);
    loglog(hh, evm_capacity(hh, p, b, S.h0), 'r-', hh, evm_capacity(hh, p, 1, S.h0), 'b--', 'LineWidth', 1.4); hold on
    loglog(S.h0, evm_capacity(S.h0, p, b, S.h0), 'ko', 'MarkerFaceColor', 'k');
    xline(S.hmin, 'k:'); grid on; xlim([S.hLow 3*S.h0]);
    xlabel('h (cm)'); ylabel('C(h) (cm^{-1})');
    legend(sprintf('\\beta = 1 + 0.1/n = %.3f: C(0^+) = 0', b), '\beta = 1: C(0^+) > 0', 'C(h_0) > 0 (dry end)', 'Location', 'southwest');
    title(sprintf('(a) KIT-022: C > 0, C \\propto h^{%.3f} as h \\rightarrow 0', min(p(2), b) - 1));
    subplot(1, 3, 2);
    pf = linspace(4, log10(S.h0), 400)';
    plot(pf, evm_theta(10.^pf, p, b, S.h0), 'r-', 'LineWidth', 1.6); hold on
    s = log(10)*S.h0*evm_capacity(S.h0, p, b, S.h0);
    plot([6 log10(S.h0)], s*[1 0], 'k--'); plot([log10(S.h0) 7.3], [0 0], 'k:', 'LineWidth', 1.6);
    xlabel('pF = log_{10}(h/cm)'); ylabel('\theta (cm^3 cm^{-3})'); grid on; xlim([4 7.3]);
    legend('Eq. (7)', sprintf('slope at h_0: -d\\theta/dpF = %.3f', s), '\theta = 0 beyond h_0 (outside Eq. 7)', 'Location', 'northeast');
    title('(b) dry end: \theta(h_0) = 0 with a finite slope');
    subplot(1, 3, 3);
    hw = logspace(log10(S.hLow), 0, 2000)';
    loglog(hw, evm_dCdh(hw, p, b, S.h0), 'r-', 'LineWidth', 1.6); grid on; xlim([S.hLow 1]);
    xlabel('h (cm)'); ylabel('dC/dh = -d^2\theta/dh^2 (cm^{-2})');
    title(sprintf('(c) dC/dh \\propto h^{%.3f} \\rightarrow \\infty: \\theta is C^1, not C^2', min(p(2), b) - 2));
    exportgraphics(f7, fullfile(outDir, 'FigS6_endpoints_continuity.png'), 'Resolution', 300); close(f7);
end

% Figure S7: capacity functions of all constrained fits
f8 = figure('Visible', 'off', 'Position', [100 100 1450 820]);
gl = ["TUDD" "UFZ" "FSUJ" "KIT" "TUBS"];
hh = logspace(-3, log10(S.h0) - 1e-9, 600)';
for k = 1:numel(gl)
    subplot(2, 3, k); hold on
    sel = find(grp == gl(k));
    for i = sel'
        p = pvec(FT, 'P2', i);
        plot(hh, evm_capacity(hh, p, 1 + S.c/p(2), S.h0), 'LineWidth', 0.5);
    end
    set(gca, 'XScale', 'log', 'YScale', 'log'); grid on; box on; xlim([1e-3 S.h0]); ylim([1e-14 1]);
    xlabel('h (cm)'); ylabel('C(h) (cm^{-1})');
    title(sprintf('(%s) %s, N = %d: C > 0 for all', char('a' + k - 1), gl(k), numel(sel)));
end
subplot(2, 3, 6);
histogram(log10(FT.P2_e0), 30); grid on
xlabel('log_{10} \epsilon, with C \propto (\alpha h)^{\epsilon} as h \rightarrow 0^+'); ylabel('soils');
title('(f) wet-end exponent \epsilon = min(n, \beta) - 1');
exportgraphics(f8, fullfile(outDir, 'FigS7_capacity_all_soils.png'), 'Resolution', 300); close(f8);
end
