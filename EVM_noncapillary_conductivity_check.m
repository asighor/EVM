function R = EVM_noncapillary_conductivity_check(dataDir, nSoilsMax)
%EVM_NONCAPILLARY_CONDUCTIVITY_CHECK  Checks the derivation of the EVM non-capillary
% conductivity (Eqs. 18-20) against the exact Mualem integral, and quantifies the error of
% the power-law form near the dry endpoint, analytically and for the 572 German soils of
% Hohenbrink et al. (2023): RetMeas.csv and CondMeas.csv.
%
%   R = EVM_noncapillary_conductivity_check              % CSV files next to this file
%   R = EVM_noncapillary_conductivity_check(dataDir)     % CSV files in dataDir
%   R = EVM_noncapillary_conductivity_check(dataDir, 20) % quick test on the first 20 soils
%
% Text S2 of the Supporting Information (derivation and accuracy of Eq. 20).
% Notation: S = Snc of Eq. (6), R = h0/ha = alpha*h0, A = R^beta, p = gamma/beta, and the
% normalized Mualem integral M(S) = int_0^S h(x)^-gamma dx / int_0^1 h(x)^-gamma dx, so that
% Knc = Ks,nc S^l M(S)^eta (Eq. 12 applied to the non-capillary domain).
%   M19(S) = (R^(gamma S) - 1)/(R^gamma - 1)         Eq. (19) inserted, integral exact (Eq. B2)
%   M20(S) = R^(-gamma (1 - S))                       power law of Eq. (20): R^gamma >> 1 and
%                                                     R^(gamma S) >> 1 in addition
%   Mex(S) = [I(tS) - I(t0)]/[1 - I(t0)]              exact inversion, Eq. (18): I = regularized
%            tS = 1/(1 + (h/ha)^beta), t0 = 1/(1+A)   incomplete beta function I_t(p, 1 - p)
% For h >> ha, Mex = M19/(p B(p, 1 - p)) with p B(p, 1 - p) = p pi/sin(p pi).
% Parts:
%   A. Closed forms against adaptive numerical quadrature on a grid of R and beta.
%   B. Error factors of the power law, M20/M19 = (1 - R^-gamma)/(1 - R^(-gamma S)), and of the
%      approximate inversion, M19/Mex, as functions of pF.
%   C. The 572 EVM fits (theta_s > theta_r >= 0, beta = 1 + 0.1/n; procedure of Script S2):
%      conductivity with the power law (Eq. 20), with M19 (Eq. B2, vanishes at h0) and with
%      Mex (same Ks,nc); RMSE of log10 K, changes at the measured suctions and toward h0.
%   D. The five Table 1 literature soils (Table1_literature_RetMeas.csv and
%      Table1_literature_CondMeas.csv, written by EVM_literature_data.m; skipped if absent), whose
%      conductivities extend to pF 5.3: the same comparison at their measured suctions.
% l = 0, eta = 2, gamma = 0.75, Ks,nc = 86.2 theta(1e5 cm) ha^-1.5, h0 = 1e7 cm.
%
% Outputs (folder ncond_results next to the data):
%   EVM_ncond_report.txt     all numbers of Text S2 and Appendix B
%   EVM_ncond_grid.csv       Part A/B: closed forms, quadrature and error factors on the grid
%   EVM_ncond_suction_errors.csv  Table S3: error factors at fixed suctions (R = 1e3, 1e5 and 1e7)
%   EVM_ncond_by_soil.csv    per soil: RMSE_K for the three forms and the dry-end diagnostics
%                            (column suffixes: _Eq20 = power law, Eq. 20; _M19 = integral with Eq. 19;
%                            _exact = exact inversion)
%   FigS5_noncapillary_conductivity.png, EVM_ncond_results.mat
% Runtime: about 3 min with 24 parallel workers.

if nargin < 1 || isempty(dataDir)
    dataDir = fileparts(mfilename('fullpath'));
end
outDir = fullfile(dataDir, 'ncond_results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
tStart = tic;

%% ---------------- settings ----------------
S.h0 = 1e7; S.c = 0.1; S.gam = 0.75; S.eta = 2; S.l = 0;
S.fMax = 0.9999; S.nMin = 1.0001;
S.pFgrid = [3 3.93 4 4.5 5 5.5 6 6.5 6.8 6.9 6.99];        % suctions of the error tables
S.pFK = [3.93 5 6 6.5 6.9 log10(1e7) - 1e-9];               % conductivity reported here
S.eps = [0.01 0.1 1 9];                                     % power law above M19 by 1%, 10%, x2, x10
opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxFunctionEvaluations', 4000, 'MaxIterations', 1000, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-10);
optsLC = optimoptions('lsqcurvefit', 'Algorithm', 'interior-point', 'Display', 'off', ...
    'MaxFunctionEvaluations', 8000, 'MaxIterations', 2000, 'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-10, 'ConstraintTolerance', 1e-9);

%% ---------------- A, B. closed forms vs quadrature, error factors ----------------
Rg = 10.^(2:7); bg = [1.0100 1.0500 1.0999];
Sg = [1e-4 1e-3 0.01 0.05 0.1 0.2 0.3 0.5 0.7 0.9 0.99 1];
rows = cell(0, 12);
for Rv = Rg
    for b = bg
        for s = Sg
            [m19, m20, mex] = mualem_forms(s, Rv, b, S);
            [n19, nex] = mualem_quadrature(s, Rv, b, S);
            p = S.gam/b;
            rows(end + 1, :) = {Rv, b, s, snc_to_pF(s, Rv, b, S), m19, n19, mex, nex, m20, ...
                (m20/m19)^S.eta, (m19/mex)^S.eta, (p*pi/sin(p*pi))^S.eta}; %#ok<AGROW>
        end
    end
end
GA = cell2table(rows, 'VariableNames', {'R', 'beta', 'Snc', 'pF', 'M19_closed', 'M19_quad', 'Mex_closed', ...
    'Mex_quad', 'M20', 'K20_over_K19', 'K19_over_Kex', 'pB_eta'});
GA.relErr19 = abs(GA.M19_closed - GA.M19_quad)./max(GA.M19_quad, realmin);
GA.relErrEx = abs(GA.Mex_closed - GA.Mex_quad)./max(GA.Mex_quad, realmin);

% error factors at fixed suctions for three representative R
rowsB = cell(0, 8);
for Rv = [1e3 1e5 1e7]
    for pf = S.pFgrid
        h = 10^pf; b = 1.05;
        if h >= S.h0, continue, end
        s = snc(h, Rv, b, S);
        [m19, m20, mex] = mualem_forms(s, Rv, b, S);
        rowsB(end + 1, :) = {Rv, b, pf, s, (m20/m19)^S.eta, (1 - (h/S.h0)^S.gam)^(-S.eta), (m19/mex)^S.eta, ...
            (m20/mex)^S.eta}; %#ok<AGROW>
    end
end
GB = cell2table(rowsB, 'VariableNames', {'R', 'beta', 'pF', 'Snc', 'K20_over_K19', 'approx_1_minus_h_over_h0', ...
    'K19_over_Kex', 'K20_over_Kex'});

%% ---------------- C. the fitted soils ----------------
ret  = readtable(fullfile(dataDir, 'RetMeas.csv'),  'TextType', 'string');
cond = readtable(fullfile(dataDir, 'CondMeas.csv'), 'TextType', 'string');
ret  = ret(isfinite(ret.pF) & isfinite(ret.theta), :);
ids  = unique(ret.Sample_ID, 'stable');
if nargin > 1 && ~isempty(nSoilsMax)
    if isnumeric(nSoilsMax), ids = ids(1:min(nSoilsMax, numel(ids)));
    else,                    ids = ids(ismember(ids, string(nSoilsMax)));
    end
end
nS = numel(ids);
[~, ir] = ismember(ret.Sample_ID, ids);
[~, ic] = ismember(cond.Sample_ID, ids);
D = struct('id', cellstr(ids), 'h', [], 'th', [], 'hK', [], 'K', [], 'Ksat', NaN);
for i = 1:nS
    D(i).h = 10.^ret.pF(ir == i); D(i).th = ret.theta(ir == i);
    isKs  = ic == i & cond.MeasType == "KSAT";
    unsat = ic == i & cond.MeasType ~= "KSAT" & cond.k > 0;
    D(i).hK = 10.^cond.pF(unsat); D(i).K = cond.k(unsat);
    if any(isKs), D(i).Ksat = cond.k(find(isKs, 1)); end
end
fprintf('Loaded %d soils, %d with measured KSAT.\n', nS, sum(isfinite([D.Ksat])));
if isempty(gcp('nocreate'))          % parfor runs serially when no pool can start
    try
        pc = parcluster('Processes');
        pc.NumWorkers = max(pc.NumWorkers, feature('numcores'));
        parpool(pc, feature('numcores'));
    catch
        try, parpool('Processes'); catch, end %#ok<NOCOM>
    end
end
out = cell(nS, 1);
parfor i = 1:nS
    out{i} = analyse_soil(D(i), S, opts, optsLC);
end
T = struct2table(cell2mat(cellfun(@(o) o.row, out, 'UniformOutput', false)));
T = [table(ids, extractBefore(ids, "_"), [D.Ksat]', 'VariableNames', {'Sample_ID', 'Group', 'Ksat_cm_d'}), T];
dK = cell2mat(cellfun(@(o) o.dlogK, out, 'UniformOutput', false));    % [19 - 20, ex - 20] at the K points
pFK = cell2mat(cellfun(@(o) o.pFK, out, 'UniformOutput', false));

chk = struct('n', 0, 'd', NaN, 'dp', NaN, 'worst', "");
f2 = fullfile(dataDir, 'positivity_results', 'EVM_positivity_fits.csv');
if exist(f2, 'file')
    F2 = readtable(f2, 'TextType', 'string');
    [tf, j] = ismember(ids, F2.Sample_ID);
    k = tf & isfinite(T.RMSEK_20);
    if any(k), chk.n = sum(k); chk.d = max(abs(T.RMSEK_20(k) - F2.P2_RMSE_logK(j(k)))); end
    a = [T.alpha T.n T.theta_r T.theta_s]; b2 = nan(size(a));
    b2(tf, :) = [F2.P2_alpha(j(tf)) F2.P2_n(j(tf)) F2.P2_theta_r(j(tf)) F2.P2_theta_s(j(tf))];
    rd = max(abs(a - b2)./max(abs(b2), 1e-6), [], 2);
    [chk.dp, w] = max(rd); chk.worst = ids(w);
end

%% ---------------- D. the Table 1 literature soils ----------------
TL = table();
fl = fullfile(dataDir, 'Table1_literature_RetMeas.csv'); fk = fullfile(dataDir, 'Table1_literature_CondMeas.csv');
if exist(fl, 'file') && exist(fk, 'file') && (nargin < 2 || isempty(nSoilsMax))
    rl = readtable(fl, 'TextType', 'string'); kl = readtable(fk, 'TextType', 'string');
    rl = rl(isfinite(rl.pF) & isfinite(rl.theta), :);
    lid = unique(kl.Sample_ID(kl.MeasType ~= "KSAT"), 'stable');
    rowsL = cell(numel(lid), 1);
    for i = 1:numel(lid)
        u = kl.Sample_ID == lid(i) & kl.MeasType ~= "KSAT" & kl.k > 0;
        s = kl.Sample_ID == lid(i) & kl.MeasType == "KSAT";
        d = struct('id', char(lid(i)), 'h', 10.^rl.pF(rl.Sample_ID == lid(i)), ...
            'th', rl.theta(rl.Sample_ID == lid(i)), 'hK', 10.^kl.pF(u), 'K', kl.k(u), 'Ksat', NaN);
        if any(s), d.Ksat = kl.k(find(s, 1)); end
        o = analyse_soil(d, S, opts, optsLC);
        rowsL{i} = table(lid(i), numel(d.K), max(log10(d.hK)), o.row.RMSEK_20, o.row.RMSEK_19, ...
            o.row.maxAbsDlogK_19, 'VariableNames', {'Sample_ID', 'nK', 'maxpF_K', 'RMSEK_Eq20', ...
            'RMSEK_M19', 'maxAbsDlogK_M19'});
    end
    TL = vertcat(rowsL{:});
end

%% ---------------- outputs ----------------
writetable(GA, fullfile(outDir, 'EVM_ncond_grid.csv'));
writetable(GB, fullfile(outDir, 'EVM_ncond_suction_errors.csv'));
vn = string(T.Properties.VariableNames);                     % file labels of the revised manuscript
vn = regexprep(vn, '_20(?=_|$)', '_Eq20'); vn = regexprep(vn, '_19(?=_|$)', '_M19');
vn = regexprep(vn, '_ex(?=_|$)', '_exact'); vn = replace(vn, '_powerlaw', '_Eq20');
Tw = T; Tw.Properties.VariableNames = cellstr(vn);
writetable(Tw, fullfile(outDir, 'EVM_ncond_by_soil.csv'));
fid = fopen(fullfile(outDir, 'EVM_ncond_report.txt'), 'w');
if fid < 0, error('Cannot write the report in %s', outDir); end
for f = [1, fid]
    print_report(f, S, GA, GB, T, dK, pFK, chk, bg, TL);
end
fclose(fid);
make_figure(outDir, S, T, D, ids, out);
R = struct('settings', S, 'grid', GA, 'suctionErrors', GB, 'bySoil', T, 'table1Soils', TL);
save(fullfile(outDir, 'EVM_ncond_results.mat'), 'R');
fprintf('Outputs written to %s (%.1f min).\n', outDir, toc(tStart)/60);
fprintf('EVM_NCOND_CHECK_DONE\n');
end

%% =====================================================================
%% local functions
%% =====================================================================
function s = snc(h, Rv, b, S)
% Eq. (6) with x = h/ha = h R/h0 and x0 = h0/ha = R
s = 1 - log1p((h*Rv/S.h0).^b)/log1p(Rv^b);
end

function pf = snc_to_pF(s, Rv, b, S)
% exact inversion, Eq. (18), expressed as pF
ha = S.h0/Rv;
pf = log10(ha*((1 + Rv^b).^(1 - s) - 1).^(1/b));
end

function [m19, m20, mex] = mualem_forms(s, Rv, b, S)
% normalized Mualem integrals of the non-capillary domain (closed forms)
g = S.gam; p = g/b; A = Rv^b;
m19 = expm1(g*s*log(Rv))/expm1(g*log(Rv));                 % (R^(gS) - 1)/(R^g - 1)
m20 = Rv.^(-g*(1 - s));                                    % power law
t0 = 1/(1 + A); tS = (1 + A).^(-(1 - s));                  % tS = 1/(1 + (h/ha)^beta)
I0 = betainc(t0, p, 1 - p);
mex = (betainc(tS, p, 1 - p) - I0)/(1 - I0);
end

function [n19, nex] = mualem_quadrature(s, Rv, b, S)
% the same integrals by adaptive quadrature in y = 1 - x; the exact integrand has the weak
% singularity y^-p at y = 0 (h -> 0), removed by y = z^q with q = 1/(1 - p)
g = S.gam; p = g/b; A = Rv^b; q = 1/(1 - p);
f19 = @(y) Rv.^(-g*y);                                     % (h/ha)^-gamma with Eq. (19)
fex = @(z) (expm1(log1p(A)*z.^q)).^(-p).*q.*z.^(q - 1);    % [(1+A)^y - 1]^-p dy, Eq. (18)
o = {'RelTol', 1e-10, 'AbsTol', 1e-300};
ws = warning; warning('off', 'MATLAB:integral:MaxIntervalCountReached'); warning('off', 'MATLAB:integral:MinStepSize');
n19 = integral(f19, 1 - s, 1, o{:})/integral(f19, 0, 1, o{:});
nex = integral(fex, (1 - s)^(1/q), 1, o{:})/integral(fex, 0, 1, o{:});
warning(ws);
end

function o = analyse_soil(d, S, opts, optsLC)
% EVM fit (theta_s > theta_r >= 0, as Script S2) and conductivity with the three forms
p = fit_soil(d, S, opts, optsLC);
al = p(1); n = p(2); thr = p(3); ths = p(4); b = 1 + S.c/n; Rv = al*S.h0;
r.alpha = al; r.n = n; r.theta_r = thr; r.theta_s = ths; r.beta = b; r.R = Rv;
pp = S.gam/b; r.pB_eta = (pp*pi/sin(pp*pi))^S.eta;
r.Knc_h0_powerlaw = 86.2*evm_theta(1e5, p, b, S.h0)*al^1.5*Rv^(-S.gam*S.eta);   % Eq. (20) at h0
r.Kv_h0 = vapor_K(S.h0, 0, ths);
% pF at which the power law exceeds M19 by 1%, 10%, x2 and x10
for e = 1:numel(S.eps)
    s = -log(1 - (1 - Rv^(-S.gam))*(1 + S.eps(e))^(-1/S.eta))/(S.gam*log(Rv));
    r.(sprintf('pF_powerlaw_err_%g', S.eps(e)*100)) = log10(S.h0/Rv*((1 + Rv^b)^(1 - s) - 1)^(1/b));
end
% quadrature check of the exact form at three suctions
er = 0;
for pf = [3 5 6.5]
    s = snc(10^pf, Rv, b, S);
    [~, ~, mex] = mualem_forms(s, Rv, b, S); [~, nex] = mualem_quadrature(s, Rv, b, S);
    er = max(er, abs(mex - nex)/nex);
end
r.maxRelErr_exact_closed_vs_quad = er;
r.RMSEK_20 = NaN; r.RMSEK_19 = NaN; r.RMSEK_ex = NaN; r.maxAbsDlogK_19 = NaN;
lk = nan(3, numel(S.pFK));
o.dlogK = zeros(0, 2); o.pFK = zeros(0, 1);
if isfinite(d.Ksat) && ~isempty(d.K)
    L = zeros(numel(d.hK), 3);
    for k = 1:3, L(:, k) = log10(evm_K(d.hK, p, b, d.Ksat, S, k)); end
    e = L - log10(d.K);
    r.RMSEK_20 = sqrt(mean(e(:, 1).^2)); r.RMSEK_19 = sqrt(mean(e(:, 2).^2)); r.RMSEK_ex = sqrt(mean(e(:, 3).^2));
    r.maxAbsDlogK_19 = max(abs(L(:, 2) - L(:, 1)));
    o.dlogK = [L(:, 2) - L(:, 1), L(:, 3) - L(:, 1)]; o.pFK = log10(d.hK);
    for k = 1:3, lk(k, :) = log10(evm_K(10.^S.pFK(:), p, b, d.Ksat, S, k))'; end
end
nm = ["20" "19" "ex"];
for k = 1:3
    for j = 1:numel(S.pFK)
        r.(sprintf('log10K_%s_pF%s', nm(k), strrep(sprintf('%.2f', min(S.pFK(j), 7)), '.', 'p'))) = lk(k, j);
    end
end
o.row = r; o.p = p;
end

function K = evm_K(h, p, beta, Ksc, S, form)
% K = Kc + Knc + Kv (Eq. 21 with Ks,c counted once); Knc with form 1 = power law (Eq. 20),
% 2 = M19 (Eq. B2, vanishes at h0), 3 = Mex (exact inversion), all with the same Ks,nc
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n; h0 = S.h0; Rv = al*h0;
G   = (1 + (al*h).^n).^(-m);
G0  = (1 + (al*h0)^n)^(-m);
Sc  = max((G - G0)/(1 - G0), 0);
Snc = max(1 - log1p((al*h).^beta)/log1p((al*h0)^beta), 0);
Kc  = Ksc*sqrt(Sc).*(1 - ((1 - G.^(1/m))/(1 - G0^(1/m))).^m).^2;   % Eq. (14)
Ksnc = 86.2*evm_theta(1e5, p, beta, h0)*al^1.5;                     % Eq. (17), ha = 1/alpha
switch form
    case 1, M = Rv.^(-S.gam*(1 - Snc));
    case 2, M = expm1(S.gam*Snc*log(Rv))/expm1(S.gam*log(Rv));
    case 3
        pp = S.gam/beta; A = Rv^beta; I0 = betainc(1/(1 + A), pp, 1 - pp);
        M = (betainc(1./(1 + (al*h).^beta), pp, 1 - pp) - I0)/(1 - I0);
end
Knc = Ksnc*Snc.^S.l.*max(M, 0).^S.eta;
th  = (ths - thr)*Sc + thr*Snc;
Kv  = vapor_K(h, th, ths);
K   = max(Kc + Knc + Kv, realmin);
end

function Kv = vapor_K(h, th, ths)
% isothermal vapour conductivity (cm/d), Appendix A, 20 C
T = 293.15; Mw = 0.018015; g = 9.81; Rg = 8.314;
rho_sv = 1e-3*exp(31.3716 - 6014.79/T - 7.92495e-3*T)/T;
rho_w  = 998.2;
Da  = 2.12e-5*(T/273.15)^2*1e4*86400;
cH  = Mw*g/(Rg*T)/100;
tau = max(ths - th, 0).^(10/3)/ths^2;
Kv  = rho_sv/rho_w*tau*Da*cH.*exp(-cH*h);
end

function th = evm_theta(h, p, beta, h0)
% Eq. (7), written as in Script S2 so that the fits are identical
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n;
G   = (1 + (al*h).^n).^(-m);
G0  = (1 + (al*h0)^n)^(-m);
Sc  = (G - G0)/(1 - G0);
Snc = 1 - log1p((al*h).^beta)/log1p((al*h0)^beta);
th  = (ths - thr)*Sc + thr*Snc;
end

function p = fit_soil(d, S, opts, optsLC)
% the constrained EVM fit of Script S2 (P1 and P2 aligned, P2x check where active)
[p1, s1] = fit_evm(d.h, d.th, 'P1', start_set('P1', d.th, false), S, opts, optsLC);
[p2, s2] = fit_evm(d.h, d.th, 'P2', start_set('P2', d.th, false), S, opts, optsLC);
if p1(3) > S.fMax*p1(4)
    [pb, sb] = fit_evm(d.h, d.th, 'P1', start_set('P1', d.th, true), S, opts, optsLC);
    if sb < s1, p1 = pb; s1 = sb; end
    [pb, sb] = fit_evm(d.h, d.th, 'P2', start_set('P2', d.th, true), S, opts, optsLC);
    if sb < s2, p2 = pb; s2 = sb; end
end
if s2 < s1*(1 - 1e-9), p1 = p2; s1 = s2; end
act = p1(3) > S.fMax*p1(4);
if ~act && s1 < s2*(1 - 1e-9), p2 = p1; s2 = s1; end
if act
    [px, sx] = fit_evm(d.h, d.th, 'P2x', start_set('P2x', d.th, false), S, opts, optsLC);
    if sx < s2*(1 - 1e-9) && px(3) < px(4), p2 = px; end
end
p = p2;
end

function [p, sseBest] = fit_evm(h, th, pr, Q0, S, opts, optsLC)
switch pr
    case 'P1',  lb = [-5 S.nMin 0 0];    ub = [0 10 1 1];
    case 'P2',  lb = [-5 S.nMin 0.01 0]; ub = [0 10 1 S.fMax];
    case 'P2x', lb = [-5 S.nMin 0 0.01]; ub = [0 10 1 1];
end
fun = @(q, hh) evm_theta(hh, q2p(q, pr), 1 + S.c/q(2), S.h0);
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

function print_report(f, S, GA, GB, T, dK, pFK, chk, bg, TL)
fprintf(f, '\n========== EVM non-capillary conductivity: Eqs. (18)-(20) vs the exact Mualem integral ==========\n');
fprintf(f, 'gamma = %.2f, eta = %d, l = %d, h0 = %.0e cm; R = h0/ha = alpha h0, p = gamma/beta\n', S.gam, S.eta, S.l, S.h0);
fprintf(f, 'M19 = (R^(gS) - 1)/(R^g - 1) | power law M20 = R^(-g(1-S)) | Mex = [I(tS) - I(t0)]/[1 - I(t0)], I = betainc(t, p, 1-p)\n');

fprintf(f, '\n--- A. Closed forms vs adaptive quadrature (R = 1e2-1e7, beta = %s, %d values of Snc) ---\n', mat2str(bg), numel(unique(GA.Snc)));
fprintf(f, '   max relative difference: M19 %.1e | Mex %.1e\n', max(GA.relErr19), max(GA.relErrEx));
for b = bg
    p = S.gam/b;
    fprintf(f, '   beta = %.4f: p = %.4f, p B(p,1-p) = %.3f, (p B)^eta = %.2f (%.2f log10 units)\n', b, p, p*pi/sin(p*pi), ...
        (p*pi/sin(p*pi))^S.eta, S.eta*log10(p*pi/sin(p*pi)));
end

fprintf(f, '\n--- B. Error factors in Knc (beta = 1.05) ---\n');
fprintf(f, '   power law / M19 = [(1 - R^-g)/(1 - R^(-g S))]^eta ~ [1 - (h/h0)^g]^-eta | M19 / Mex -> (p B)^eta for h >> ha\n');
for Rv = unique(GB.R)'
    g = GB(GB.R == Rv, :);
    fprintf(f, '   R = %.0e (alpha = %.0e cm-1):\n', Rv, Rv/S.h0);
    for i = 1:height(g)
        fprintf(f, '      pF %.2f  Snc %.4f | power law/M19 %8.4f (approx %8.4f) | M19/Mex %6.3f | power law/Mex %8.3f\n', ...
            g.pF(i), g.Snc(i), g.K20_over_K19(i), g.approx_1_minus_h_over_h0(i), g.K19_over_Kex(i), g.K20_over_Kex(i));
    end
end
fprintf(f, '   general (R^g >> 1): power law exceeds M19 by 1%% at h0/h = %.0f (pF %.2f), 10%% at h0/h = %.1f (pF %.2f), x2 at %.2f (pF %.2f), x10 at %.2f (pF %.2f)\n', ...
    thr_ratio(0.01, S), 7 - log10(thr_ratio(0.01, S)), thr_ratio(0.1, S), 7 - log10(thr_ratio(0.1, S)), thr_ratio(1, S), ...
    7 - log10(thr_ratio(1, S)), thr_ratio(9, S), 7 - log10(thr_ratio(9, S)));

fprintf(f, '\n--- C. The %d EVM fits (theta_s > theta_r >= 0) ---\n', height(T));
k = isfinite(T.RMSEK_20);
fprintf(f, '   mean RMSE of log10 K (%d soils with KSAT): power law (Eq. 20) %.4f | M19 (Eq. B2), vanishes at h0 %.4f | Mex, same Ks,nc %.4f\n', ...
    sum(k), mean(T.RMSEK_20(k)), mean(T.RMSEK_19(k)), mean(T.RMSEK_ex(k)));
fprintf(f, '   M19 vs power law at the %d measured K points (pF %.2f-%.2f): max |d log10 K| %.2e, max per-soil |d RMSE_K| %.2e\n', ...
    size(dK, 1), min(pFK), max(pFK), max(abs(dK(:, 1))), max(abs(T.RMSEK_19(k) - T.RMSEK_20(k))));
fprintf(f, '   Mex vs power law at the measured K points: d log10 K median %+.3f (5-95%%: %+.3f to %+.3f)\n', ...
    median(dK(:, 2)), quantile(dK(:, 2), 0.05), quantile(dK(:, 2), 0.95));
fprintf(f, '   (p B)^eta over the fits: median %.2f (range %.2f-%.2f)\n', median(T.pB_eta), min(T.pB_eta), max(T.pB_eta));
for e = S.eps
    v = T.(sprintf('pF_powerlaw_err_%g', e*100));
    fprintf(f, '   pF at which the power law exceeds M19 by %5g%%: median %.2f (range %.2f-%.2f)\n', e*100, median(v), min(v), max(v));
end
fprintf(f, '   at h0: power-law Knc median %.2e cm/d, Kv median %.2e cm/d, power-law share of K %.2f (median)\n', ...
    median(T.Knc_h0_powerlaw), median(T.Kv_h0), median(T.Knc_h0_powerlaw./(T.Knc_h0_powerlaw + T.Kv_h0)));
nm = ["20" "19" "ex"]; lab = ["power law" "M19      " "Mex      "];
for k2 = 1:3
    fprintf(f, '   median log10 K, %s:', lab(k2));
    for j = 1:numel(S.pFK)
        v = T.(sprintf('log10K_%s_pF%s', nm(k2), strrep(sprintf('%.2f', min(S.pFK(j), 7)), '.', 'p')));
        fprintf(f, ' pF %.2f %.2f |', min(S.pFK(j), 7), median(v, 'omitnan'));
    end
    fprintf(f, '\n');
end
fprintf(f, '   quadrature check of Mex for every fit (pF 3, 5, 6.5): max relative difference %.1e\n', max(T.maxRelErr_exact_closed_vs_quad));
if chk.n > 0
    fprintf(f, '   power-law RMSE_K vs Script S2 (%d soils): max difference %.1e | parameters: max relative difference %.1e (%s)\n', ...
        chk.n, chk.d, chk.dp, chk.worst);
end
if ~isempty(TL)
    fprintf(f, '\n--- D. Table 1 literature soils with conductivity data (theta_s > theta_r >= 0) ---\n');
    for i = 1:height(TL)
        fprintf(f, '   %-13s %3d K data to pF %.2f | RMSE of log10 K: power law (Eq. 20) %.4f, M19 %.4f | max |d log10 K| %.2e\n', ...
            TL.Sample_ID(i), TL.nK(i), TL.maxpF_K(i), TL.RMSEK_Eq20(i), TL.RMSEK_M19(i), TL.maxAbsDlogK_M19(i));
    end
    fprintf(f, '   max |d RMSE_K| %.2e | max |d log10 K| %.2e\n', max(abs(TL.RMSEK_M19 - TL.RMSEK_Eq20)), max(TL.maxAbsDlogK_M19));
end
end

function x = thr_ratio(e, S)
% h0/h at which [1 - (h/h0)^gamma]^-eta = 1 + e
x = (1 - (1 + e)^(-1/S.eta))^(-1/S.gam);
end

function make_figure(outDir, S, T, D, ids, out)
col = [0.15 0.45 0.80; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56];
f14 = figure('Visible', 'off', 'Position', [100 100 1450 900]);
subplot(2, 2, 1);
Rv = 1e5; b = 1.05; s = logspace(-5, 0, 400)';
M = zeros(numel(s), 3);
for i = 1:numel(s), [M(i, 1), M(i, 2), M(i, 3)] = mualem_forms(s(i), Rv, b, S); end
loglog(s, M(:, 2), '-', 'Color', col(2, :), 'LineWidth', 1.8); hold on
loglog(s, M(:, 1), '-', 'Color', col(1, :), 'LineWidth', 1.8);
loglog(s, M(:, 3), '--', 'Color', col(3, :), 'LineWidth', 1.8);
xline(1/(S.gam*log(Rv)), 'k:', 'S_{nc} = 1/(\gamma ln R)', 'LabelVerticalAlignment', 'top', ...
    'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
grid on; xlabel('S_{nc}'); ylabel('normalized Mualem integral M(S_{nc})');
legend('Eq. (20): R^{-\gamma(1-S)}', 'M_{19} = (R^{\gammaS}-1)/(R^\gamma-1)', ...
    'M_{ex}: exact inversion, Eq. (18)', 'Location', 'southeast');
title(sprintf('(a) R = h_0/h_a = 10^5, \\beta = %.2f, \\gamma = %.2f', b, S.gam));
subplot(2, 2, 2); hold on
pf = linspace(2, 6.999, 500)';
for j = 1:3
    Rv = 10^(2*j + 1); r1 = zeros(size(pf)); r2 = r1;
    for i = 1:numel(pf)
        sv = snc(10^pf(i), Rv, b, S); [m19, m20, mex] = mualem_forms(sv, Rv, b, S);
        r1(i) = (m20/m19)^S.eta; r2(i) = (m19/mex)^S.eta;
    end
    v = pf >= log10(S.h0/Rv);
    if j == 2
        plot(pf(v), r1(v), 'k-', 'LineWidth', 1.8, 'DisplayName', 'Eq. (20) / M_{19} (same for every R)');
    end
    plot(pf(v), r2(v), '--', 'Color', col(j, :), 'LineWidth', 1.4, 'DisplayName', sprintf('M_{19} / M_{ex}, R = 10^%d', 2*j + 1));
end
xline(3.93, 'k-', 'driest K datum', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
set(gca, 'YScale', 'log'); ylim([0.9 1e3]); xlim([2 7]); grid on; box on
xlabel('pF'); ylabel('ratio of K_{nc}'); legend('Location', 'northwest', 'FontSize', 7);
title(sprintf('(b) error factors in K_{nc} (\\eta = %d)', S.eta));
subplot(2, 2, 3);
labs = strings(0, 1); dv = [];
for j = 1:numel(S.pFK)
    tg = strrep(sprintf('%.2f', min(S.pFK(j), 7)), '.', 'p');
    v = T.(sprintf('log10K_19_pF%s', tg)) - T.(sprintf('log10K_20_pF%s', tg)); v = v(isfinite(v));
    labs = [labs; repmat(string(sprintf('%.2f', min(S.pFK(j), 7))), numel(v), 1)]; dv = [dv; v]; %#ok<AGROW>
end
boxchart(categorical(labs, unique(labs, 'stable')), dv, 'BoxFaceColor', col(1, :), 'MarkerStyle', '.');
yline(0, 'k-'); grid on; xlabel('pF'); ylabel('\Delta log_{10} K, M_{19} minus Eq. (20)');
title(sprintf('(c) total K of the %d fits with K_{sat}: change toward oven dryness', sum(isfinite(T.RMSEK_20))));
subplot(2, 2, 4);
i = find(ids == "KIT_022", 1);
if isempty(i), i = find(isfinite([D.Ksat]), 1); end
d = D(i); p = out{i}.p; b = 1 + S.c/p(2); hh = logspace(-2, log10(S.h0) - 1e-9, 600)';
al = p(1); Rv = al*S.h0; Ksnc = 86.2*evm_theta(1e5, p, b, S.h0)*al^1.5;
Snc = max(1 - log1p((al*hh).^b)/log1p(Rv^b), 0);
pp = S.gam/b; I0 = betainc(1/(1 + Rv^b), pp, 1 - pp);
Knc20 = Ksnc*Rv.^(-S.gam*S.eta*(1 - Snc)); Knc19 = Ksnc*(expm1(S.gam*Snc*log(Rv))/expm1(S.gam*log(Rv))).^S.eta;
Kncex = Ksnc*max((betainc(1./(1 + (al*hh).^b), pp, 1 - pp) - I0)/(1 - I0), 0).^S.eta;
th = evm_theta(hh, p, b, S.h0); Kv = vapor_K(hh, th, p(4));
semilogy(log10(hh), Knc20, '-', 'Color', col(2, :), 'LineWidth', 1.6); hold on
semilogy(log10(hh), Knc19, '-', 'Color', col(1, :), 'LineWidth', 1.6);
semilogy(log10(hh), Kncex, '--', 'Color', col(3, :), 'LineWidth', 1.4);
semilogy(log10(hh), Kv, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.4);
if isfinite(d.Ksat)
    semilogy(log10(hh), evm_K(hh, p, b, d.Ksat, S, 1), 'k-', 'LineWidth', 0.8);
    semilogy(log10(hh), evm_K(hh, p, b, d.Ksat, S, 2), 'k--', 'LineWidth', 0.8);
    semilogy(log10(d.hK), d.K, 'ko', 'MarkerSize', 3);
end
xlim([-2 7]); ylim([1e-13 1e4]); grid on; xlabel('pF'); ylabel('K (cm d^{-1})');
legend('K_{nc}, Eq. (20)', 'K_{nc}, M_{19}', 'K_{nc}, M_{ex} (same K_{s,nc})', 'K_v', ...
    'K total, Eq. (20)', 'K total, M_{19}', 'measured K', 'Location', 'southwest', 'FontSize', 7);
title(sprintf('(d) %s: components toward oven dryness', strrep(ids(i), "_", "-")));
exportgraphics(f14, fullfile(outDir, 'FigS5_noncapillary_conductivity.png'), 'Resolution', 300); close(f14);
end
