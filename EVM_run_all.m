function EVM_run_all(dataDir, nSoilsMax)
%EVM_RUN_ALL  Runs Scripts S1-S8 in order and lists what each one produced.
%
%   EVM_run_all              % data files next to this file (about 56 min with 24 workers)
%   EVM_run_all(dataDir)     % data files in dataDir
%   EVM_run_all(dataDir, 20) % quick test with the first 20 German soils in every script
%
% Data: RetMeas.csv and CondMeas.csv (Hohenbrink et al., 2023) and, for Scripts S6 and S8,
% Table1_literature_RetMeas.csv and Table1_literature_CondMeas.csv. Script S6 compares its fits
% with the outputs of Scripts S2 and S4, Script S7 with those of Script S6, and Script S8 with
% those of Scripts S1, S6 and S7, so they run last. See README_reproducibility.txt.

if nargin < 1 || isempty(dataDir), dataDir = fileparts(mfilename('fullpath')); end
if nargin < 2, nSoilsMax = []; end
steps = {'S1', @EVM_beta_closure_analysis,           'beta_closure_results', 'Eq. (6) smoothness, soil-specific beta'; ...
         'S2', @EVM_positivity_constraint_analysis,  'positivity_results',   'positivity and theta_s > theta_r >= 0'; ...
         'S3', @EVM_suction_range_validation,        'range_results',        'validation by suction range'; ...
         'S4', @EVM_VG_comparison,                   'vg_results',           'comparison with van Genuchten'; ...
         'S5', @EVM_noncapillary_conductivity_check, 'ncond_results',        'Eq. (20) and Appendix B'; ...
         'S6', @EVM_reproduce_manuscript,            'manuscript_results',   'Tables 2-3, Fig. 1 (revised fits)'; ...
         'S7', @EVM_original_code_recheck,           'recheck_results',      'earlier procedure: Table 6, Figs. 2-3'; ...
         'S8', @EVM_beta_closure_recheck,            'beta_recheck_results', 'closure beta: sensitivity to c'};
nS = size(steps, 1); tMin = nan(nS, 1); ok = false(nS, 1); msg = strings(nS, 1);
for k = 1:nS
    fprintf('\n===== Script %s: %s =====\n', steps{k, 1}, func2str(steps{k, 2}));
    t0 = tic;
    try
        steps{k, 2}(dataDir, nSoilsMax);
        ok(k) = true;
    catch e
        msg(k) = string(e.message);
        fprintf(2, 'Script %s failed: %s\n', steps{k, 1}, e.message);
    end
    tMin(k) = toc(t0)/60;
end
fprintf('\n%-3s %-38s %-22s %-42s %6s\n', '', 'function', 'output folder', 'content', 'min');
for k = 1:nS
    st = ''; if ~ok(k), st = ['  FAILED: ' char(msg(k))]; end
    fprintf('%-3s %-38s %-22s %-42s %6.1f%s\n', steps{k, 1}, func2str(steps{k, 2}), steps{k, 3}, steps{k, 4}, tMin(k), st);
end
fprintf('EVM_RUN_ALL_DONE (%d of %d scripts completed, %.0f min)\n', sum(ok), nS, sum(tMin));
end
