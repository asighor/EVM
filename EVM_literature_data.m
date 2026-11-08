function EVM_literature_data(srcDir, outDir, clayLoamFile)
%EVM_LITERATURE_DATA  Builds Data Set S2, the five literature soils of Table 1, in the format of
% the German database (RetMeas.csv and CondMeas.csv: Sample_ID, MeasType, pF, theta or k).
%
%   EVM_literature_data(srcDir)          % srcDir holds the downloaded data sets; output next to this file
%   EVM_literature_data(srcDir, outDir)
%   EVM_literature_data(srcDir, outDir, clayLoamFile)
%
% Sources: Data Sets S11, S12, S14 and S15 of Rudiyanto et al. (2015), A complete soil hydraulic model
% accounting for capillary and adsorptive water retention, capillary and film conductivity, and
% hysteresis, Water Resour. Res., 51, 8757-8772, doi:10.1002/2015WR017703. In each file, columns 1-2
% are suction h (cm) and theta, columns 3-4 suction h (cm) and K (cm/d), from row 4. The Pachepsky
% clay loam (Pachepsky et al., 1984) is not among these data sets; its data are read from the
% authors' data file, clayLoamFile (default: clay_loam_Pachepsky_1984.csv in srcDir), whose columns 1-2 are
% suction (m) and theta and columns 3-4 suction (m) and K (m/s). The KSAT rows (pF = -2) are the
% capillary Ks of Table 1 of the manuscript.
% Output: Table1_literature_RetMeas.csv and Table1_literature_CondMeas.csv (MeasType "Literature").

if nargin < 2 || isempty(outDir), outDir = fileparts(mfilename('fullpath')); end
if nargin < 3 || isempty(clayLoamFile), clayLoamFile = fullfile(srcDir, 'clay_loam_Pachepsky_1984.csv'); end
L = {'Berlin_Sand',  'wrcr21760-sup-0011-2015wr017703-ds11.xls', 673;     % Peters (2013)
     'Pachepsky_SL', 'wrcr21760-sup-0015-2015wr017703-ds15.xls', 12.8;    % Pachepsky et al. (1984)
     'Gilat_Loam',   'wrcr21760-sup-0012-2015wr017703-ds12.xls', 17.28;   % Mualem (1976b)
     'Pachapa_FSC',  'wrcr21760-sup-0014-2015wr017703-ds14.xls', 12.1};   % Mualem (1976b)
R = table(); K = table();
num = @(x) isnumeric(x) && isscalar(x) && isfinite(x);
for i = 1:size(L, 1)
    c = readcell(fullfile(srcDir, L{i, 2}));
    c = c(4:end, :);                                   % skip the title and the two header rows
    for r = 1:size(c, 1)
        if num(c{r, 1}) && num(c{r, 2})
            R = [R; table(string(L{i, 1}), "Literature", log10(c{r, 1}), c{r, 2}, ...
                'VariableNames', {'Sample_ID', 'MeasType', 'pF', 'theta'})]; %#ok<AGROW>
        end
        if size(c, 2) >= 4 && num(c{r, 3}) && num(c{r, 4})
            K = [K; table(string(L{i, 1}), "Literature", log10(c{r, 3}), c{r, 4}, ...
                'VariableNames', {'Sample_ID', 'MeasType', 'pF', 'k'})]; %#ok<AGROW>
        end
    end
    K = [K; table(string(L{i, 1}), "KSAT", -2, L{i, 3}, 'VariableNames', {'Sample_ID', 'MeasType', 'pF', 'k'})]; %#ok<AGROW>
end
% Pachepsky clay loam: suction m -> cm, K m/s -> cm/d; Ks of Table 1
c = readcell(clayLoamFile);
for r = 1:size(c, 1)
    if num(c{r, 1}) && num(c{r, 2})
        R = [R; table("Pachepsky_CL", "Literature", log10(100*c{r, 1}), c{r, 2}, ...
            'VariableNames', {'Sample_ID', 'MeasType', 'pF', 'theta'})]; %#ok<AGROW>
    end
    if size(c, 2) >= 4 && num(c{r, 3}) && num(c{r, 4})
        K = [K; table("Pachepsky_CL", "Literature", log10(100*c{r, 3}), 8.64e6*c{r, 4}, ...
            'VariableNames', {'Sample_ID', 'MeasType', 'pF', 'k'})]; %#ok<AGROW>
    end
end
K = [K; table("Pachepsky_CL", "KSAT", -2, 28.94, 'VariableNames', {'Sample_ID', 'MeasType', 'pF', 'k'})];
writetable(R, fullfile(outDir, 'Table1_literature_RetMeas.csv'));
writetable(K, fullfile(outDir, 'Table1_literature_CondMeas.csv'));
fprintf('%d retention and %d conductivity rows written to %s\n', height(R), height(K), outDir);
end
