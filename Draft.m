%% Proxy SVAR for US fiscal policy (Mertens & Ravn, 2014) - MAIN SCRIPT

% Required folder structure:
%   Draft.m
%   CK_RESTUD_DATASET.xlsx
%   additional files/
%       doPVAR.m          estimation + narrative (proxy) identification
%       doPVARbs.m        wild bootstrap
%       doPVARbs_MBB.m    moving block bootstrap
%       testResiduals.m   residual diagnostics

clear; close all; clc;
clear doPVAR   % resets the persistent flag so first-stage diagnostics print on every run

%% Paths
rootDir = fileparts(mfilename('fullpath'));
if isempty(rootDir)           
    rootDir = pwd;
end
addpath(fullfile(rootDir, 'additional files'));
dataFile = fullfile(rootDir, 'CK_RESTUD_DATASET.xlsx');

%% Settings
rng(1234, 'twister');   % fixed seed: bootstrap bands are identical on every run
nboot  = 1000;          % bootstrap replications (the paper uses 10000)
clevel = 95;            % confidence level (%)
L      = 8;             % block length for the moving block bootstrap (quarters)

%% Data (sheet FINAL, 1950Q1-2006Q4)
% DATES        : year.quarter (.00 = Q1, .25 = Q2, ...)
% TAX, G, GDP  : federal tax revenues, government spending, output
%                (logs, real, per capita)
% TAXNARRATIVE : Mertens & Ravn (2012) unanticipated narrative tax shocks
D = readtable(dataFile, 'Sheet', 'FINAL');

dates        = D.DATES;
VARNA.vars   = [D.TAX D.G D.GDP];
VARNA.p      = 4;
T            = size(VARNA.vars, 1);
VARNA.DET    = [ones(T,1) (1:T)' ((1:T).^2)' dates==1975.25];  % const, trend, trend^2, 1975Q2 dummy
VARNA.irhor  = 20;
VARNA.taxshocks = D.TAXNARRATIVE;

% Demean the narrative shocks (non-zero observations only)
for j = 1:size(VARNA.taxshocks, 2)
    nz = VARNA.taxshocks(:,j) ~= 0;
    VARNA.taxshocks(nz,j) = VARNA.taxshocks(nz,j) - mean(VARNA.taxshocks(nz,j));
end

% Shock size: 1% of GDP
DAT.TRY = 0.1822;   % average federal tax revenues / GDP
DAT.GY  = 0.2048;   % average federal expenditures / GDP
VARNA.tshocksize = 0.01/DAT.TRY;
VARNA.gshocksize = 0.01/DAT.GY;

%% Estimation and identification (point estimate)
VARNA = doPVAR(VARNA);

fprintf('Parameter estimates\n');
fprintf('  thetaG = %f\n', VARNA.thetaG);
fprintf('  thetaY = %f\n', VARNA.thetaY);
fprintf('  gammaT = %f\n', VARNA.gammaT);
fprintf('  zetaT  = %f\n', VARNA.zetaT);
fprintf('  zetaG  = %f\n', VARNA.zetaG);

%% Residual diagnostics
fprintf('\n========================================\n');
fprintf('VAR RESIDUAL DIAGNOSTICS\n');
fprintf('========================================\n');
residualTests = testResiduals(VARNA);

%% Bootstrap confidence bands
VARNAbs     = doPVARbs(VARNA, DAT, nboot, clevel);         % wild bootstrap
VARNAbs_MBB = doPVARbs_MBB(VARNA, DAT, nboot, clevel, L);  % moving block bootstrap

fprintf('\nMBB discards: allzero=%d, weakNZ=%d, failPVAR=%d, smallD11=%d\n', ...
    VARNAbs_MBB.discard_allzero, VARNAbs_MBB.discard_weakNZ, ...
    VARNAbs_MBB.discard_failPVAR, VARNAbs_MBB.discard_smallD11);

%% Impulse responses
% Note: quantile fields are named irsH = lower bound, irsL = upper bound
% (naming inherited from the original code).
h        = (0:VARNA.irhor-1)';                     % horizon 0 = impact
varNames = {'Tax Revenues', 'Government Spending', 'Output'};
yLims    = {[-10 10], [-10 10], [-2 5.5]};         % by variable
panelVar = [3 2 1];                                % plot order: Output, G, Tax
col      = [0 0 0.5];

shocks(1).title = 'Tax cut of 1% of GDP';
shocks(1).irf   = VARNA.irs;
shocks(1).wLo   = VARNAbs.irsH;      shocks(1).wHi = VARNAbs.irsL;
shocks(1).mLo   = VARNAbs_MBB.irsH;  shocks(1).mHi = VARNAbs_MBB.irsL;

shocks(2).title = 'Spending increase of 1% of GDP';
shocks(2).irf   = VARNA.irsg;
shocks(2).wLo   = VARNAbs.irsgH;     shocks(2).wHi = VARNAbs.irsgL;
shocks(2).mLo   = VARNAbs_MBB.irsgH; shocks(2).mHi = VARNAbs_MBB.irsgL;

for s = 1:numel(shocks)
    figure('Name', shocks(s).title, 'Position', [100 100 1200 380]);
    for k = 1:numel(panelVar)
        v = panelVar(k);
        subplot(1, 3, k); hold on; box on;
        p1 = plot(h, shocks(s).irf(:,v), '-', 'LineWidth', 2, 'Color', col);
        p2 = plot(h, [shocks(s).wLo(:,v) shocks(s).wHi(:,v)], '--', 'LineWidth', 1, 'Color', col);
        p3 = plot(h, [shocks(s).mLo(:,v) shocks(s).mHi(:,v)], ':', 'LineWidth', 1.5, 'Color', col);
        yline(0, 'k-');
        xlim([0 VARNA.irhor-1]); ylim(yLims{v});
        title(varNames{v}, 'FontSize', 14);
        xlabel('quarters', 'FontSize', 12); ylabel('percent', 'FontSize', 12);
        if k == 1
            legend([p1 p2(1) p3(1)], {'Point estimate', ...
                sprintf('Wild bootstrap %d%%', clevel), ...
                sprintf('MBB %d%%', clevel)}, 'Location', 'best');
        end
    end
    sgtitle(shocks(s).title, 'FontSize', 16);
end
