% PROJECT:      WP1b - assessment of visual attention on a tablet device 
%               (But can be used for any rm ANOVA
% AUTHOR:       Lulu Wang
% INSTITUTION:  KU Leuven
% CONTENT:      Check data quality
%               1. Check trials, exposure durations, response performance
%               2. Plot figures
% Input:        var 1 - ttva_rawdata table format
%               var 2 - fitted data table format with added cols
%               var 3 - stage: raw, fit, all
% NF: requires at least 2 factors!
% -------------------------------------------------------------------------
function [ranovatbl,rm] = calculate_rmANOVA(rmtab,factors)
fprintf('************** calculate_rmANOVA **************\n')
rm = fitrm(rmtab,sprintf('%s-%s~1',rmtab.Properties.VariableNames{[1 end]}),'WithinDesign',factors);
[ranovatbl] = ranova(rm, 'WithinModel',...
    [sprintf('%s*',factors.Properties.VariableNames{1:end-1}) factors.Properties.VariableNames{end}]);

mau = mauchly(rm);
if and(mau.pValue<0.05,size(factors,2)>2)
    eps = epsilon(rm);
    disp('Mauchly test significant, calculate epsilon');disp([mau eps])
    ranovatbl.DF = ranovatbl.DF * eps.GreenhouseGeisser;
    ranovatbl = ranovatbl(:,[1:4 6:end-2]);
else
    ranovatbl = ranovatbl(:,[1:5 7:end-2]);
end
ranovatbl.eta2 = zeros(height(ranovatbl),1);
ranovatbl.eta2(1:2:end) = ...
    ranovatbl.SumSq(1:2:end)./(ranovatbl.SumSq(1:2:end)+ranovatbl.SumSq(2:2:end));
disp('Repeated Measures ANOVA:')
disp([ranovatbl(3:2:end,4:end) table([ranovatbl.DF(3:2:end) ranovatbl.DF(4:2:end)],...
    'VariableName',{'DF'})])

rmtabrows = 3:2:height(ranovatbl);
if size(factors,1)/size(factors,2)>2
if any(ranovatbl{rmtabrows,5}<0.05),disp('Multiple comparisons:'); 
for i = 1:length(rmtabrows)
if ranovatbl{rmtabrows(i),5}<0.05
    if i<length(factors.Properties.VariableNames)+1
        mc = multcompare(rm,factors.Properties.VariableNames{i}); 
    else
        tmp = strsplit(ranovatbl.Properties.RowNames{i},':');
        mc = multcompare(rm,tmp{2},'by',tmp{3}); 
    end
    disp(mc(mc.Difference>0,:)); 
end
end
end
end
fprintf('**************  **************  **************\n')
    
    
    
    