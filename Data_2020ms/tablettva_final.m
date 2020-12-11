% PROJECT:      WP1b - assessment of visual attention on a tablet device 
% AUTHOR:       Lulu Wang
% INSTITUTION:  KU Leuven
% CONTENT:      Statistical analysis and visualisation for manuscript
%               submission (and thesis)
%   1. General info:   Demographics, Performance, compare tablet-desktop
%   2. Fitting info:   Correlations observed-predicted
%   3. Parameter info: Show table of average, range, etc
%                      Check variations across datasets with 2 way ANOVA
%                      Compute correlations for Parallel-version reliability
%                      and Test-retest reliability
% Input: tablettva_raw.txt and tablettva_fit_awCu_fix.txt
% -------------------------------------------------------------------------
% Last edited: 2020.08.07

% -------------------------------------------------------------------------
%% Initialisation
% -------------------------------------------------------------------------
clc; clearvars; close all;

% Set paths
fn_dataraw = '/Users/Lulu/Documents/Experiments/WP1b_TVA_tablet_new/Data/tablettva_rawdata.txt';
fn_datafit = '/Users/Lulu/Documents/Experiments/WP1b_TVA_tablet_new/Code/2020_manuscript/tablettva_fit_awCu_fixed.txt';
fn_datafitsplit = {'/Users/Lulu/Documents/Experiments/WP1b_TVA_tablet_new/Code/2020_manuscript/tablettva_fit_awCu_A.txt'
    '/Users/Lulu/Documents/Experiments/WP1b_TVA_tablet_new/Code/2020_manuscript/tablettva_fit_awCu_B.txt'};
addpath('/Users/Lulu/Documents/MATLAB/simple_mixed_anova')
addpath('/Users/Lulu/Documents/Experiments/WP1b_TVA_tablet_new/Code/Functions')

% Read and process raw data
rawdata = readtable(fn_dataraw);
rawdata = rawdata(rawdata.ispractice==0,:); % remove practice trials
rawdata.sesmark = findgroups(rawdata.session,rawdata.device);
pinfo = unique(rawdata(:,1:7),'rows'); % without sep sessions
np = length(unique(rawdata.p_num));
[ntrialsp,~]=hist(rawdata.p_num,pinfo.p_num);
pinfo.ntrials = ntrialsp';

% Read and process fitted data
fitdata = readtable(fn_datafit); 
fitdata.p_num = cellfun(@(x) str2double(x(2:3)),fitdata.ID);
fitdata.w_index = fitdata.w2./(fitdata.w2+fitdata.w1); % Left/(left+Right)
fitdata.sesmark = findgroups(cellfun(@(x) x(end-2:end),fitdata.ID,'UniformOutput',false));
fitdata = fitdata(1:236,:);
fitsplit = cellfun(@readtable,fn_datafitsplit,'UniformOutput',false);
fitsplit = cellfun(@(x) {x(1:236,:)},fitsplit);

% Set groups and visualisation settings
vis.parnames = strsplit('K t0 C alpha w_index');
vis.parunits = strsplit('Elements ms Elements/s [-] [-]');
vis.parcols = cellfun(@(x) find(strcmp(fitdata.Properties.VariableNames,x)),vis.parnames);
vis.sesmark = {1,'D';1,'T';2,'D';2,'T'};
vis.expdur = [17 33 67 100 150 200];
vis.predcols = 245:268;
vis.obscols = 221:244;

% rmANOVA settings
factors = cell2table(vis.sesmark,'VariableNames',strsplit('ses dev'));
factors.ses = categorical(factors.ses); factors.dev = categorical(factors.dev);

fprintf('Ready!\n')

% -------------------------------------------------------------------------
%% General information
% -------------------------------------------------------------------------

% -- Materials and Methods - Participants
pincl = pinfo(pinfo.ntrials==4*324,:);
[n,~]=hist(findgroups(pincl.gender),1:2);
fprintf('# Participants:%d\nAge: %.1f ± %.1f, range %d-%d\n%d female\n',...
    np,mean(pincl.age),std(pincl.age),min(pincl.age),max(pincl.age),n(1))

% -------------------------------------------------------------------------
%% Results - General performance
% -------------------------------------------------------------------------
clc
% -- Get performance
if exist('perf','var')==0
perf = zeros(59,4);
for p = 1:59
    for s = 1:4
        data = rawdata(and(rawdata.p_num==p,rawdata.sesmark==s),:);
        n_responses = sum(cellfun(@length,data.response))-sum(cellfun(@(x) contains(x,'-'),data.response));
        n_correct = sum(cellfun(@(x,t) length(intersect(x,t)),data.response,data.targets));
        perf(p,s) = 100*n_correct/n_responses;
    end
end
end
fprintf('Mean performance: %.1f ± %.1f\n',mean(perf(:)),std(perf(:)))
fprintf('%.1f ± %.1f\n',[mean(perf); std(perf)])
fprintf('S1: %.1f ± %.1f\tS2:%.1f ± %.1f\n',[mean(reshape(perf(:,1:2),59*2,1))...
    std(reshape(perf(:,1:2),59*2,1)) mean(reshape(perf(:,3:4),59*2,1)) std(reshape(perf(:,3:4),59*2,1))])
fprintf('D: %.1f ± %.1f\tT:%.1f ± %.1f\n',[mean(reshape(perf(:,[1 3]),59*2,1))...
    std(reshape(perf(:,[1 3]),59*2,1)) mean(reshape(perf(:,[2 4]),59*2,1)) std(reshape(perf(:,[2 4]),59*2,1))])

% 2-way rmANOVA on performance
rmtab = array2table(perf); 
[ranovatbl,rm] = calculate_rmANOVA(rmtab,factors);
disp(ranovatbl)


%% -- Visualisation of raw performance with TVA parameters

% Select and process data
p = 6; % selected participant for manuscript: 6. strange alphas: 11 13 16 18 23 33 38 48 57
data = rawdata(and(rawdata.p_num==p,rawdata.sesmark<3),:); % select session 1
perf_p = [data.trial_code, data.sesmark, ... % trial conditions, device, total/left/right
    cellfun(@(x,t) length(intersect(x,t)),data.response,data.targets)...
    cellfun(@(x,t) length(intersect(x,t(4:6))),data.response,data.targets)...
    cellfun(@(x,t) length(intersect(x,t(1:3))),data.response,data.targets)]; 

% Visualise data
close all; font = 'Arial';
figure('pos',[63   480   935   272],'DefaultTextFontName', font, 'DefaultAxesFontName', font)
% Whole report
subplot(1,5,[1:3]); hold on;
s=1;tmp = splitapply(@mean,perf_p(perf_p(:,2)==s,3),perf_p(perf_p(:,2)==s,1));
scatter(vis.expdur,tmp(1:6),'kd','filled')
x = [fitdata.t0(and(fitdata.p_num==p,fitdata.sesmark==s)) vis.expdur];
v = [0 fitdata{and(fitdata.p_num==p,fitdata.sesmark==s),vis.predcols(1:6)}];
xq = 0:1:210; vq = interp1(x,v,xq,'pchip','extrap'); plot(xq(11:end),vq(11:end),'k-')

s=2;tmp = splitapply(@mean,perf_p(perf_p(:,2)==s,3),perf_p(perf_p(:,2)==s,1));
scatter(vis.expdur,tmp(1:6),'k^','filled')
x = [fitdata.t0(and(fitdata.p_num==p,fitdata.sesmark==s)) vis.expdur];
v = [0 fitdata{and(fitdata.p_num==p,fitdata.sesmark==s),vis.predcols(1:6)}];
vq = interp1(x,v,xq,'pchip','extrap'); plot(xq,vq,'k--')
legend('Desktop (data)','Desktop (model)','Tablet (data)','Tablet (model)','location','southeast');legend boxoff
xlabel('Exposure duration (ms)'); ylabel('Correctly reported letters')
% Add parameters to figure
presentresults = [fitdata(and(fitdata.p_num==p,fitdata.sesmark==1),vis.parcols);...
    fitdata(and(fitdata.p_num==p,fitdata.sesmark==2),vis.parcols)];
C_x = [presentresults.t0(1) 60]; C_y = [0 diff(C_x)/1000*presentresults.C(1)];
K_x = [120 210]; K_y = presentresults.K([1 1]);
plot(C_x,C_y,'-','Color',ones(1,3)*.5)
plot(K_x,K_y,'-','Color',ones(1,3)*.5)
text(presentresults.t0(1)-3,0-.1,'\it{t_{0}}','HorizontalAlignment', 'right','VerticalAlignment','top')
text(C_x(2),C_y(2),'\it{C}','HorizontalAlignment', 'left','VerticalAlignment','bottom')
text(K_x(1),K_y(1),'\it{K}','HorizontalAlignment', 'right','VerticalAlignment','middle')
text(min(xlim)-5,max(ylim),'A','FontWeight','bold','FontSize',14,'HorizontalAlignment', 'right','VerticalAlignment','bottom')

xlim([0 210]);yticks([0:4]);ylim([0 max(ylim)]); allylim = ylim;

% Partial report
subplot(1,5,4); hold on
s=1;tmp = splitapply(@mean,perf_p(perf_p(:,2)==s,3),perf_p(perf_p(:,2)==s,1));% 2T: 7:15, 2T4D: 16:24
perfpart = [mean(tmp(7:15)) mean(tmp(16:24))];
s=2;tmp = splitapply(@mean,perf_p(perf_p(:,2)==s,3),perf_p(perf_p(:,2)==s,1));
perfpart = [perfpart; [mean(tmp(7:15)) mean(tmp(16:24))]];
b = bar(perfpart');
b(1).FaceColor = 'k';b(2).FaceColor = ones(1,3)*.9; xlim([.5 2.5])
legend('Desktop','Tablet');xticks([1 2]);legend boxoff;xticklabels({'2T','2T4D'});xlabel('Condition');
yticks([0:4]);ylim(allylim)
presentresults = [presentresults array2table([perfpart perfpart(:,1)-perfpart(:,2)],...
    'VariableNames', strsplit('raw2T raw2T4D r2t4Ddiv2T'))];

% Hemifield
subplot(1,5,5); hold on
perfpart = [splitapply(@mean,perf_p(:,4),perf_p(:,2)) splitapply(@mean,perf_p(:,5),perf_p(:,2))];
b = bar(perfpart');
b(1).FaceColor = 'k';b(2).FaceColor = ones(1,3)*.9; xlim([.5 2.5])
legend('Desktop','Tablet');xticks([1 2]);legend boxoff;xticklabels({'Left','Right'});xlabel('Hemifield');
yticks([0:4]);ylim(allylim)

% Show parameters
% presentresults = [presentresults array2table([perfpart perfpart(:,1)./(perfpart(:,1)+perfpart(:,2))],...
%     'VariableNames', strsplit('rawL rawR rawLdivSum'))];
presentresults.Properties.RowNames = strsplit('D T');
fprintf('Participant %d presentresults:\n',p);disp(presentresults);
%% Check alphas
close all;
alphaDvsT = fitdata.alpha(fitdata.sesmark==1)-fitdata.alpha(fitdata.sesmark==2);
tdDvsT = zeros(size(alphaDvsT)); perfvariance = zeros(59,4); predms = zeros(59,2);
for p = 1:59
    data = rawdata(and(rawdata.p_num==p,rawdata.sesmark<3),:);
    perf_p = [data.trial_code, data.sesmark, ... % trial conditions, device, total/left/right
    cellfun(@(x,t) length(intersect(x,t)),data.response,data.targets)...
    cellfun(@(x,t) length(intersect(x,t(4:6))),data.response,data.targets)...
    cellfun(@(x,t) length(intersect(x,t(1:3))),data.response,data.targets)]; 
s=1;tmp = splitapply(@mean,perf_p(perf_p(:,2)==s,3),perf_p(perf_p(:,2)==s,1));% 2T: 7:15, 2T4D: 16:24
perfpart = mean(tmp(7:15))-mean(tmp(16:24));
s=2;tmp = splitapply(@mean,perf_p(perf_p(:,2)==s,3),perf_p(perf_p(:,2)==s,1));
tdDvsT(p) = perfpart - (mean(tmp(7:15))-mean(tmp(16:24)));
s=1;tmp = splitapply(@std,perf_p(perf_p(:,2)==s,3),perf_p(perf_p(:,2)==s,1)); perfvariance(p,1:2) = [mean(tmp(7:15)) mean(tmp(16:24))];
predms(p,s) = mean(fitdata{and(fitdata.p_num==p,fitdata.sesmark==s),251:259})-mean(fitdata{and(fitdata.p_num==p,fitdata.sesmark==s),260:268});
s=2;tmp = splitapply(@std,perf_p(perf_p(:,2)==s,3),perf_p(perf_p(:,2)==s,1)); perfvariance(p,3:4) = [mean(tmp(7:15)) mean(tmp(16:24))];
predms(p,s) = mean(fitdata{and(fitdata.p_num==p,fitdata.sesmark==s),251:259})-mean(fitdata{and(fitdata.p_num==p,fitdata.sesmark==s),260:268});

end
figure;scatter(alphaDvsT,tdDvsT);[R,P] = corrcoef(alphaDvsT,tdDvsT);title(sprintf('R=%.3f,p=%.4f',R(1,2),P(1,2)))
hold on; plot([0 0],ylim);plot(xlim,[0 0])
[alphaDvsT tdDvsT predms(:,1)-predms(:,2)]
% -------------------------------------------------------------------------
%% Results - TVA parameters (ANOVA)
% -------------------------------------------------------------------------
clc;

% Check parameter values
data = fitdata{[1:59*4],[vis.parcols end]};
meandata = splitapply(@mean,data(:,1:5),data(:,end))';
stddata = splitapply(@std,data(:,1:5),data(:,end))';
fprintf('Mean (SD) parameter values\nS1D\ts1T\ts2D\ts2T\n')
fprintf('%.3f (%.3f)\t%.3f (%.3f)\t%.3f (%.3f)\t%.3f (%.3f)\n',reshape([meandata;stddata],5,8)')

% Check fit performance
perf_fit = zeros(59,4,2); % participants, versions, R and P
for p = 1:59
    for s = 1:4
        [R,P] = corrcoef(fitdata{and(fitdata.p_num==p,fitdata.sesmark==s),vis.obscols(1:6)},...
            fitdata{and(fitdata.p_num==p,fitdata.sesmark==s),vis.predcols(1:6)});
        perf_fit(p,s,1) = R(1,2);perf_fit(p,s,2) = P(1,2);
    end
end
fprintf('\nFitting performance:%.3f ± %.3f\np-value < %.3f\n',...
    mean(reshape(perf_fit(:,:,1),59*4,1)),std(reshape(perf_fit(:,:,1),59*4,1)),max(reshape(perf_fit(:,:,2),59*4,1)))

% 2-way rmANOVA on fitting performance
rmtab = array2table(perf_fit(:,:,1)); 
[ranovatbl,rm] = calculate_rmANOVA(rmtab,factors);


% 2-way rmANOVA on parameters
for i = 1:5
    fprintf('\n%s\n',vis.parnames{i})
    data = splitapply(@(x) {x},fitdata(1:59*4,vis.parcols(i)),fitdata.sesmark(1:59*4));
    data = horzcat(data{:});
    rmtab = array2table(data); 
    [ranovatbl,rm] = calculate_rmANOVA(rmtab,factors);
    subplot(2,3,i)
    scatter(rmtab.data1,rmtab.data3,'b^');hold on
    scatter(rmtab.data2,rmtab.data4,'ro');
end


% -------------------------------------------------------------------------
%% Results - TVA parameters (Correlations)
% Steiger's Z: https://psych.unl.edu/psycrs/statpage/biv_corr_comp_eg.pdf
% -------------------------------------------------------------------------
clc;
parallel = zeros(5,8);
testretest = zeros(5,10);
internal = zeros(5,12);
SEZD = sqrt(1/(59 -3) + 1/(59-3));

for i = 1:5
    data = splitapply(@(x) {x},fitdata(1:59*4,vis.parcols(i)),fitdata.sesmark(1:59*4));
    data = horzcat(data{:});
    % Parallel version
    [R,P,RLO,RUP]=corrcoef(data(:,1),data(:,2));
    parallel(i,1:4) = [R(1,2),P(1,2),RLO(1,2),RUP(1,2)];
    [R,P,RLO,RUP]=corrcoef(data(:,3),data(:,4));
    parallel(i,5:8) = [R(1,2),P(1,2),RLO(1,2),RUP(1,2)];    
    % Test-retest
    [R,P,RLO,RUP]=corrcoef(data(:,1),data(:,3));
    testretest(i,1:4) = [R(1,2),P(1,2),RLO(1,2),RUP(1,2)];
    z1 = atanh(R(1,2));
    [R,P,RLO,RUP]=corrcoef(data(:,2),data(:,4));
    testretest(i,5:8) = [R(1,2),P(1,2),RLO(1,2),RUP(1,2)];  
    z2 = atanh(R(1,2));
    % Steiger's Z
    testretest(i,9) = (z1-z2)/SEZD; 
    testretest(i,10) = 2*normcdf(-abs(testretest(i,9)));%0.5 * erfc(-abs(testretest(i,9)) ./ sqrt(2));%
    %normcdf(abs(testretest(i,9)));
    % 
    
    % Internal reliability
    intcols = [2 10 13 15 19];
    for s = 1:4
        [R,P] = corrcoef(fitsplit{1}{contains(fitsplit{1}.ID,sprintf('_%d_',s)),intcols(i)},...
            fitsplit{2}{contains(fitsplit{2}.ID,sprintf('_%d_',s)),intcols(i)});
        rcorr = 2*R(1,2)/(1+R(1,2));
        internal(i,[1:2]+(s-1)*2) = [R(1,2) rcorr];
    end
    internal(i,9) = (atanh(internal(i,1))-atanh(internal(i,3)))/SEZD;
    internal(i,10) = 2*normcdf(-abs(internal(i,9)));
    internal(i,11) = (atanh(internal(i,5))-atanh(internal(i,7)))/SEZD;
    internal(i,12) = 2*normcdf(-abs(internal(i,11)));
end

fprintf('\nParallel-version reliability\n')
fprintf('%.3f(%.4f)\t%.3f, %.3f\t%.3f(%.4f)\t%.3f, %.3f\n',round(parallel,2)')
fprintf('\nTest-retest reliability\n')
fprintf('%.3f(%.4f)\t%.3f, %.3f\t%.3f(%.4f)\t%.3f, %.3f\t%.3f(%.3f)\n',round(testretest,2)')
fprintf('\nInternal reliability\n')
fprintf('%.3f(%.3f)\t%.3f(%.3f)\t%.3f(%.4f)\t%.3f(%.3f)\t%.3f(%.3f)\t%.3f(%.4f)\n',round(internal(:,[1:4 9:10 5:8 11:12])',2))

% -------------------------------------------------------------------------
%% Supplementary Results - Performance in the two groups
% -------------------------------------------------------------------------
%clc
for t = 3:4
    datamat = zeros(59,2);
    for p = 1:59
        data = rawdata(and(rawdata.p_num==p,rawdata.trial_code==t),:);
        % Number of correctly reported letters, not as a percentage of the
        % total reports
%         idx = ismember(data.sesmark,[1 3]);
%         datamat(p,1) = mean(cellfun(@(x,t) length(intersect(x,t)),data.response(idx),data.targets(idx)));
%         idx = ismember(data.sesmark,[2 4]);
%         datamat(p,2) = mean(cellfun(@(x,t) length(intersect(x,t)),data.response(idx),data.targets(idx)));

        % percentage of total reports
        idx = ismember(data.sesmark,[1 3]);
        datamat(p,1) = 100*sum(cellfun(@(x,t) length(intersect(x,t)),data.response(idx),data.targets(idx)))/...
            sum(cellfun(@(x) length(x)-length(intersect(x,'-')),data.response(idx)));
        idx = ismember(data.sesmark,[2 4]);
        datamat(p,2) = 100*sum(cellfun(@(x,t) length(intersect(x,t)),data.response(idx),data.targets(idx)))/...
            sum(cellfun(@(x) length(x)-length(intersect(x,'-')),data.response(idx)));       
    end
    fprintf('\nT = %d\n',t-2)
    between_factors = [zeros(24,1);ones(35,1)];
    [tbl,rm] = simple_mixed_anova(datamat, between_factors, {'ver'}, {'group'});disp(tbl)
end

fprintf('3-way mixed ANOVA\n')
datamat = zeros(59,2,2);
for t = 3:4
    for p = 1:59
        data = rawdata(and(rawdata.p_num==p,rawdata.trial_code==t),:);
        % Number of correctly reported letters, not as a percentage of the total reports
%         idx = ismember(data.sesmark,[1 3]);
%         datamat(p,1,t-2) = mean(cellfun(@(x,t) length(intersect(x,t)),data.response(idx),data.targets(idx)));
%         idx = ismember(data.sesmark,[2 4]);
%         datamat(p,2,t-2) = mean(cellfun(@(x,t) length(intersect(x,t)),data.response(idx),data.targets(idx)));
        
        % percentage of total reports
        idx = ismember(data.sesmark,[1 3]);
        datamat(p,1,t-2) = 100*sum(cellfun(@(x,t) length(intersect(x,t)),data.response(idx),data.targets(idx)))/...
            sum(cellfun(@(x) length(x)-length(intersect(x,'-')),data.response(idx)));
        idx = ismember(data.sesmark,[2 4]);
        datamat(p,2,t-2) = 100*sum(cellfun(@(x,t) length(intersect(x,t)),data.response(idx),data.targets(idx)))/...
            sum(cellfun(@(x) length(x)-length(intersect(x,'-')),data.response(idx)));       
    end
end
between_factors = [zeros(24,1);ones(35,1)];
[tbl,rm] = simple_mixed_anova(datamat, between_factors, {'ver','exp'}, {'group'});disp(tbl)

fprintf('2-way mixed ANOVA on fitting performance\n')
[tbl,rm] = simple_mixed_anova(perf_fit(:,1:2,1), between_factors, {'ver'}, {'group'});disp(tbl)






