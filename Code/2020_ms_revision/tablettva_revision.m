%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TITLE:        Encouraging digital technology in neuropsychology: 
%               The Theory of Visual Attention on tablet devices
% AUTHOR:       Lulu Wang
% INSTITUTION:  KU Leuven
% PROJECT:      WP1b - assessment of visual attention on a tablet device 
% PURPOSE:      Analyze data and visualize results
% -------------------------------------------------------------------------
% Revision:     ver. 2020.10.21
% -------------------------------------------------------------------------
% FUNCTIONS:    visualiseTVA.m
%               calculate_rmANOVA.m
%               f_CCC.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =========================================================================
%% Initialisation
% =========================================================================
clc; clearvars; close all;

% -------------------------------------------------------------------------
% SET PATHS
% -------------------------------------------------------------------------

dirs.main = '/Users/Lulu/Documents/Experiments/WP1b_TVA_tablet_new';
cd([dirs.main '/Code/2020_ms_revision'])
addpath([dirs.main '/Code/Functions'])

% -------------------------------------------------------------------------
% GET DATA
% -------------------------------------------------------------------------
dirs.data   = [dirs.main '/Code/2020_manuscript/'];
fn.raw      = [dirs.data 'tablettva_rawdata.txt'];
fn.fit      = [dirs.data 'tablettva_fit_awCu_fixed.txt'];
fn.fitsplit = {[dirs.data 'tablettva_fit_awCu_A.txt']; [dirs.data 'tablettva_fit_awCu_B.txt']};
disp('Loading raw data...')
d.raw       = readtable(fn.raw);
disp('Loading fitted data...')
d.fit       = readtable(fn.fit);
disp('Loading split data...')
d.split     = cellfun(@readtable,fn.fitsplit,'UniformOutput',false);

% -------------------------------------------------------------------------
% PREPARE DATA/ANALYSIS SETTINGS
% -------------------------------------------------------------------------

% Add visualisation settings
vis.parnames = strsplit('K t0 C alpha w2');
vis.parunits = strsplit('Elements ms Elements/s [-] [-]');
vis.parcols  = cellfun(@(x) find(strcmp(d.fit.Properties.VariableNames,x)),vis.parnames,'UniformOutput',false);
vis.parcols  = [vis.parcols{:}];
vis.sesmark  = {1,'D';1,'T';2,'D';2,'T'};
vis.expdur   = [17 33 67 100 150 200];
vis.predcols = 245:268;
vis.obscols  = 221:244;

% Add rmANOVA settings
factors      = cell2table(vis.sesmark,'VariableNames',strsplit('ses dev'));
factors.ses  = categorical(factors.ses); factors.dev = categorical(factors.dev);

fprintf('Ready!\n\n')

% =========================================================================
%% Data pre-processing
% =========================================================================

% -------------------------------------------------------------------------
% CLEAN DATA
% -------------------------------------------------------------------------

% Add data markers
d.raw.sesmark = findgroups(d.raw.session,d.raw.device);
d.fit.p_num   = cellfun(@(x) str2double(x(2:3)),d.fit.ID);
d.fit.sesmark = cellfun(@(x) str2double(x(5)),d.fit.ID);
for s = 1:2
    d.split{s}.p_num   = cellfun(@(x) str2double(x(2:3)),d.split{s}.ID);
    d.split{s}.sesmark = cellfun(@(x) str2double(x(5)),d.split{s}.ID);
end

% Get participant information
d.info         = unique(d.raw(:,1:9),'rows'); 
d.info.sesmark = findgroups(d.info.session,d.info.device);
fprintf('Total participants: %d\n',length(unique(d.info.p_num)))

% Remove incomplete datasets
[n,~]   = hist(d.info.p_num,unique(d.info.p_num)); 
d.incl  = unique(d.info(:,1:7),'rows'); d.incl = d.incl(n>2,:);
d.split = cellfun(@(x) {x(1:236,:)},d.split); 
d.fit   = d.fit(1:236,:);
d.raw   = d.raw(ismember(d.raw.p_num,d.incl.p_num),:);

% Sample information
fprintf('\nCompleted datasets: %d\n',length(unique(d.incl.p_num)))
[n,~] = hist(findgroups(d.incl.gender),[1 2]); fprintf(' %d women, %d men\n',n(1),n(2));
fprintf(' age: %.1f±%.1f, range %d-%d\n',mean(d.incl.age),std(d.incl.age),min(d.incl.age),max(d.incl.age))


% -------------------------------------------------------------------------
% CHECK FOR OUTLIERS
% -------------------------------------------------------------------------

fprintf('\nParameter values (over all datasets):\n')
close all;
for i = 1:length(vis.parcols)
    tmp = d.fit{:,vis.parcols(i)};
    subplot(2,3,i); hist(tmp);title(vis.parnames{i})
    fprintf('%s:\t%.2f±%.2f, range %.2f-%.2f (%.2f/%.2f)\n',vis.parnames{i},mean(tmp),std(tmp),...
        min(tmp),max(tmp),(mean(tmp)-min(tmp))/std(tmp),(max(tmp)-mean(tmp))/std(tmp))
end

% -------------------------------------------------------------------------
% REMOVE ALPHA OUTLIER
% -------------------------------------------------------------------------

% close all;
% visualiseTVA(52,1,d,vis)

d.incl = d.incl([1:51 53:end],:);
fprintf('\nCompleted datasets: %d\n',length(unique(d.incl.p_num)))
[n,~] = hist(findgroups(d.incl.gender),[1 2]); fprintf(' %d women, %d men\n',n(1),n(2));
fprintf(' age: %.1f±%.1f, range %d-%d\n',mean(d.incl.age),std(d.incl.age),min(d.incl.age),max(d.incl.age))

% Remove incomplete datasets
d.fit = d.fit(ismember(d.fit.p_num,d.incl.p_num),:);
for s = 1:2, d.split{s} = d.split{s}(ismember(d.split{s}.p_num,d.incl.p_num),:);end
d.raw = d.raw(ismember(d.raw.p_num,d.incl.p_num),:);

% Check for outliers again
fprintf('\nParameter values (over all datasets):\n')

for i = 1:length(vis.parcols)
    tmp = d.fit{:,vis.parcols(i)};
%     subplot(2,3,i);hist(tmp);title(vis.parnames{i})
    fprintf('%s:\t%.2f±%.2f, range %.2f-%.2f (%.2f/%.2f)\n',vis.parnames{i},mean(tmp),std(tmp),...
        min(tmp),max(tmp),(mean(tmp)-min(tmp))/std(tmp),(max(tmp)-mean(tmp))/std(tmp))
    clear tmp;
end

N = height(d.incl);
fprintf('\nReady!\n\n')

% =========================================================================
%% Analyze raw response accuracy
% =========================================================================
close all;

if exist('percorrect','var')==0
percorrect = zeros(height(d.incl),4);
for p = 1:height(d.incl)
    for s = 1:4
        tmp = d.raw(and(d.raw.p_num == d.incl.p_num(p),and(d.raw.ispractice==0,d.raw.sesmark==s)),:);
        percorrect(p,s) = sum(cellfun(@(xt,xr) length(intersect(xt,xr)),tmp.targets,tmp.response))./...
            sum(cellfun(@(xr) length(xr)-length(intersect(xr,'-')),tmp.response))*100;
    end
end
end

fprintf('***\nMean performance across all participants and datasets: %.1f ± %.1f, range %.1f-%.1f\n',...
    mean(percorrect(:)),std(percorrect(:)),min(percorrect(:)),max(percorrect(:)))

tmp = reshape(percorrect(:,1:2),1,N*2); 
fprintf('\nSession 1:%.1f ± %.1f, %.1f-%.1f\n',mean(tmp), std(tmp),min(tmp),max(tmp))
tmp = reshape(percorrect(:,3:4),1,N*2); 
fprintf('Session 2:%.1f ± %.1f, %.1f-%.1f\n',mean(tmp), std(tmp),min(tmp),max(tmp))
tmp = reshape(percorrect(:,[1 3]),1,N*2); 
fprintf('Desktop:  %.1f ± %.1f, %.1f-%.1f\n',mean(tmp), std(tmp),min(tmp),max(tmp))
tmp = reshape(percorrect(:,[2 4]),1,N*2); 
fprintf('Tablet:   %.1f ± %.1f, %.1f-%.1f\n',mean(tmp), std(tmp),min(tmp),max(tmp))

fprintf('Per dataset:\n')
fprintf('%.1f ± %.1f, %.1f-%.1f\n',[mean(percorrect);std(percorrect);min(percorrect);max(percorrect)])

% 2-way rmANOVA
rmtab = array2table(percorrect); 
calculate_rmANOVA(rmtab,factors);

% =========================================================================
%% Analyze TVA-fit performance
% =========================================================================

if exist('tvafit','var')==0
tvafit = zeros(height(d.incl),4); tvafitp = zeros(height(d.incl),4);
for p = 1:height(d.incl)
    for s = 1:4
        tmp = d.fit(and(d.fit.p_num==d.incl.p_num(p),d.fit.sesmark==s),[221:226 245:250]);
        [R,P] = corrcoef(tmp{1,1:6},tmp{1,7:12});
        tvafit(p,s) = R(1,2);
        tvafitp(p,s) = P(1,2);
    end
end
end

fprintf('***\nMean fit across all participants and datasets: %.3f ± %.3f, range %.3f-%.3f, p < %.3f\n',...
    mean(tvafit(:)),std(tvafit(:)),min(tvafit(:)),max(tvafit(:)),max(tvafitp(:)))
fprintf('Per dataset:\n')
fprintf('%.3f ± %.3f, %.3f-%.3f\n',[mean(tvafit);std(tvafit);min(tvafit);max(tvafit)])

% 2-way rmANOVA
rmtab = array2table(tvafit); 
calculate_rmANOVA(rmtab,factors);

% =========================================================================
%% Analyze TVA-parameters 
% =========================================================================
clc; close all;
fprintf('***\nTVA parameter values\n\n')

% Get descriptive statistics
fprintf('\nName\ts1D\ts1T\ts2D\ts2T\n')
for i = 1:5
    fprintf('%s\t',vis.parnames{i})
    for j = 1:4
        tmp = d.fit{d.fit.sesmark==j,vis.parcols(i)};
        if i == 5,tmp = tmp./(1+tmp);end
        fprintf('%.2f (%.2f)\t',mean(tmp),std(tmp))
    end
    fprintf('\n')
end

% -------------------------------------------------------------------------
% RMANOVA SESSION/VERSION
% -------------------------------------------------------------------------

for i = 1:5
    fprintf('\n%s\n',vis.parnames{i})
    rmtab = array2table(reshape(d.fit{:,vis.parcols(i)},4,height(d.fit)/4)'); 
    [ranovatbl,rm] = calculate_rmANOVA(rmtab,factors);
end


% -------------------------------------------------------------------------
% RELIABILITY RETEST & PARALLEL
% -------------------------------------------------------------------------

fprintf('\n\n******************* TEST RELIABILITY ***********************\n')
SEZD     = sqrt(1/(N -3) + 1/(N-3));

retest   = zeros(5,2);
parallel = zeros(5,2);
tvatitles = {'VSTM capacity','perceptual threshold','processing speed',...
             'attentional selectivity','spatial attention bias'};

for istr = 0:1 % 1: test-retest; 0: parallel-version
    
if istr, cols = [1 3; 2 4]; fprintf('\n***\nTEST-RETEST\n')
else,    cols = [1 2; 3 4]; fprintf('\n***\nPARALLEL-VERSION\n')
end

% Create visualisation
figure('pos',[728   350   713   448])
    
for i = 1:length(vis.parcols)
    subplot(2,3,i); hold on
    disp(vis.parnames{i})
    
    % Calculate correlation coefficients
    clear s p R P RLO RUP; CCCrho = [0 0];
    for j = 1:2
        tmp = [d.fit{d.fit.sesmark==cols(j,1),vis.parcols(i)},d.fit{d.fit.sesmark==cols(j,2),vis.parcols(i)}];
        s{j} = scatter(tmp(:,1),tmp(:,2),'.');
        p{j} = plot([min(tmp(:,1)); max(tmp(:,1))],[min(tmp(:,2)); max(tmp(:,2))],'.-','color',s{j}.CData);
        [R{j},P{j},RLO{j},RUP{j}]=corrcoef(tmp(:,1),tmp(:,2));

        CCC = f_CCC(tmp,0.05);
        fprintf('CCC = %.2f (%.2f to %.2f)\tPearson R = %.2f (%.2f to %.2f)\n',...
            CCC{1}.est,CCC{1}.confInterval,R{j}(1,2),RLO{j}(1,2),RUP{j}(1,2))
        CCCrho(j) = CCC{1}.est;
    end
    
    % Calculate Steiger's Z
    if istr
        retest(i,1) = R{1}(1,2);retest(i,2) = R{2}(1,2);
        zPear    = cellfun(@(x) atanh(x(1,2)),R);
        zCCC     = atanh(CCCrho);
        steigerr = [(zCCC(1)-zCCC(2))/SEZD;   2*normcdf(-abs((zCCC(1)-zCCC(2))/SEZD));...
                    (zPear(1)-zPear(2))/SEZD; 2*normcdf(-abs((zPear(1)-zPear(2))/SEZD))];
        fprintf('Steigers CCC = %.2f (%.3f)\tPearson %.3f (%.3f)\n',steigerr)
    else
        parallel(i,1) = R{1}(1,2);parallel(i,2) = R{2}(1,2);
    end
    fprintf('\n')
    
    % Visualisation settings
    axis('square'); 
    xlim([min([xlim ylim]) max([xlim ylim])]);
    ylim([min([xlim ylim]) max([xlim ylim])])
    title(tvatitles{i})
    if istr
        xlabel('session 1');ylabel('session 2');
        l = legend([p{1},p{2}],'Desktop','Tablet','location','southeast');
    else
        xlabel('Desktop');ylabel('Tablet');
        l = legend([p{1},p{2}],'session 1','session 2','location','southeast');
    end
    legend boxoff;
%     text(diff(xlim)*.05,diff(ylim)*.95,sprintf('D: %.2f (%.3f)\nT: %.2f (%.3f)',R{1}(1,2),P{1}(1,2),...
%         R{2}(1,2),P{2}(1,2)), 'HorizontalAlignment','left','VerticalAlignment','top')
    plot(xlim,ylim,':','color',ones(1,3)*.5)
end
end


% -------------------------------------------------------------------------
% INTERNAL RELIABILITY
% -------------------------------------------------------------------------
fprintf('\n\nInternal reliability\n\nName\ts1D\ts1T\ts2D\ts2T\tSteiger1\tSteiger2\n')
internal  = zeros(5,4,4);
for i = 1:5
    for s = 1:4
        % Calculate correlation coefficients
        [R,P] = corrcoef(d.split{1}{d.split{1}.sesmark==s,vis.parcols(i)},...
                         d.split{2}{d.split{2}.sesmark==s,vis.parcols(i)});
        internal(i,s,1) = R(1,2);
        internal(i,s,2) = P(1,2);
        internal(i,s,3) = 2*R(1,2)/(1+R(1,2));
        % Steiger's Z-test
        if mod(s,2)==0 
            internal(i,s-1,4) = (atanh(internal(i,s,1))-atanh(internal(i,s-1,1)))/SEZD;
            internal(i,s,4)   = 2*normcdf(-abs(internal(i,s-1,4)));
        end
    end
    fprintf('%s\t %.2f (%.2f)\t%.2f (%.2f)\t%.2f (%.2f)\t%.2f (%.2f)\t%.2f (%.3f)\t%.2f (%.3f)\n',...
        vis.parnames{i},reshape([internal(i,:,1);internal(i,:,3);],8,1),internal(i,:,4))
end


% =========================================================================
%% Additional checks for the Reviewer
% =========================================================================
% 2020.10.27

% -------------------------------------------------------------------------
% Correction for Attenuation
% -------------------------------------------------------------------------
% Info:
% https://onlinelibrary-wiley-com.kuleuven.ezproxy.kuleuven.be/doi/full/10.1111/test.12087
% http://www.med.mcgill.ca/epidemiology/hanley/bios601/Surveys/AttenuationOfCorrelations-2016-Teaching_Statistics.pdf

% retest = [0.7543    0.7282;    0.6400    0.6759;    0.6379    0.7331;    0.7783    0.6723;    0.6055    0.7750];
% parallel = [0.7819    0.9272;    0.7594    0.7743;    0.7734    0.7685;    0.7501    0.6683;    0.6974    0.8742];

rt = parallel./sqrt(retest(:,1).*retest(:,2));
disp('Correct parallel-version reliability for attenuation')
disp(array2table([parallel rt],'VariableNames',strsplit('S1raw S2raw S1corrected S2corrected'),'RowNames',vis.parnames))


% -------------------------------------------------------------------------
% Time-on-task effects for sustained attention
% -------------------------------------------------------------------------

if exist('percorrect_block','var')==0
percorrect_block = zeros(height(d.incl),4,9);
for p = 1:height(d.incl)
    for s = 1:4
        for b = 1:9
        tmp = d.raw(and(and(d.raw.p_num == d.incl.p_num(p),d.raw.blocknr == b),...
            and(d.raw.ispractice==0,d.raw.sesmark==s)),:);
        percorrect_block(p,s,b) = sum(cellfun(@(xt,xr) length(intersect(xt,xr)),tmp.targets,tmp.response))./...
            sum(cellfun(@(xr) length(xr)-length(intersect(xr,'-')),tmp.response))*100;
        end
    end
end
end

% Visualise performance trend across runs
figure
for s = 1:4
    subplot(2,2,s)
    errorbar(1:9,mean(squeeze(percorrect_block(:,s,:))),std(squeeze(percorrect_block(:,s,:))));
    title(sprintf('s%d%s',vis.sesmark{s,:}))
end

% Compare performance in first 3 blocks versus final 3 blocks
disp('Check for time-on-task effects on response accuracy')
for s = 1:4
    fprintf('s%d%s\t',vis.sesmark{s,:})
    [H,P,CI,STATS] = ttest(mean(squeeze(percorrect_block(:,s,1:3)),2),mean(squeeze(percorrect_block(:,s,7:9)),2));
    fprintf('H=%d, p = %.3f,CI %.3f - %.3f, t(%d) = %.2f\n',H,P,CI(1),CI(2),STATS.df,STATS.tstat)
end

[H,P,CI,STATS] = ttest(mean(squeeze(mean(percorrect_block(:,:,1:3),2)),2),...
    mean(squeeze(mean(percorrect_block(:,:,7:9),2)),2));
fprintf('All\tH=%d, p = %.3f,CI %.3f - %.3f, t(%d) = %.2f\n',H,P,CI(1),CI(2),STATS.df,STATS.tstat)

fprintf('\nReady!\n')



