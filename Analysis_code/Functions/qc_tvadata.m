% PROJECT:      WP1b - assessment of visual attention on a tablet device 
% AUTHOR:       Lulu Wang
% INSTITUTION:  KU Leuven
% CONTENT:      Check data quality
%               1. Check trials, exposure durations, response performance
%               2. Plot figures
% Input:        var 1 - ttva_rawdata table format
%               var 2 - fitted data table format with added cols
%               var 3 - stage: raw, fit, all
% -------------------------------------------------------------------------

function qc_tvadata(rawd,fitd,stage)
if nargin==2, stage = 'fit'; elseif nargin==1, stage = 'raw'; end

%% Initialisation

pnums = unique(rawd.p_num);
psesm = unique(rawd.sesmark);

% plot settings
sp = [length(psesm) 4];
session_markers = {1,'D';1,'T';2,'D';2,'T'};
expdur = [17 33 50 83 150 200];

%% 1. Check trials and exposure durations

for pp = 1:length(pnums)
    p = pnums(pp);
    fprintf('\n***************\nP%02d\n',p)
    figure%('pos',[1         378        1440         420])
    
    for ss = 1:length(psesm)
        s = psesm(ss);
        fprintf('----------------\nSes: %d\tDev: %s\n',session_markers{s,:})
        % prep data, remove practice trials
        data = rawd(and(rawd.p_num==p,rawd.sesmark == s),:);
        data = data(data.ispractice==0,:);
        if ~strcmp(stage,'raw'),fdata = fitd(and(fitd.p_num==p,fitd.p_sesmark==s),:);end
        
        if or(strcmp(stage,'raw'),strcmp(stage,'all'))
        % Check trials
        n_trials = height(data);
        [n,h]=hist(double(data.trial_code),unique(data.trial_code)');
        if isequal(n,[repmat(27,1,6) repmat(9,1,18)]),disp('Trial count OK');
        else, disp('ERROR in trial count!');disp('Trials per condition:'); disp([h;n]);end

        % Check exposure duration
        stim_dur = table(strsplit('1 2 3 4 5 6 7-24')',...
            [splitapply(@mean,data.stim_dur(data.trial_code<7),...
            data.trial_code(data.trial_code<7));mean(data.stim_dur(data.trial_code>6))],...
                   [splitapply(@std,data.stim_dur(data.trial_code<7),...
                   data.trial_code(data.trial_code<7));std(data.stim_dur(data.trial_code>6))],...
                   [splitapply(@min,data.stim_dur(data.trial_code<7),...
                   data.trial_code(data.trial_code<7));min(data.stim_dur(data.trial_code>6))],...
                   [splitapply(@max,data.stim_dur(data.trial_code<7),...
                   data.trial_code(data.trial_code<7));max(data.stim_dur(data.trial_code>6))],...
                   'VariableNames',strsplit('trialcode mean std min max'));
        if and(max(stim_dur.max-stim_dur.min)<5,max(stim_dur.std)<2),disp('Timing OK');
        else,disp('Timing off');disp(stim_dur);end
        
        % Check performance
        performance = [[1:length(unique(data.blocknr))]' zeros(length(unique(data.blocknr)),1)];
        for b = 1:length(unique(data.blocknr))
            tmp = data(data.blocknr==b,:);
            performance(b,2) = 100*sum(cellfun(@(t,r) length(intersect(t,r)),...
                tmp.targets,tmp.response))/sum(cellfun(@(r) length(r),tmp.response));
        end
        if any(abs(performance(:,2)-85)>15),disp('Performance low');
            disp(array2table(performance,'VariableNames',strsplit('block performance')))
        else, disp('Performance OK');end
        end
        
        % Check fit performance Whole report
        if or(strcmp(stage,'fit'),strcmp(stage,'all'))
            ms = fdata{:,find(strcmp(fdata.Properties.VariableNames,'MeanScoreC1'))+[0:5]};
            pms = fdata{:,find(strcmp(fdata.Properties.VariableNames,'PredMeanScoreC1'))+[0:5]};
            [R,P] = corrcoef(ms,pms);
            if P(1,2)<0.05
                fprintf('\nlibTVA fitting performance OK: r^2 = %.3f, p = %.4f\n',R(1,2),P(1,2))
            else
                fprintf('\nlibTVA fitting performance low: r^2 = %.3f, p = %.4f\n',R(1,2),P(1,2))
            end
            fprintf('K:  %.2f\nt0: %.1f ± %.1f\nC:  %.1f ± %.1f\na:  %.3f ± %.3f\nw:  %.3f ± %.3f\n',...
                fdata.K, fdata.t0, fdata.t0_err, fdata.C, fdata.C_err, fdata.alpha,...
                fdata.alpha_err,fdata.w_index, fdata.w2_err)
        end


        %% 2. Plot figures

        % Prepare left-right responses
        n_responses = cellfun(@length,data.response);
        n_corr = zeros(size(n_responses)); n_distr = zeros(size(n_responses)); 
        n_corrLR = zeros(length(n_responses),2);
        for t = 1:n_trials
            n_corr(t) = sum(cellfun(@(r) contains(data.targets{t},r),cellstr(data.response{t}')));
            n_distr(t) = sum(cellfun(@(r) contains(data.distractors{t},r),...
                cellstr(data.response{t}')));
            n_corrLR(t,:) = [sum(cell2mat(cellfun(@(r) strfind(data.targets{t},r),...
                cellstr(data.response{t}'),'UniformOutput',false))<4),...
                sum(cell2mat(cellfun(@(r) strfind(data.targets{t},r),...
                cellstr(data.response{t}'),'UniformOutput',false))>3)];

        end
        n_corrLR = n_corrLR(:,[2 1]);
    
        % Whole report responses
        subplot(sp(1),sp(2),[1 2]+(ss-1)*sp(2))
        if strcmp(stage,'raw')
            plot(expdur,splitapply(@mean, n_corr(data.trial_code<7), ...
                data.trial_code(data.trial_code<7)),'o-')
        else
            plot(expdur,ms,'o');hold on
            plot(expdur,pms,'-');
            legend(strsplit('observed predicted'))
        end
        xlim([0 210]); ylim([0 6]); 
        if s==max(psesm),xlabel('exposure duration (ms)');end
        ylabel('mean correct responses');
        title(sprintf('P%02d s%d%s Whole report 6T',p,session_markers{s,:}))

        subplot(sp(1),sp(2),3+(ss-1)*sp(2))
        stmp = [mean(n_corr(ismember(data.trial_code,7:15)))...
            mean(n_corr(ismember(data.trial_code,16:24)))];
        if strcmp(stage,'raw')
            mtmp = [mean(n_corr(ismember(data.trial_code,7:15)))...
                mean(n_corr(ismember(data.trial_code,16:24)))];
            bar(mtmp);hold on; errorbar(1:2,mtmp,stmp,'k.')
        else
            mtmp = [mean(n_corr(ismember(data.trial_code,7:15)))...
                mean(n_corr(ismember(data.trial_code,16:24))); ...
                mean(fdata{:,find(strcmp(fdata.Properties.VariableNames,'PredMeanScoreC7'))+[0:8]}) ...
                mean(fdata{:,find(strcmp(fdata.Properties.VariableNames,'PredMeanScoreC16'))+[0:8]})];

            bar(mtmp');hold on
            errorbar([1:2]-.15,mtmp(1,:),stmp,'k.')
        end
        ylim([0 max(ylim)])
        xticklabels({'2T','2T4D'});xlim([.5 2.5]);title('Whole vs. partial')

        subplot(sp(1),sp(2),4+(ss-1)*sp(2))
        if strcmp(stage,'raw')
            bar((mean(n_corrLR(data.trial_code<7,:)))); hold on;
            errorbar(1:2,(mean(n_corrLR(data.trial_code<7,:))),...
                (std(n_corrLR(data.trial_code<7,:))),'k.')
            xticklabels({'Left','Right'});xlim([.5 2.5]);
        else
            bar([mean(n_corrLR(data.trial_code<7,:)) fdata.w_index]); hold on;
            errorbar(1:3,[mean(n_corrLR(data.trial_code<7,:)) fdata.w_index],...
                [std(n_corrLR(data.trial_code<7,:)) fdata.w2_err],'k.')
            xticklabels({'Left','Right','w_{index}'});xlim([.5 3.5]);
        end
        ylim([0 max(ylim)])
        title('Left vs Right (6T)')
        
    end
end
