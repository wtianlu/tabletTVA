%% -- Visualisation of raw performance with TVA parameters
% 20201021 Copied from tablettva_final in 2020_manuscript

function visualiseTVA(p,ses,d,vis)
rawdata = d.raw; %rawdata.sesmark = findgroups(rawdata.session,rawdata.device);
fitdata = d.fit; %fitdata.sesmark = cellfun(@(x) str2double(x(5)),fitdata.ID);
                 %fitdata.p_num = cellfun(@(x) str2double(x(2:3)),fitdata.ID);

% Select and process data
% p = 6; % selected participant for manuscript: 6. strange alphas: 11 13 16 18 23 33 38 48 57
if ses == 1
data = rawdata(and(rawdata.p_num==p,rawdata.sesmark<3),:); % select session 1
else
    data = rawdata(and(rawdata.p_num==p,rawdata.sesmark>2),:); % select session 1
%     data.sesmark = data.sesmark-2;
end
perf_p = [data.trial_code, data.sesmark, ... % trial conditions, device, total/left/right
    cellfun(@(x,t) length(intersect(x,t)),data.response,data.targets)...
    cellfun(@(x,t) length(intersect(x,t(4:6))),data.response,data.targets)...
    cellfun(@(x,t) length(intersect(x,t(1:3))),data.response,data.targets)]; 

% Visualise data
close all; font = 'Arial';
figure('pos',[63   480   935   272],'DefaultTextFontName', font, 'DefaultAxesFontName', font)
% Whole report
subplot(1,5,[1:3]); hold on;
s=1+(ses-1)*2;
tmp = splitapply(@mean,perf_p(perf_p(:,2)==s,3),perf_p(perf_p(:,2)==s,1));
scatter(vis.expdur,tmp(1:6),'kd','filled')
x = [fitdata.t0(and(fitdata.p_num==p,fitdata.sesmark==s)) vis.expdur];
v = [0 fitdata{and(fitdata.p_num==p,fitdata.sesmark==s),vis.predcols(1:6)}];
xq = 0:1:210; vq = interp1(x,v,xq,'pchip','extrap'); plot(xq(11:end),vq(11:end),'k-')

s=2+(ses-1)*2;
tmp = splitapply(@mean,perf_p(perf_p(:,2)==s,3),perf_p(perf_p(:,2)==s,1));
scatter(vis.expdur,tmp(1:6),'k^','filled')
x = [fitdata.t0(and(fitdata.p_num==p,fitdata.sesmark==s)) vis.expdur];
v = [0 fitdata{and(fitdata.p_num==p,fitdata.sesmark==s),vis.predcols(1:6)}];
vq = interp1(x,v,xq,'pchip','extrap'); plot(xq,vq,'k--')
legend('Desktop (data)','Desktop (model)','Tablet (data)','Tablet (model)','location','southeast');legend boxoff
xlabel('Exposure duration (ms)'); ylabel('Correctly reported letters')

% Add parameters to figure
presentresults = [fitdata(and(fitdata.p_num==p,fitdata.sesmark==1+(ses-1)*2),vis.parcols);...
    fitdata(and(fitdata.p_num==p,fitdata.sesmark==2+(ses-1)*2),vis.parcols)];
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
s=1+(ses-1)*2;
tmp = splitapply(@mean,perf_p(perf_p(:,2)==s,3),perf_p(perf_p(:,2)==s,1));% 2T: 7:15, 2T4D: 16:24
perfpart = [mean(tmp(7:15)) mean(tmp(16:24))];
s=2+(ses-1)*2;
tmp = splitapply(@mean,perf_p(perf_p(:,2)==s,3),perf_p(perf_p(:,2)==s,1));
perfpart = [perfpart; [mean(tmp(7:15)) mean(tmp(16:24))]];
b = bar(perfpart');
b(1).FaceColor = 'k';b(2).FaceColor = ones(1,3)*.9; xlim([.5 2.5])
legend('Desktop','Tablet');xticks([1 2]);legend boxoff;xticklabels({'2T','2T4D'});xlabel('Condition');
yticks([0:4]);ylim(allylim)
presentresults = [presentresults array2table([perfpart perfpart(:,1)-perfpart(:,2)],...
    'VariableNames', strsplit('raw2T raw2T4D r2t4Ddiv2T'))];

% Hemifield
subplot(1,5,5); hold on
perfpart = [splitapply(@mean,perf_p(:,4),findgroups(perf_p(:,2))) splitapply(@mean,perf_p(:,5),findgroups(perf_p(:,2)))];
b = bar(perfpart');
b(1).FaceColor = 'k';b(2).FaceColor = ones(1,3)*.9; xlim([.5 2.5])
legend('Desktop','Tablet');xticks([1 2]);legend boxoff;xticklabels({'Left','Right'});xlabel('Hemifield');
yticks([0:4]);ylim(allylim)

% Show parameters
% presentresults = [presentresults array2table([perfpart perfpart(:,1)./(perfpart(:,1)+perfpart(:,2))],...
%     'VariableNames', strsplit('rawL rawR rawLdivSum'))];
presentresults.Properties.RowNames = strsplit('D T');
fprintf('Participant %d presentresults:\n',p);disp(presentresults);