% PROJECT:      WP1b - assessment of visual attention on a tablet device 
% AUTHOR:       Lulu Wang
% INSTITUTION:  KU Leuven
% CONTENT:      Read raw data (unity 2019 format)
% Input:        Full path to data file
% -------------------------------------------------------------------------


function data = read_tva_unity2019(fn)
    fid = fopen(fn,'r');
    tline = fgetl(fid);tlines = {};
    while ischar(tline)
        tlines{end+1} = tline;
        tline = fgetl(fid);
    end
    fclose(fid); tlines = tlines';
    clear fid tline;

    tlines = tlines(~cellfun(@isempty,tlines));

    rawdata = tlines(cellfun(@(x) ~ismember(x(1),{'#','-'}),tlines));
    ntrials = rawdata{1};
    rawdata = cellfun(@(x) strsplit(x(1:strfind(x,'kp'))),rawdata(2:end),'UniformOutput',false); rawdata = vertcat(rawdata{:});

    % dots became comma's?
    rawdata = cellfun(@(x) strrep(x,',','.'),rawdata,'UniformOutput',false);

    data = cell2table(rawdata(:,1:5),'VariableNames',...
        {'trial_code','stim_dur','targets','distractors','response'}); 
    data.trial_code = str2double(data.trial_code); data.stim_dur = str2double(data.stim_dur);
end