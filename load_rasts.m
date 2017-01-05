function [rast_mrg, files, rast_bsl_mrg]=load_rasts(f)





% Concatenate all rast files (but only the relevant units, TTLs and conditions)
rast_mrg = {};
rast_bsl_mrg = {};
files = cell(0,6);
for a = 1:length(f)
    clear rast
    load(f{a})

    rast_rsp = rast(ismember(cell2mat(rast(:,9)),7) | ismember(cell2mat(rast(:,9)),43) | ismember(cell2mat(rast(:,9)),44),:); %take all rast for evetn 7 ==image || event 43==likert-2AFC || event 44= likert rating
    rast_mrg = [rast_mrg;rast_rsp]; %includes not all the information of the rasts, so patient, experiment, PSTH...
    files(end+1,1:5) = rast(1,1:5);
    if nargout>2
        rast_bsl = rast(ismember(cell2mat(rast(:,9)),7) & ismember(cell2mat(rast(:,10)),7),:); % choice A1 trials as baseline for the B1 analysis only ,as B and A2 have no baseline time: condition 7,7 == A1
        rast_bsl_mrg = [rast_bsl_mrg;rast_bsl];
    end
end

rast_mrg(all(cellfun(@isempty,rast_mrg),2),:) = [];       % remove empty rows

if nargout>2
    rast_bsl_mrg(all(cellfun(@isempty,rast_bsl_mrg),2),:) = [];       % remove empty rows
end

