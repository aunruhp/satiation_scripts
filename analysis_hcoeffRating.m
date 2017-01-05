 function analysis_hcoeffRating



% % % General settings
% % rsp_2AFC               = [200,1000];   %remember, the longer the period the less sensitive, but of course, an answer could be later than 1 second
% % rsp_rating             = [200,1000];
% % rsp_choice             = [200,1000];
% % % bsl_time            = [-500,0];  	% The time before presentation of the A1 pic, which is used a reference baseline for A1, b and A2 in all (and only) the rank-sum analysis

% Rating Paradigm Settings

rsp_rating             = [200,1000];


shiftwindows=[3,5]; 

% shiftwindows variable creates the as many 500 ms shift wondows as specified, with [3,5] the data time windows analysed would look 
% something like this
% rsp_rating1             = [200,1000];
% rsp_rating2             = [200,1000]+500;
% rsp_rating3             = [200,1000]+1000;
% rsp_choice1             = [200,1000];
% rsp_choice2             = [200,1000]+500;
% rsp_choice3             = [200,1000]+1000;
% rsp_choice4             = [200,1000]+1500;
% rsp_choice5             = [200,1000]+2000;

% this allows for faster analysis as it computes one single baseline for
% each cell, and compares it with all you want. Moreover, it is
% necessary they have  exactly the same duration for this


fast=1;


%% -------  Part A: Load the Rast files and provide a name for this analysis
dbstop if error
anatitle = inputdlg('Please enter a meaningful title (file name) for this analysis','Analysis title',1,{''});
ana_title= char(anatitle{1,1});
saveffolder='/media/Projects/Alex/New Analysis/';
mkdir([saveffolder, ana_title]);
savehere = [saveffolder, ana_title];
f = {};
FileName = {};
while iscell(FileName)
    try
        cd(saveffolder);      %/Volumes/MTL/CS4 output')   %orginal form michael '/Volumes/Barnacle/Data/Satiation_data/Satiation_session_data'
    catch
        cd('/media/Projects/Alex/');      %/Volumes/MTL/CS4 output')
    end
    [FileName,PathName] = uigetfile('.mat','Multiselect','on');
    if iscell(FileName)
        disp(anatitle)
        for a = 1:length(FileName)
            f{end+1,1} = [PathName,FileName{1,a}];
        end
    elseif FileName ~= 0
        disp(savehere)
        f{end+1,1} = [PathName,FileName];
        FileName = {};
    end
end
kosher_files = true;
for a = 1:length(f)
    if isempty(strfind(f{a,1},'/Rast___'));
        kosher_files = false;
    end
end
if ~kosher_files
    errordlg('Make sure all yor files are raster files (file names start with "Rast...")!','Try again...');
    return
end
% Concatenate all rast files (but only the relevant units, TTLs and conditions)
rast_mrg = {};
rast_bsl_mrg = {};
files_rank = cell(0,6);
files_satcue = cell(0,7);


for a = 1:length(f)
    clear rast
    load(f{a,1})
    pat_ID=f{a,1}(35:37);
    [pathstr,~,~] = fileparts(f{a,1});
    modpath = strrep(pathstr,'Rast___','Files__');
    cd(modpath);
    files_dir = dir;
    for b = 1:length(files_dir)
        if ~isempty(regexp(files_dir(b,1).name,'new_sat_output_'))
%             if isempty(find(cellfun(@(x) isequal(x,[modpath,'/',files_dir(b,1).name]),files_satcue(:,6))))
                files_satcue(end+1,1:5) = rast(1,1:5);
                files_satcue(end,6) = {[modpath,'/',files_dir(b,1).name]};
                files_satcue(end,7) = {[modpath,'/alexevents_p', pat_ID, '.mat']}; % load also the events
                
%             end
        end
    end
    rast_rsp = rast(ismember(cell2mat(rast(:,9)),7) | ismember(cell2mat(rast(:,9)),43) | ismember(cell2mat(rast(:,9)),44),:); %take all rast for evetn 7 ==image || event 43==likert-2AFC
    rast_mrg = [rast_mrg;rast_rsp]; %includes not all the information of the rasts, so patient, experiment, PSTH...
%     rast_bsl = rast(ismember(cell2mat(rast(:,9)),7) & ismember(cell2mat(rast(:,10)),7),:); % choice A1 trials as baseline for the B1 analysis only ,as B and A2 have no baseline time
%     rast_bsl_mrg = [rast_bsl_mrg;rast_bsl];
end
rast_mrg(all(cellfun(@isempty,rast_mrg),2),:) = [];       % remove empty rows
% rast_bsl_mrg(all(cellfun(@isempty,rast_bsl_mrg),2),:) = [];       % remove empty rows






%% H coefficient analysis

% Each session combines all the Rasts of a session, of all regions, so all
% regions saved in same place.... change it
cd(savehere)

for sess = 1:size(files_satcue,1)
    tic
    clear rank satcue
    load(files_satcue{sess,6})
    load(files_satcue{sess,7})
    
    pat_ID=files_satcue{sess,2};
    
    % region name
    [~,regionname,~] = fileparts(f{sess,1});
    kurz=regionname(1:5);
    len=length(regexp(kurz, '[A-Z]'));
    region=regionname(1:len);
    

    

    
    % Now run over all units from this recording session:
    
    regnum=(regexp(rast_mrg(:,6),region));
    for ff=1:length(regnum)
        
        if isempty(regnum{ff})
            regnum{ff}=0;
        end
    end
    regnum=cell2mat(regnum);  % this allows to separate each region, so it has to fit the name of the session, AND it has to fit the number of the region
    
    
% % % %     n_ses_units = sum(ismember(cell2mat(rast_mrg(:,2:5)),cell2mat(files_satcue(sess,2:5)),'rows') & cell2mat(rast_mrg(:,9)) == 43 & cell2mat(rast_mrg(:,10)) == 0); 
% % % %     rast_mrg_sess = rast_mrg(ismember(cell2mat(rast_mrg(:,2:5)),cell2mat(files_satcue(sess,2:5)),'rows'),:);     


    n_ses_units = sum(ismember(cell2mat(rast_mrg(:,2:5)),cell2mat(files_satcue(sess,2:5)),'rows') & cell2mat(rast_mrg(:,9)) == 43 & cell2mat(rast_mrg(:,10)) == 0 & regnum); % Unique number of units in rast_mrg recorded in this session and this region
    rast_mrg_sess = rast_mrg(ismember(cell2mat(rast_mrg(:,2:5)),cell2mat(files_satcue(sess,2:5)),'rows') & regnum,:); % Units in rast_mrg recorded in this session and this region
    
    rast_mrg_sess_T7 = rast_mrg_sess(cell2mat(rast_mrg_sess(:,9)) == 7 & cell2mat(rast_mrg_sess(:,10)) == 0,:);  % all T7 responses, so all stimulus
    rast_mrg_sess_T43 = rast_mrg_sess(cell2mat(rast_mrg_sess(:,9)) == 43 & cell2mat(rast_mrg_sess(:,10)) == 0,:); % all T43 responses, so all 2AFC likerts
    rast_mrg_sess_T44 = rast_mrg_sess(cell2mat(rast_mrg_sess(:,9)) == 44 & cell2mat(rast_mrg_sess(:,10)) == 0,:); % all T44 responses, so all rating likerts

    
        
% %     ttl_n_cond = cell2mat(rast_mrg(:,9:10)); %icludes when the event (7 or 43 or 44) and the criteria (0-20)
    ttl_n_cond = cell2mat(rast_mrg_sess(:,9:10)); %icludes when the event (7 or 43 or 44) and the criteria (0-20)
%     judge1 = find((ttl_n_cond(:,1) == 7) & (ttl_n_cond(:,2) == 3));     %  First block judge pics
%     judge2 = find((ttl_n_cond(:,1) == 7) & (ttl_n_cond(:,2) == 4));     %  Second block judge pics
% %     judge =  find((ttl_n_cond(:,1) == 7) & (ttl_n_cond(:,2) == 1));    %  All judge pics
    judge =  find((ttl_n_cond(:,1) == 7) & (ttl_n_cond(:,2) == 1));
    
    pic_A1 = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 7)));   % Pic A1
    pic_B  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 8)));   % Pic B
    pic_A2  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 9)));

    
    
    if length(ts_index_missmatch_between_ttls_and_events_ts) == 3
        stop_index=ts_index_missmatch_between_ttls_and_events_ts(2); %point of cutting experiment (food time)
    else
        error('the timestamps to calculate the <<food time>> are not correct, an error happened, shit!!')   
    end

    for unit = 1:n_ses_units %loop over each cluster in this m
        
        % Rating
        
%         perittl_judge1=rast_mrg_sess{judge1(unit),12};
%         perittl_judge2=rast_mrg_sess{judge2(unit),12};
        perittl_judge= rast_mrg_sess{judge(unit), 12}; 
        perittl_rating=rast_mrg_sess_T44{unit,12};
        
        
        %2AFC
        
% % % % % %         perittl_A1=rast_mrg_sess{pic_A1(unit),12};
% % % % % %  tl_choice=r       perittl_B= rast_mrg_sess{pic_B(unit),12};
% % % % % %         perittl_A2=rast_mrg_sess{pic_A2(unit),12};
% % % % % %         peritast_mrg_sess_T43{:,12};

        
        
        raw2= rast_mrg_sess_T7{unit, 13};  %this contains the raw data
%         raw2=rast{(23*unit)-22, 13}; %every cluster has 23 rasts, as it has 1 for each condition in experiment. 24th column are the timestamps all spikes in this cluster
        raw= (raw2 - raw2(1))'; %raw_data have to be in a row vector and start from zero to duration of experiment
        food_time= diff(events_ts(stop_index:stop_index+1))/1e6; %time to eat, around the index of cutting, in seconds
        events_time=(events_ts(end)-events_ts(1))/1e6; %in seconds
        total_time= raw(end)/1e3;
        
%         total_time/60
%         food_time/60
%         events_time/60
%         
%         if total_time >= events_time
            final_time= total_time - food_time; % in seconds
            mfr=(length(raw))/final_time; %in Hz
%             nstims_2AFC=length(perittl_A1);
            nstims_rat= length(perittl_judge);
%             final_time/60
            
% % % % %             [response_A1{sess, unit}, response_B{sess, unit}, response_A2{sess, unit}, response_choice{sess, unit}]=hcoeff(raw, mfr, nstims_2AFC, rsp_2AFC, fast, perittl_A1,  perittl_B, perittl_A2, perittl_choice);
            [response_judge1{1, unit},response_judge2{1, unit},response_judge3{1, unit}, response_rating1{1, unit}, response_rating2{1, unit},response_rating3{1, unit},response_rating4{1, unit},response_rating5{1, unit}]=hcoeffRating(raw, mfr, nstims_rat, rsp_rating, fast, shiftwindows, perittl_judge, perittl_rating); 
            
            
%         else
%             error('the length of the recording is not longer than the length of the events... and should be so bro!')
%             
%         end
        
        
    end
    
    
    
cd(savehere)
savingname=['Hcoeff_responses_', pat_ID , '_', region];
save (savingname, 'response_judge1', 'response_judge2', 'response_judge3', 'response_rating1', 'response_rating2', 'response_rating3', 'response_rating4', 'response_rating5')
clear ('response_judge1', 'response_judge2', 'response_judge3', 'response_rating1', 'response_rating2', 'response_rating3', 'response_rating4', 'response_rating5')
toc
end

end



