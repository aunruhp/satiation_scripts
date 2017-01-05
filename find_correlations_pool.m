function find_correlations_pool

tic
dbstop if error

%% settings 


rsp_t_AFC= [0,1000]; %use at the beginning this time period to find the corrlations, this could change
rsp_t_likertAFC= [0,1000]; %use at the beginning this time period to find the corrlations, this could change
rsp_t_rating= [0,1000]; %use at the beginning this time period to find the corrlations, this could change
rsp_t_likertRat= [0,1000]; %use at the beginning this time period to find the corrlations, this could change
rsp_t_bsl= [-400,0];

value_signal={'All'}; % this is the value signal used in this analysis, thus, Tot, Rat, Rank or AFC.  If 'All', then it loop over the four types of signlas and saves all of them 
typeCorr= 'Spearman'; % here choose which correlation to use 



%% folder and subfolders to save results

anatitle = inputdlg('Please enter a meaningful title (file name) for this analysis','Analysis title',1,{''});
ana_title= char(anatitle{1,1});
mkdir(sprintf('/media/Projects/Alex/New Analysis/%s',ana_title));
folder_to_save = sprintf('/media/Projects/Alex/New Analysis/%s',ana_title);

% subfolders 

path_for_corr= [folder_to_save, '/Spearman_Values']; 
mkdir(path_for_corr)
path_for_ancova= [folder_to_save, '/Ancova_Values']; 
mkdir(path_for_ancova)
path_for_anovakruskal= [folder_to_save, '/Anova_Kruskal_Detection']; 
mkdir(path_for_anovakruskal)
path_for_2wayancova= [folder_to_save, '/2_way_Ancova']; 
mkdir(path_for_2wayancova)
path_for_changeVal= [folder_to_save, '/Change_value']; 
mkdir(path_for_changeVal)
 

%% 

namefiles = {};

satiationfolder = '/media/Projects/Alex/New Analysis';

cd (satiationfolder);      %/Volumes/MTL/CS4 output')   %orginal form michael '/Volumes/Barnacle/Data/Satiation_data/Satiation_session_data'



DirName=dir;
behavDir=[];
for a=1:length(DirName)
    if regexp(DirName(a).name, '[0-9][0-9][0-9]')
        if isempty(regexp(DirName(a).name, '039', 'once')) % don't analyze pat 39, he did not eat.
            if isempty(regexp(DirName(a).name, '037', 'once')) % don't analyze pat 37, he ate the wrong product
                
              
                cd (DirName(a).name)
                
                
                pathfolders=dir;
                for i=1:length(pathfolders)
                    if regexp(pathfolders(i).name, 'Rast__')
                        
                        cd (pathfolders(i).name)
                        filesfolder=dir;
                        for j=1:length(filesfolder)
                            if ~isempty(regexp(filesfolder(j).name, 'RA___', 'once')) || ~isempty(regexp(filesfolder(j).name, 'LA___', 'once'))

                                behavDir{end+1}=DirName(a).name;
%                                 load (filesfolder(j).name)
                                namefiles{end+1}=[satiationfolder, '/', DirName(a).name, '/', pathfolders(i).name '/', filesfolder(j).name];
                                
                            end
                        end
                    end
                    
                end
                cd ../..
                
            end
        end
    end
end

% behavDir(:,cellfun(@isempty,behavDir)) = []; 
% 
% [rasts, outfiles]=load_rasts(namefiles);

load amygdala_rasts_no37nor39 



if isequal(value_signal, {'All'}) 
    value_signal=[{'Rank'}, {'Rat'}, {'Tot'}, {'AFC'}];
    
    for hh=1:length(value_signal)
        main(value_signal(hh), rasts, rsp_t_AFC, rsp_t_rating, rsp_t_likertRat, rsp_t_likertAFC, rsp_t_bsl, folder_to_save, typeCorr, namefiles, outfiles, behavDir, path_for_corr, path_for_ancova, path_for_anovakruskal, path_for_2wayancova, path_for_changeVal)
    end
    
else 
    main(value_signal, rasts, rsp_t_AFC, rsp_t_rating, rsp_t_likertRat, rsp_t_likertAFC, rsp_t_bsl, folder_to_save, typeCorr, namefiles, outfiles, behavDir, path_for_corr, path_for_ancova, path_for_anovakruskal, path_for_2wayancova, path_for_changeVal)

end





function main(value_signal, rasts, rsp_t_AFC, rsp_t_rating, rsp_t_likertRat, rsp_t_likertAFC, rsp_t_bsl, folder_to_save, typeCorr, namefiles, outfiles, behavDir, path_for_corr, path_for_ancova, path_for_anovakruskal, path_for_2wayancova, path_for_changeVal)

if length(namefiles) == length(outfiles)
    
    
    % just preallocate to speed up this supoer slow stuff 
    
%     n_ses_units=size(rasts,1)/23;
%     CorrRat= nan(n_ses_units, 8);
%     CorrLikRat= nan(n_ses_units, 8);
%     CorrA1= nan(n_ses_units, 30);
%     CorrB= nan(n_ses_units, 30);
%     CorrA2= nan(n_ses_units, 30);
%     CorrLikAFC= nan(n_ses_units, 30);
%     CorrA1_B= nan(n_ses_units, 30);
%     CorrB_A2= nan(n_ses_units, 30);
%     CorrA1_Babs= nan(n_ses_units, 30);
%     CorrB_A2abs= nan(n_ses_units, 30);
%     CorrA1_A2= nan(n_ses_units, 30);
%     CorrA1_A2abs= nan(n_ses_units, 30);
%     
%     pvalRat= nan(n_ses_units, 8);
%     pvalLikRat= nan(n_ses_units, 8);
%     pvalA1= nan(n_ses_units, 30);
%     pvalB= nan(n_ses_units, 30);
%     pvalA2= nan(n_ses_units, 30);
%     pvalLikAFC= nan(n_ses_units, 30);
%     pvalA1_B= nan(n_ses_units, 30);
%     pvalB_A2= nan(n_ses_units, 30);
%     pvalA1_Babs= nan(n_ses_units, 30);
%     pvalB_A2abs= nan(n_ses_units, 30);
%     pvalA1_A2= nan(n_ses_units, 30);
%     pvalA1_A2abs= nan(n_ses_units, 30);
%     
%     CorrRatbsl= nan(n_ses_units, 8);
%     CorrLikRatbsl= nan(n_ses_units, 8);
%     CorrA1bsl= nan(n_ses_units, 30);
%     CorrBbsl= nan(n_ses_units, 30);
%     CorrA2bsl= nan(n_ses_units, 30);
%     CorrLikAFCbsl= nan(n_ses_units, 30);
%     CorrA1_Bbsl= nan(n_ses_units, 30);
%     CorrB_A2bsl= nan(n_ses_units, 30);
%     CorrA1_Babsbsl= nan(n_ses_units, 30);
%     CorrB_A2absbsl= nan(n_ses_units, 30);
%     CorrA1_A2bsl= nan(n_ses_units, 30);
%     CorrA1_A2absbsl= nan(n_ses_units, 30);
%     
%     pvalRatbsl= nan(n_ses_units, 8);
%     pvalLikRatbsl= nan(n_ses_units, 8);
%     pvalA1bsl= nan(n_ses_units, 30);
%     pvalBbsl= nan(n_ses_units, 30);
%     pvalA2bsl= nan(n_ses_units, 30);
%     pvalLikAFCbsl= nan(n_ses_units, 30);
%     pvalA1_Bbsl= nan(n_ses_units, 30);
%     pvalB_A2bsl= nan(n_ses_units, 30);
%     pvalA1_Babsbsl= nan(n_ses_units, 30);
%     pvalB_A2absbsl= nan(n_ses_units, 30);
%     pvalA1_A2bsl= nan(n_ses_units, 30);
%     pvalA1_A2absbsl= nan(n_ses_units, 30);
%     
%     PancovaRat= nan(n_ses_units, 8);
%     PancovaLikRat= nan(n_ses_units, 8);
%     PancovaA1= nan(n_ses_units, 30);
%     PancovaB= nan(n_ses_units, 30);
%     PancovaA2= nan(n_ses_units, 30);
%     PancovaLikAFC= nan(n_ses_units, 30);
%     PancovaA1_B= nan(n_ses_units, 30);
%     PancovaB_A2= nan(n_ses_units, 30);
%     PancovaA1_Babs= nan(n_ses_units, 30);
%     PancovaB_A2abs= nan(n_ses_units, 30);
%     PancovaA1_A2= nan(n_ses_units, 30);
%     PancovaA1_A2abs= nan(n_ses_units, 30);
% 
%     
%     PancovaRatbsl= nan(n_ses_units, 8);
%     PancovaLikRatbsl= nan(n_ses_units, 8);
%     PancovaA1bsl= nan(n_ses_units, 30);
%     PancovaBbsl= nan(n_ses_units, 30);
%     PancovaA2bsl= nan(n_ses_units, 30);
%     PancovaLikAFCbsl= nan(n_ses_units, 30);
%     PancovaA1_Bbsl= nan(n_ses_units, 30);
%     PancovaB_A2bsl= nan(n_ses_units, 30);
%     PancovaA1_Babsbsl= nan(n_ses_units, 30);
%     PancovaB_A2absbsl= nan(n_ses_units, 30);
%     PancovaA1_A2bsl= nan(n_ses_units, 30);
%     PancovaA1_A2absbsl= nan(n_ses_units, 30);

    %% regional analysis
    
    unitnum=1;
    
    for sess=1:length(namefiles)
        
        [~,regionname,~] = fileparts(namefiles{sess});
        kurz=regionname(1:5);
        len=length(regexp(kurz, '[A-Z]'));
        region=regionname(1:len);
        
        regnum=(regexp(rasts(:,6),region));
        for ff=1:length(regnum)

            if isempty(regnum{ff})
                regnum{ff}=0;
            end
        end
        regind=cell2mat(regnum);  % this allows to separate each region, so it has to fit the name of the session, AND it has to fit the number of the region
        unitind=ismember(cell2mat(rasts(:,2:5)),cell2mat(outfiles(sess,2:5)),'rows'); % check if the unit is part of the channel, compares the name of the patient, etc


        n_ses_units = sum(unitind & cell2mat(rasts(:,9)) == 43 & cell2mat(rasts(:,10)) == 0 & regind); % Unique number of units in rast_mrg recorded in this session and this region
        rast_mrg_sess = rasts(unitind & regind,:); % Units in rast_mrg recorded in this session and this region
% %         rast_mrg_sess_T7 = rast_mrg_sess(cell2mat(rast_mrg_sess(:,9)) == 7 & cell2mat(rast_mrg_sess(:,10)) == 0,:);  % all T7 responses, so all stimulus
% %         rast_mrg_sess_T43 = rast_mrg_sess(cell2mat(rast_mrg_sess(:,9)) == 43 & cell2mat(rast_mrg_sess(:,10)) == 0,:); % all T43 responses, so all 2AFC likerts
% %         rast_mrg_sess_T44 = rast_mrg_sess(cell2mat(rast_mrg_sess(:,9)) == 44 & cell2mat(rast_mrg_sess(:,10)) == 0,:); % all T44 responses, so all rating likerts


        ttl_n_cond = cell2mat(rast_mrg_sess(:,9:10)); %icludes when the event (7 or 43 or 44) and the criteria (0-20)

%         judge =  find((ttl_n_cond(:,1) == 7) & (ttl_n_cond(:,2) == 1));    %  All judge pics
%         pic_A1 = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 7)));   % Pic A1
%         pic_B  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 8)));   % Pic B
%         pic_A2  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 9)));
        pic_rat1 = find((ttl_n_cond(:,1) == 7) & (ttl_n_cond(:,2) == 3));     %  First block judge pics
        pic_rat2 = find((ttl_n_cond(:,1) == 7) & (ttl_n_cond(:,2) == 4));     %  Second block judge pics
        pic_A11 = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 11)));   % Pic A1
        pic_A12 = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 12)));   % Pic A1
        pic_B1  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 13)));   % Pic B
        pic_B2  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 14)));   % Pic B
        pic_A21  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 15)));
        pic_A22  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 16)));
        lik_AFC1 = find((ttl_n_cond(:,1) == 43) & ((ttl_n_cond(:,2) == 17)));
        lik_AFC2 = find((ttl_n_cond(:,1) == 43) & ((ttl_n_cond(:,2) == 18)));
        lik_Rat1 = find((ttl_n_cond(:,1) == 44) & ((ttl_n_cond(:,2) == 19)));        
        lik_Rat2 = find((ttl_n_cond(:,1) == 44) & ((ttl_n_cond(:,2) == 20)));
        
        baseline1= find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 11)));  % baseline extraction for analysis of change of firing rate, instead of mean firing rate
        baseline2 = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 12)));
        
        [satcue,ranking]=load_behavioural_data(behavDir{sess});
        
        
        % I extract here the time-series of the value. 4 input values
        % and the real output value are coded here, moreover th chosen and
        % unchosen value
            
            
            
        
        for jj=[3,6]
            
            
            % Output values (real answers)
            val_rat(:,jj/3)= cell2mat(satcue{jj-1}(:,14));          
            val_2AFC(:,jj/3)= cell2mat(satcue{jj}(:,14));

            
            
            % Input values
            
            % Rating paradigm
            
            val_rat_Rat(:,jj/3)=nan(length(satcue{jj-1}(:,14)),1);
            val_rat_AFC(:,jj/3)=nan(length(satcue{jj-1}(:,14)),1);
            val_rat_Rank(:,jj/3)=nan(length(satcue{jj-1}(:,14)),1);
            val_rat_Tot(:,jj/3)=nan(length(satcue{jj-1}(:,14)),1);
            
            
            for i=1:20
                
                ind3=find(cell2mat(satcue{jj-1}(:,2))==i);
                valRat=  ranking{jj-1}(i,4); % ranking
                val2AFC= ranking{jj}(i,4);
                valRank= ranking{jj}(i,3);
                valTotal= ranking{jj}(i,5);
                val_rat_Rat(ind3,jj/3)= valRat;
                val_rat_AFC(ind3,jj/3)= val2AFC;
                val_rat_Rank(ind3,jj/3)= valRank;
                val_rat_Tot(ind3,jj/3)= valTotal;
                
                stimuliRat(ind3,jj/3)= i;
            end
            
            
            % 2AFC paradigm
            
            val_A_Rat(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
            val_B_Rat(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on B given in previous rating trials (-300 to 300)
            %                 val_left_rat=  nan(length(satcue{jj}),1);              % value of image on the left given in previous rating trials (-300 to 300)
            %                 val_right_rat= nan(length(satcue{jj}),1);              % value of image on the right given in previous rating trials (-300 to 300)
            val_A_AFC(:,jj/3) = nan(length(satcue{jj}(:,14)),1);            % value of image on A given in previous whole 2AFC trials (0 to 1900)
            val_B_AFC(:,jj/3)= nan(length(satcue{jj}(:,14)),1);            % value of image on B given in previous whole 2AFC trials (0 to 1900)
            val_A_Rank(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
            val_B_Rank(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20)
            val_A_Tot(:,jj/3) = nan(length(satcue{jj}(:,14)),1);       % value of image on A given in previous whole 2AFC trials (0 to 1900)
            val_B_Tot(:,jj/3) = nan(length(satcue{jj}(:,14)),1);       % value of image on B given in previous whole 2AFC trials (0 to 1900)
            
            
            
            % chosen_value (C), unchosen_value (U)
            
            val_C_Rat(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
            val_U_Rat(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on B given in previous rating trials (-300 to 300)
            val_C_AFC(:,jj/3) = nan(length(satcue{jj}(:,14)),1);            % value of image on A given in previous whole 2AFC trials (0 to 1900)
            val_U_AFC(:,jj/3) = nan(length(satcue{jj}(:,14)),1);            % value of image on B given in previous whole 2AFC trials (0 to 1900)
            val_C_Rank(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
            val_U_Rank(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20)
            val_C_Tot(:,jj/3) = nan(length(satcue{jj}(:,14)),1);       % value of image on A given in previous whole 2AFC trials (0 to 1900)
            val_U_Tot(:,jj/3) = nan(length(satcue{jj}(:,14)),1);       % value of image on B given in previous whole 2AFC trials (0 to 1900)
            
            
            for i=1:20
                
                indl = find(cell2mat(satcue{jj}(:,2))==i);            % when this stimulus was on left
                indr = find(cell2mat(satcue{jj}(:,3))==i);           % when this stimulus was on right
                lr_to_ab = cell2mat(satcue{jj}(indl,8));              % conversion factor from left_rigth to A_B, this factor saved on the column 8
                lr_to_ab2 = cell2mat(satcue{jj}(indr,8));            % conversion factor from left_rigth to A_B
                
                lr_to_cu = cell2mat(satcue{jj}(:,14))>0;              % conversion factor from left_rigth to A_B, just if positive or negative the output, 1 if right chosen, zero if left chosen
                
                
                indA = cat(1,((lr_to_ab) .* indl), (~lr_to_ab2).* indr);  indA(indA==0)= [];  %indices for A
                indB = cat(1,((~lr_to_ab) .* indl), (lr_to_ab2).* indr);  indB(indB==0)= [];  %indices for B
                
                indC = cat(1,indl(~lr_to_cu(indl)), indr(lr_to_cu(indr)));   % indeces for chosen, 1 if right chosen, so 1*indr, and zero if left chosen, so ~0 * indl
                indU = cat(1,indl(lr_to_cu(indl)), indr(~lr_to_cu(indr)));
                
                
                valRat = ranking{jj-1}(i,4);                        % rating value (-300 to 300)
                valRank = ranking{jj}(i,3);                         % ranking (0 to 20)
                val2AFC = ranking{jj}(i,4);                         % choice rating (0-1900)
                val2AFCtotal = ranking{jj}(i,5);                         % choice rating (-1900-1900)
                
                
                val_A_Rat(indA, jj/3) = valRat;
                val_B_Rat(indB, jj/3) = valRat;
                %         val_left_rat(ind)= valRat;
                %         val_right_rat(ind2)= valRat;
                val_A_Rank(indA, jj/3) = valRank;
                val_B_Rank(indB, jj/3) = valRank;
                val_A_AFC(indA, jj/3) = val2AFC;
                val_B_AFC(indB, jj/3) = val2AFC;
                val_A_Tot(indA, jj/3) = val2AFCtotal;
                val_B_Tot(indB, jj/3) = val2AFCtotal;
                
                val_C_Rat(indC, jj/3) = valRat;
                val_U_Rat(indU, jj/3) = valRat;
                val_C_Rank(indC, jj/3) = valRank;
                val_U_Rank(indU, jj/3) = valRank;
                val_C_AFC(indC, jj/3) = val2AFC;
                val_U_AFC(indU, jj/3) = val2AFC;
                val_C_Tot(indC, jj/3) = val2AFCtotal;
                val_U_Tot(indU, jj/3) = val2AFCtotal;
                
                stimuliA(indA,jj/3)= i;
                stimuliB(indB,jj/3)= i;
                
            end
                      
        end
        
        

        
        for unit=1:n_ses_units
            
            % load neural responses indeces, before and after lunch
            
            ratResp1 = rasts{pic_rat1(unit),12};
            ratResp2 = rasts{pic_rat2(unit),12};
            likRat1 = rasts{lik_Rat1(unit),12};
            likRat2 = rasts{lik_Rat2(unit),12};  
            A11 = rasts{pic_A11(unit),12};
            A12 = rasts{pic_A12(unit),12};
            B1 = rasts{pic_B1(unit),12};
            B2 = rasts{pic_B2(unit),12};
            A21 = rasts{pic_A21(unit),12};
            A22 = rasts{pic_A22(unit),12};
            likAFC1 = rasts{lik_AFC1(unit),12};
            likAFC2 = rasts{lik_AFC2(unit),12};
            bsl1 = rasts{baseline1(unit),12};
            bsl2 = rasts{baseline2(unit),12};
            
            
            % Mean Firing rates at each time window
            
            for stim=1:190
                
                mfrB1(stim,1) = length(B1{stim}(1,B1{stim} > rsp_t_AFC(1,1) & B1{stim} < rsp_t_AFC(1,2))) * (1000/abs(rsp_t_AFC(1,2) - rsp_t_AFC(1,1))); % Hz of this unit at each stimulus at this time window
                mfrB2(stim,1) = length(B2{stim}(1,B2{stim} > rsp_t_AFC(1,1) & B2{stim} < rsp_t_AFC(1,2))) * (1000/abs(rsp_t_AFC(1,2) - rsp_t_AFC(1,1))) ; % Hz of this unit at each stimulus at this time window
                mfrA11(stim,1) = length(A11{stim}(1,A11{stim} > rsp_t_AFC(1,1) & A11{stim} < rsp_t_AFC(1,2))) * (1000/abs(rsp_t_AFC(1,2) - rsp_t_AFC(1,1))) ; % Hz of this unit at each stimulus at this time window
                mfrA12(stim,1) = length(A12{stim}(1,A12{stim} > rsp_t_AFC(1,1) & A12{stim} < rsp_t_AFC(1,2))) * (1000/abs(rsp_t_AFC(1,2) - rsp_t_AFC(1,1))); % Hz of this unit at each stimulus at this time window
                mfrA21(stim,1) = length(A21{stim}(1,A21{stim} > rsp_t_AFC(1,1) & A21{stim} < rsp_t_AFC(1,2))) * (1000/abs(rsp_t_AFC(1,2) - rsp_t_AFC(1,1))) ; % Hz of this unit at each stimulus at this time window
                mfrA22(stim,1) = length(A22{stim}(1,A22{stim} > rsp_t_AFC(1,1) & A22{stim} < rsp_t_AFC(1,2))) * (1000/abs(rsp_t_AFC(1,2) - rsp_t_AFC(1,1))) ; % Hz of this unit at each stimulus at this time window
                mfrlikAFC1(stim,1) = length(likAFC1{stim}(1,likAFC1{stim} > rsp_t_likertAFC(1,1) & likAFC1{stim} < rsp_t_likertAFC(1,2))) * (1000/abs(rsp_t_likertAFC(1,2) - rsp_t_likertAFC(1,1))) ; % Hz of this unit at each stimulus at this time window
                mfrlikAFC2(stim,1) = length(likAFC2{stim}(1,likAFC2{stim} > rsp_t_likertAFC(1,1) & likAFC2{stim} < rsp_t_likertAFC(1,2))) * (1000/abs(rsp_t_likertAFC(1,2) - rsp_t_likertAFC(1,1))); % Hz of this unit at each stimulus at this time window
                
                MFRbsl1(stim,1) = length(bsl1{stim}(1,bsl1{stim} > rsp_t_bsl(1,1) & bsl1{stim} < rsp_t_bsl(1,2))) * (1000/abs(rsp_t_bsl(1,2) - rsp_t_bsl(1,1))); % Hz of this unit at each stimulus at this time window
                MFRbsl2(stim,1) = length(bsl2{stim}(1,bsl2{stim} > rsp_t_bsl(1,1) & bsl2{stim} < rsp_t_bsl(1,2))) * (1000/abs(rsp_t_bsl(1,2) - rsp_t_bsl(1,1))); % Hz of this unit at each stimulus at this time window
                
                
            end
           
            for stim=1:60
               
                mfrRat1(stim,1) = length(ratResp1{stim}(1,ratResp1{stim} > rsp_t_rating(1,1) & ratResp1{stim} < rsp_t_rating(1,2))) * (1000/abs(rsp_t_rating(1,2) - rsp_t_rating(1,1))); % Hz of this unit at each stimulus at this time window
                mfrRat2(stim,1) = length(ratResp2{stim}(1,ratResp2{stim} > rsp_t_rating(1,1) & ratResp2{stim} < rsp_t_rating(1,2))) * (1000/abs(rsp_t_rating(1,2) - rsp_t_rating(1,1))); % Hz of this unit at each stimulus at this time window
                mfrlikRat1(stim,1) = length(likRat1{stim}(1,likRat1{stim} > rsp_t_likertRat(1,1) & likRat1{stim} < rsp_t_likertRat(1,2))) * (1000/abs(rsp_t_likertRat(1,2) - rsp_t_likertRat(1,1))); % Hz of this unit at each stimulus at this time window
                mfrlikRat2(stim,1) = length(likRat2{stim}(1,likRat2{stim} > rsp_t_likertRat(1,1) & likRat2{stim} < rsp_t_likertRat(1,2))) * (1000/abs(rsp_t_likertRat(1,2) - rsp_t_likertRat(1,1))); % Hz of this unit at each stimulus at this time window
                
                bslRat1(stim,1) = length(ratResp1{stim}(1,ratResp1{stim} > rsp_t_bsl(1,1) & ratResp1{stim} < rsp_t_bsl(1,2))) * (1000/abs(rsp_t_bsl(1,2) - rsp_t_bsl(1,1)));
                bslRat2(stim,1) = length(ratResp2{stim}(1,ratResp2{stim} > rsp_t_bsl(1,1) & ratResp2{stim} < rsp_t_bsl(1,2))) * (1000/abs(rsp_t_bsl(1,2) - rsp_t_bsl(1,1)));
            end
            
            
             
            % change in fiiring rate compared to baseline
            
            
            mfrB1bsl = mfrB1 - MFRbsl1;
            mfrB2bsl = mfrB2 - MFRbsl2;
            mfrA11bsl = mfrA11 - MFRbsl1;
            mfrA12bsl = mfrA12 - MFRbsl2;
            mfrA21bsl = mfrA21 - MFRbsl1;
            mfrA22bsl = mfrA22 - MFRbsl2;
            mfrlikAFC1bsl = mfrlikAFC1 - MFRbsl1;
            mfrlikAFC2bsl = mfrlikAFC2 - MFRbsl2;
            

            mfrRat1bsl = mfrRat1 - bslRat1;
            mfrRat2bsl = mfrRat2 - bslRat2;
            mfrlikRat1bsl =mfrlikRat1 - bslRat1;
            mfrlikRat2bsl =mfrlikRat2 - bslRat2;
            
            
            

                
      %%
            
            
            
            % select the values to use
            
            eval(['val_A = val_A_', cell2mat(value_signal)]);
            eval(['val_B = val_B_', cell2mat(value_signal)]);
            eval(['val_C = val_C_', cell2mat(value_signal)]);
            eval(['val_U = val_U_', cell2mat(value_signal)]);
            eval(['val_rat_in = val_rat_', cell2mat(value_signal)]);

            
            
            % make correlation matrices for each trial (190 for AFC, and 60 for rating)
            
            [CorrRat(unitnum, :), pvalRat(unitnum, :)] = find_corr_rating (mfrRat1, mfrRat2, val_rat, val_rat_in, typeCorr);
            [CorrLikRat(unitnum, :), pvalLikRat(unitnum, :)] = find_corr_rating (mfrlikRat1, mfrlikRat2, val_rat, val_rat_in, typeCorr);
            [CorrA1(unitnum, :), pvalA1(unitnum, :)] = find_corr_2AFC  (mfrA11, mfrA12, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrB(unitnum, :), pvalB(unitnum, :)]   = find_corr_2AFC  (mfrB1,  mfrB2, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrA2(unitnum, :), pvalA2(unitnum, :)]  = find_corr_2AFC  (mfrA21, mfrA22, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrLikAFC(unitnum, :), pvalLikAFC(unitnum, :)] = find_corr_2AFC  (mfrlikAFC1, mfrlikAFC2, val_2AFC, val_A, val_B, val_C_Rank, val_U, typeCorr);
            
            [CorrA1_B(unitnum, :), pvalA1_B(unitnum, :)] = find_corr_2AFC  (mfrA11-mfrB1, mfrA12-mfrB2, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrB_A2(unitnum, :), pvalB_A2(unitnum, :)] = find_corr_2AFC  (mfrB1-mfrA21, mfrB2-mfrA22, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrA1_Babs(unitnum, :), pvalA1_Babs(unitnum, :)]  = find_corr_2AFC  (abs(mfrA11-mfrB1), abs(mfrA12-mfrB2), val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrB_A2abs(unitnum, :), pvalB_A2abs(unitnum, :) ] = find_corr_2AFC  (abs(mfrB1-mfrA21), abs(mfrB2-mfrA22), val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            
            [CorrA1_A2(unitnum, :), pvalA1_A2(unitnum, :)] = find_corr_2AFC  (mfrA11-mfrA21, mfrA12-mfrA22, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrA1_A2abs(unitnum, :), pvalA1_A2abs(unitnum, :) ] = find_corr_2AFC  (abs(mfrA11-mfrA21), abs(mfrA12-mfrA22), val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            
           
            % baseline subtraction
            [CorrRatbsl(unitnum, :), pvalRatbsl(unitnum, :)] = find_corr_rating (mfrRat1bsl, mfrRat2bsl, val_rat, val_rat_in, typeCorr);
            [CorrLikRatbsl(unitnum, :), pvalLikRatbsl(unitnum, :)] = find_corr_rating (mfrlikRat1bsl, mfrlikRat2bsl, val_rat, val_rat_in, typeCorr);
            [CorrA1bsl(unitnum, :), pvalA1bsl(unitnum, :)] = find_corr_2AFC  (mfrA11bsl, mfrA12bsl, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrBbsl(unitnum, :), pvalBbsl(unitnum, :)]   = find_corr_2AFC  (mfrB1bsl,  mfrB2bsl, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrA2bsl(unitnum, :), pvalA2bsl(unitnum, :)]  = find_corr_2AFC  (mfrA21bsl, mfrA22bsl, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrLikAFCbsl(unitnum, :), pvalLikAFCbsl(unitnum, :)] = find_corr_2AFC  (mfrlikAFC1bsl, mfrlikAFC2bsl, val_2AFC, val_A, val_B, val_C_Rank, val_U, typeCorr);
            
            [CorrA1_Bbsl(unitnum, :), pvalA1_Bbsl(unitnum, :)] = find_corr_2AFC  (mfrA11bsl-mfrB1bsl, mfrA12bsl-mfrB2bsl, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrB_A2bsl(unitnum, :), pvalB_A2bsl(unitnum, :)] = find_corr_2AFC  (mfrB1bsl-mfrA21bsl, mfrB2bsl-mfrA22bsl, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrA1_Babsbsl(unitnum, :), pvalA1_Babsbsl(unitnum, :)]  = find_corr_2AFC  (abs(mfrA11bsl-mfrB1bsl), abs(mfrA12bsl-mfrB2bsl), val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrB_A2absbsl(unitnum, :), pvalB_A2absbsl(unitnum, :) ] = find_corr_2AFC  (abs(mfrB1bsl-mfrA21bsl), abs(mfrB2bsl-mfrA22bsl), val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            
            [CorrA1_A2bsl(unitnum, :), pvalA1_A2bsl(unitnum, :)] = find_corr_2AFC  (mfrA11bsl-mfrA21bsl, mfrA12bsl-mfrA22bsl, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrA1_A2absbsl(unitnum, :), pvalA1_A2absbsl(unitnum, :) ] = find_corr_2AFC  (abs(mfrA11bsl-mfrA21bsl), abs(mfrA12bsl-mfrA22bsl), val_2AFC, val_A, val_B, val_C, val_U, typeCorr);

            
            %% ancova analysis
            
            
            [PancovaRat(unitnum,:)] = find_ancova_rating (mfrRat1, mfrRat2, val_rat, val_rat_in);
            [PancovaLikRat(unitnum,:)] = find_ancova_rating (mfrlikRat1, mfrlikRat2, val_rat, val_rat_in);
            [PancovaA1(unitnum,:)] = find_ancova_2AFC  (mfrA11, mfrA12, val_2AFC, val_A, val_B, val_C, val_U);
            [PancovaB(unitnum,:)]   = find_ancova_2AFC  (mfrB1,  mfrB2, val_2AFC, val_A, val_B, val_C, val_U);
            [PancovaA2(unitnum,:)]  = find_ancova_2AFC  (mfrA21, mfrA22, val_2AFC, val_A, val_B, val_C, val_U);
            [PancovaLikAFC(unitnum,:)] = find_ancova_2AFC  (mfrlikAFC1, mfrlikAFC2, val_2AFC, val_A, val_B, val_C_Rank, val_U);
            
            [PancovaA1_B(unitnum,:)] = find_ancova_2AFC  (mfrA11-mfrB1, mfrA12-mfrB2, val_2AFC, val_A, val_B, val_C, val_U);
            [PancovaB_A2(unitnum,:)] = find_ancova_2AFC  (mfrB1-mfrA21, mfrB2-mfrA22, val_2AFC, val_A, val_B, val_C, val_U);
            [PancovaA1_Babs(unitnum,:)]  = find_ancova_2AFC  (abs(mfrA11-mfrB1), abs(mfrA12-mfrB2), val_2AFC, val_A, val_B, val_C, val_U);
            [PancovaB_A2abs(unitnum,:)] = find_ancova_2AFC  (abs(mfrB1-mfrA21), abs(mfrB2-mfrA22), val_2AFC, val_A, val_B, val_C, val_U);
            
            [PancovaA1_A2(unitnum,:)]  = find_ancova_2AFC  (mfrA11-mfrA21, mfrA12-mfrA22, val_2AFC, val_A, val_B, val_C, val_U);
            [PancovaA1_A2abs(unitnum,:)] = find_ancova_2AFC  (abs(mfrA11-mfrA21), abs(mfrA12-mfrA22), val_2AFC, val_A, val_B, val_C, val_U);
            
            
           
            % Baseline subtraction
            
            [PancovaRatbsl(unitnum,:)] = find_ancova_rating (mfrRat1bsl, mfrRat2bsl, val_rat, val_rat_in);
            [PancovaLikRatbsl(unitnum,:)] = find_ancova_rating (mfrlikRat1bsl, mfrlikRat2bsl, val_rat, val_rat_in);
            [PancovaA1bsl(unitnum,:)] = find_ancova_2AFC  (mfrA11bsl, mfrA12bsl, val_2AFC, val_A, val_B, val_C, val_U);
            [PancovaBbsl(unitnum,:)]   = find_ancova_2AFC  (mfrB1bsl,  mfrB2bsl, val_2AFC, val_A, val_B, val_C, val_U);
            [PancovaA2bsl(unitnum,:)]  = find_ancova_2AFC  (mfrA21bsl, mfrA22bsl, val_2AFC, val_A, val_B, val_C, val_U);
            [PancovaLikAFCbsl(unitnum,:)] = find_ancova_2AFC  (mfrlikAFC1bsl, mfrlikAFC2bsl, val_2AFC, val_A, val_B, val_C_Rank, val_U);
            
            [PancovaA1_Bbsl(unitnum,:)] = find_ancova_2AFC  (mfrA11bsl-mfrB1bsl, mfrA12bsl-mfrB2bsl, val_2AFC, val_A, val_B, val_C, val_U);
            [PancovaB_A2bsl(unitnum,:)] = find_ancova_2AFC  (mfrB1bsl-mfrA21bsl, mfrB2bsl-mfrA22bsl, val_2AFC, val_A, val_B, val_C, val_U);
            [PancovaA1_Babsbsl(unitnum,:)]  = find_ancova_2AFC  (abs(mfrA11bsl-mfrB1bsl), abs(mfrA12bsl-mfrB2bsl), val_2AFC, val_A, val_B, val_C, val_U);
            [PancovaB_A2absbsl(unitnum,:)] = find_ancova_2AFC  (abs(mfrB1bsl-mfrA21bsl), abs(mfrB2bsl-mfrA22bsl), val_2AFC, val_A, val_B, val_C, val_U);
            
            [PancovaA1_A2bsl(unitnum,:)]  = find_ancova_2AFC  (mfrA11bsl-mfrA21bsl, mfrA12bsl-mfrA22bsl, val_2AFC, val_A, val_B, val_C, val_U);
            [PancovaA1_A2absbsl(unitnum,:)] = find_ancova_2AFC  (abs(mfrA11bsl-mfrA21bsl), abs(mfrA12bsl-mfrA22bsl), val_2AFC, val_A, val_B, val_C, val_U);
            
            
            
            %% 2 way ancovas, just in case 
            
            
            [Panova2wayA1(unitnum,:)]= find_2wayancova_2AFC (mfrA11, mfrA12, val_2AFC, val_A, val_B);
            [Panova2wayB(unitnum,:)]= find_2wayancova_2AFC (mfrB1, mfrB2, val_2AFC, val_A, val_B);
            [Panova2wayA2(unitnum,:)]= find_2wayancova_2AFC (mfrA21, mfrA22, val_2AFC, val_A, val_B);
            [Panova2wayLikAFC(unitnum,:)]= find_2wayancova_2AFC (mfrlikAFC1, mfrlikAFC2, val_2AFC, val_A, val_B);
            
            [Panova2wayA1bsl(unitnum,:)]= find_2wayancova_2AFC (mfrA11bsl, mfrA12bsl, val_2AFC, val_A, val_B);
            [Panova2wayBbsl(unitnum,:)]= find_2wayancova_2AFC (mfrB1bsl, mfrB2bsl, val_2AFC, val_A, val_B);
            [Panova2wayA2bsl(unitnum,:)]= find_2wayancova_2AFC (mfrA21bsl, mfrA22bsl, val_2AFC, val_A, val_B);
            [Panova2wayLikAFCbsl(unitnum,:)]= find_2wayancova_2AFC (mfrlikAFC1bsl, mfrlikAFC2bsl, val_2AFC, val_A, val_B);
            
            
            
            %% Anova and kruskal-wallis detection of stimuli-modulated cells
            
            [Panova(unitnum,:), Pkruskal(unitnum,:)]= anova_kruskall_detection (stimuliRat, stimuliA, stimuliB, mfrRat1, mfrRat2, mfrlikRat1, mfrlikRat2, mfrA11, mfrA12, mfrB1, mfrB2, mfrA21, mfrA22, mfrlikAFC1, mfrlikAFC2);
            
            
            
            %%  Change in value analysis (before vs after lunch)
            
            
            % correlations 
            
            [CorrRat_change(unitnum,:), pvalRat_change(unitnum,:)]=change_valueRatingCorr(mfrRat1, mfrRat2, val_rat_in, typeCorr);
            [CorrLikRat_change(unitnum,:), pvalLikRat_change(unitnum,:)]=change_valueRatingCorr(mfrlikRat1, mfrlikRat2, val_rat_in, typeCorr);
            [CorrA1_change(unitnum,:), pvalA1_change(unitnum,:)]= change_valueAFCCorr(mfrA11, mfrA12, val_A, val_B, typeCorr);
            [CorrB_change(unitnum,:), pvalB_change(unitnum,:)]= change_valueAFCCorr(mfrB1, mfrB2, val_A, val_B, typeCorr);
            [CorrA2_change(unitnum,:), pvalA2_change(unitnum,:)]= change_valueAFCCorr(mfrA21, mfrA22, val_A, val_B, typeCorr);
            [CorrLikAFC_change(unitnum,:), pvalLikAFC_change(unitnum,:)]= change_valueAFCCorr(mfrlikAFC1, mfrlikAFC2, val_A, val_B, typeCorr);
            
            [CorrRat_changebsl(unitnum,:), pvalRat_changebsl(unitnum,:)]=change_valueRatingCorr(mfrRat1bsl, mfrRat2bsl, val_rat_in, typeCorr);
            [CorrLikRat_changebsl(unitnum,:), pvalLikRat_changebsl(unitnum,:)]=change_valueRatingCorr(mfrlikRat1bsl, mfrlikRat2bsl, val_rat_in, typeCorr);
            [CorrA1_changebsl(unitnum,:), pvalA1_changebsl(unitnum,:)]= change_valueAFCCorr(mfrA11bsl, mfrA12bsl, val_A, val_B, typeCorr);
            [CorrB_changebsl(unitnum,:), pvalB_changebsl(unitnum,:)]= change_valueAFCCorr(mfrB1bsl, mfrB2bsl, val_A, val_B, typeCorr);
            [CorrA2_changebsl(unitnum,:), pvalA2_changebsl(unitnum,:)]= change_valueAFCCorr(mfrA21bsl, mfrA22bsl, val_A, val_B, typeCorr);
            [CorrLikAFC_changebsl(unitnum,:), pvalLikAFC_changebsl(unitnum,:)]= change_valueAFCCorr(mfrlikAFC1bsl, mfrlikAFC2bsl, val_A, val_B, typeCorr);
            
            %ancovas
            
            [PancovaRat_change(unitnum,:)]=change_valueRatingAncova(mfrRat1, mfrRat2, val_rat_in);
            [PancovaLikRat_change(unitnum,:)]=change_valueRatingAncova(mfrlikRat1, mfrlikRat2, val_rat_in);
            [PancovaA1_change(unitnum,:)]= change_valueAFCAncova(mfrA11, mfrA12, val_A, val_B);
            [PancovaB_change(unitnum,:)]= change_valueAFCAncova(mfrB1, mfrB2, val_A, val_B);
            [PancovaA2_change(unitnum,:)]= change_valueAFCAncova(mfrA21, mfrA22, val_A, val_B);
            [PancovaLikAFC_change(unitnum,:)]= change_valueAFCAncova(mfrlikAFC1, mfrlikAFC2, val_A, val_B);
            
            [PancovaRat_changebsl(unitnum,:)]=change_valueRatingAncova(mfrRat1bsl, mfrRat2bsl, val_rat_in);
            [PancovaLikRat_changebsl(unitnum,:)]=change_valueRatingAncova(mfrlikRat1bsl, mfrlikRat2bsl, val_rat_in);
            [PancovaA1_changebsl(unitnum,:)]= change_valueAFCAncova(mfrA11bsl, mfrA12bsl, val_A, val_B);
            [PancovaB_changebsl(unitnum,:)]= change_valueAFCAncova(mfrB1bsl, mfrB2bsl, val_A, val_B);
            [PancovaA2_changebsl(unitnum,:)]= change_valueAFCAncova(mfrA21bsl, mfrA22bsl, val_A, val_B);
            [PancovaLikAFC_changebsl(unitnum,:)]= change_valueAFCAncova(mfrlikAFC1bsl, mfrlikAFC2bsl, val_A, val_B);
            
            
            
            unitnum = unitnum+1;
        end
        
        
        
        
    end

% cd (folder_to_save)



name_to_saveCorr=[path_for_corr, '/', cell2mat(value_signal), '_pool_' typeCorr '_' num2str(rsp_t_AFC(1)) '_' num2str(rsp_t_AFC(2)) '_' num2str(rsp_t_likertAFC(1)) '_' num2str(rsp_t_likertAFC(2)) '_' num2str(rsp_t_rating(1)) '_' num2str(rsp_t_rating(2)), '_' num2str(rsp_t_likertRat(1)) '_' num2str(rsp_t_likertRat(2))];
name_to_saveCorrbsl=[name_to_saveCorr '_' num2str(rsp_t_bsl(1)) '_' num2str(rsp_t_bsl(2)), 'bsl'];
name_to_saveAncova=[path_for_ancova, '/',cell2mat(value_signal), '_pool_ancova_' num2str(rsp_t_AFC(1)) '_' num2str(rsp_t_AFC(2)) '_' num2str(rsp_t_likertAFC(1)) '_' num2str(rsp_t_likertAFC(2)) '_' num2str(rsp_t_rating(1)) '_' num2str(rsp_t_rating(2)), '_' num2str(rsp_t_likertRat(1)) '_' num2str(rsp_t_likertRat(2))];
name_to_saveAncovabsl=[name_to_saveAncova '_' num2str(rsp_t_bsl(1)) '_' num2str(rsp_t_bsl(2)), 'bsl'];

name_to_saveAnovaKruskal=[path_for_anovakruskal, '/', cell2mat(value_signal), '_pool_anovakruskal_' num2str(rsp_t_AFC(1)) '_' num2str(rsp_t_AFC(2)) '_' num2str(rsp_t_likertAFC(1)) '_' num2str(rsp_t_likertAFC(2)) '_' num2str(rsp_t_rating(1)) '_' num2str(rsp_t_rating(2)), '_' num2str(rsp_t_likertRat(1)) '_' num2str(rsp_t_likertRat(2))];
name_to_save2wayAncova=[path_for_2wayancova, '/', cell2mat(value_signal), '_pool_ancova2way_' num2str(rsp_t_AFC(1)) '_' num2str(rsp_t_AFC(2)) '_' num2str(rsp_t_likertAFC(1)) '_' num2str(rsp_t_likertAFC(2)) '_' num2str(rsp_t_rating(1)) '_' num2str(rsp_t_rating(2)), '_' num2str(rsp_t_likertRat(1)) '_' num2str(rsp_t_likertRat(2)) '_' num2str(rsp_t_bsl(1)) '_' num2str(rsp_t_bsl(2)), 'bsl'];

name_to_savechangeValCorr=[path_for_changeVal, '/',cell2mat(value_signal), '_pool_change_' typeCorr '_' num2str(rsp_t_AFC(1)) '_' num2str(rsp_t_AFC(2)) '_' num2str(rsp_t_likertAFC(1)) '_' num2str(rsp_t_likertAFC(2)) '_' num2str(rsp_t_rating(1)) '_' num2str(rsp_t_rating(2)), '_' num2str(rsp_t_likertRat(1)) '_' num2str(rsp_t_likertRat(2)) '_' num2str(rsp_t_bsl(1)) '_' num2str(rsp_t_bsl(2)), 'bsl'];
name_to_savechangeValAncova=[path_for_changeVal, '/', cell2mat(value_signal), '_pool_change_ancova_' num2str(rsp_t_AFC(1)) '_' num2str(rsp_t_AFC(2)) '_' num2str(rsp_t_likertAFC(1)) '_' num2str(rsp_t_likertAFC(2)) '_' num2str(rsp_t_rating(1)) '_' num2str(rsp_t_rating(2)), '_' num2str(rsp_t_likertRat(1)) '_' num2str(rsp_t_likertRat(2)) '_' num2str(rsp_t_bsl(1)) '_' num2str(rsp_t_bsl(2)), 'bsl'];



save (name_to_saveCorr, 'CorrRat', 'CorrLikRat', 'CorrA1', 'CorrB', 'CorrA2', 'CorrLikAFC', 'CorrA1_B', 'CorrB_A2', 'CorrA1_Babs', 'CorrB_A2abs',  'pvalRat', 'pvalLikRat', 'pvalA1', 'pvalB', 'pvalA2', 'pvalLikAFC', 'pvalA1_B', 'pvalB_A2', 'pvalA1_Babs', 'pvalB_A2abs', 'CorrA1_A2', 'pvalA1_A2', 'CorrA1_A2abs', 'pvalA1_A2abs');
save (name_to_saveCorrbsl, 'CorrRatbsl', 'CorrLikRatbsl', 'CorrA1bsl', 'CorrBbsl', 'CorrA2bsl', 'CorrLikAFCbsl', 'CorrA1_Bbsl', 'CorrB_A2bsl', 'CorrA1_Babsbsl', 'CorrB_A2absbsl', 'pvalRatbsl', 'pvalLikRatbsl', 'pvalA1bsl', 'pvalBbsl', 'pvalA2bsl', 'pvalLikAFCbsl', 'pvalA1_Bbsl', 'pvalB_A2bsl', 'pvalA1_Babsbsl', 'pvalB_A2absbsl', 'CorrA1_A2bsl', 'pvalA1_A2bsl', 'CorrA1_A2absbsl', 'pvalA1_A2absbsl');
save (name_to_saveAncova, 'PancovaRat', 'PancovaLikRat', 'PancovaA1', 'PancovaB', 'PancovaA2', 'PancovaLikAFC', 'PancovaA1_B', 'PancovaB_A2', 'PancovaA1_Babs', 'PancovaB_A2abs', 'PancovaA1_A2', 'PancovaA1_A2abs');
save (name_to_saveAncovabsl, 'PancovaRatbsl', 'PancovaLikRatbsl', 'PancovaA1bsl', 'PancovaBbsl', 'PancovaA2bsl', 'PancovaLikAFCbsl', 'PancovaA1_Bbsl', 'PancovaB_A2bsl', 'PancovaA1_Babsbsl', 'PancovaB_A2absbsl', 'PancovaA1_A2bsl', 'PancovaA1_A2absbsl');

save (name_to_saveAnovaKruskal, 'Panova', 'Pkruskal');
save (name_to_save2wayAncova, 'Panova2wayA1', 'Panova2wayB', 'Panova2wayA2', 'Panova2wayLikAFC', 'Panova2wayA1bsl', 'Panova2wayBbsl', 'Panova2wayA2bsl', 'Panova2wayLikAFCbsl');

save (name_to_savechangeValCorr, 'CorrRat_change', 'CorrLikRat_change', 'CorrA1_change', 'CorrB_change', 'CorrA2_change', 'CorrLikAFC_change', 'CorrRat_changebsl', 'CorrLikRat_changebsl', 'CorrA1_changebsl', 'CorrB_changebsl', 'CorrA2_changebsl', 'CorrLikAFC_changebsl', 'pvalRat_change', 'pvalLikRat_change', 'pvalA1_change', 'pvalB_change', 'pvalA2_change', 'pvalLikAFC_change', 'pvalRat_changebsl', 'pvalLikRat_changebsl', 'pvalA1_changebsl', 'pvalB_changebsl', 'pvalA2_changebsl', 'pvalLikAFC_changebsl')
save (name_to_savechangeValAncova, 'PancovaRat_change', 'PancovaLikRat_change', 'PancovaA1_change', 'PancovaB_change', 'PancovaA2_change', 'PancovaLikAFC_change', 'PancovaRat_changebsl', 'PancovaLikRat_changebsl', 'PancovaA1_changebsl', 'PancovaB_changebsl', 'PancovaA2_changebsl', 'PancovaLikAFC_changebsl')



%% saves a vectro with tha names of the regions of each cell, and a second one with the patient ID

regionname= rasts(:,6);
regnum=(regexp(rasts(:,6),region));
for ff=1:length(regnum)
    
    if isempty(regnum{ff})
        regnum{ff}=0;
    end
end
regvector=cell2mat(regnum);

patientID=rasts(:,2);


save ([folder_to_save, '/PatientRegion_info'], 'regionname', 'regvector', 'patientID', 'region')

%% 

% cd ..
toc      
    
else
    error('Not the same length of behavioural and neural data')
end




function [CorrRat, pvalRat]= find_corr_rating (mfrRat1, mfrRat2, val_rat_out, val_rat_in, typeCorr)

% val out is always the answer from the patient 
% val in is the input, so the expected value from the stimulus, and can be
% coded as as 4 types, ranting, ranking, AFC and total

CorrRat=nan(1, 8);
pvalRat=nan(1, 8);

    [CorrRat(1, 1), pvalRat(1, 1)] =  corr([mfrRat1(:);mfrRat2(:)], val_rat_out(:), 'type', typeCorr);
    [CorrRat(1, 2), pvalRat(1, 2)] =  corr([mfrRat2(:);mfrRat1(:)], val_rat_out(:), 'type', typeCorr);
    [CorrRat(1, 3), pvalRat(1, 3)] =  corr([mfrRat1(:);mfrRat2(:)], val_rat_in(:), 'type', typeCorr);
    [CorrRat(1, 4), pvalRat(1, 4)] =  corr([mfrRat2(:);mfrRat1(:)], val_rat_in(:), 'type', typeCorr);
    
    [CorrRat(1, 5), pvalRat(1, 5)] =  corr([mfrRat1(:);mfrRat2(:)], abs(val_rat_out(:)), 'type', typeCorr);
    [CorrRat(1, 6), pvalRat(1, 6)] =  corr([mfrRat2(:);mfrRat1(:)], abs(val_rat_out(:)), 'type', typeCorr);
    [CorrRat(1, 7), pvalRat(1, 7)] =  corr([mfrRat1(:);mfrRat2(:)], abs(val_rat_in(:)), 'type', typeCorr);
    [CorrRat(1, 8), pvalRat(1, 8)] =  corr([mfrRat2(:);mfrRat1(:)], abs(val_rat_in(:)), 'type', typeCorr);
    




function  [Corr2AFC, pval2AFC]= find_corr_2AFC (mfr1, mfr2, val_2AFCout, val_A, val_B, val_C, val_U, typeCorr)

% val_out is always the answer from the patient 
% val A and B is the input, so the expected value from the stimulus, and can be
% coded as as 4 types, ranting, ranking, AFC and total
% val C is the iexpected value (input) of the chosen stimulus, and U of the
% unchosen one


Corr2AFC = nan(1, 30);
pval2AFC = nan(1, 30);

    
    
    % corr output value, so the revealed preferences
    [Corr2AFC(1, 1), pval2AFC(1, 1)] = corr([mfr1(:);mfr2(:)], val_2AFCout(:), 'type', typeCorr);
    [Corr2AFC(1, 2), pval2AFC(1, 2)] = corr([mfr2(:);mfr1(:)], val_2AFCout(:), 'type', typeCorr);
    
    % corr input value during A
    [Corr2AFC(1, 3), pval2AFC(1, 3)] = corr([mfr1(:);mfr2(:)], val_A(:), 'type', typeCorr);
    [Corr2AFC(1, 4), pval2AFC(1, 4)] = corr([mfr2(:);mfr1(:)], val_A(:), 'type', typeCorr);
    
    % corr input value during B
    [Corr2AFC(1, 5), pval2AFC(1, 5)] = corr([mfr1(:);mfr2(:)], val_B(:), 'type', typeCorr);
    [Corr2AFC(1, 6), pval2AFC(1, 6)] = corr([mfr2(:);mfr1(:)], val_B(:), 'type', typeCorr);
    
    % corr absolute preference
    [Corr2AFC(1, 7), pval2AFC(1, 7)] = corr([mfr1(:);mfr2(:)], abs(val_2AFCout(:)), 'type', typeCorr);
    [Corr2AFC(1, 8), pval2AFC(1, 8)] = corr([mfr2(:);mfr1(:)], abs(val_2AFCout(:)), 'type', typeCorr);
    
    % corr salience of A
    [Corr2AFC(1, 9), pval2AFC(1, 9)] = corr([mfr1(:);mfr2(:)], abs(val_A(:)), 'type', typeCorr);
    [Corr2AFC(1, 10), pval2AFC(1, 10)] = corr([mfr2(:);mfr1(:)], abs(val_A(:)), 'type', typeCorr);
    
    % corr salience of B
    [Corr2AFC(1, 11), pval2AFC(1, 11)] = corr([mfr1(:);mfr2(:)], abs(val_B(:)), 'type', typeCorr);
    [Corr2AFC(1, 12), pval2AFC(1, 12)] = corr([mfr2(:);mfr1(:)], abs(val_B(:)), 'type', typeCorr);
    
    % corr difference value A-B
    [Corr2AFC(1, 13), pval2AFC(1, 13)] = corr([mfr1(:);mfr2(:)], val_A(:) - val_B(:), 'type', typeCorr);
    [Corr2AFC(1, 14), pval2AFC(1, 14)] = corr([mfr2(:);mfr1(:)], val_A(:) - val_B(:), 'type', typeCorr);

    % corr absolute difference value A-B
    [Corr2AFC(1, 15), pval2AFC(1, 15)] = corr([mfr1(:);mfr2(:)], abs(val_B(:) - val_A(:)), 'type', typeCorr);
    [Corr2AFC(1, 16), pval2AFC(1, 16)] = corr([mfr2(:);mfr1(:)], abs(val_B(:) - val_A(:)), 'type', typeCorr);
    
    % corr sum A-B
    [Corr2AFC(1, 17), pval2AFC(1, 17)] = corr([mfr1(:);mfr2(:)], (val_B(:) + val_A(:)), 'type', typeCorr);
    [Corr2AFC(1, 18), pval2AFC(1, 18)] = corr([mfr2(:);mfr1(:)], (val_B(:) + val_A(:)), 'type', typeCorr);
    
    % corr sum of salience A-B
    [Corr2AFC(1, 19), pval2AFC(1, 19)] = corr([mfr1(:);mfr2(:)], abs(val_B(:)) + abs(val_A(:)), 'type', typeCorr);
    [Corr2AFC(1, 20), pval2AFC(1, 20)] = corr([mfr2(:);mfr1(:)], abs(val_B(:)) + abs(val_A(:)), 'type', typeCorr);
    
    % corr absolute of sum A-B
    [Corr2AFC(1, 21), pval2AFC(1, 21)] = corr([mfr1(:);mfr2(:)], abs(val_B(:) + val_A(:)), 'type', typeCorr);
    [Corr2AFC(1, 22), pval2AFC(1, 22)] = corr([mfr2(:);mfr1(:)], abs(val_B(:) + val_A(:)), 'type', typeCorr);
    
    % corr chosen value 
    [Corr2AFC(1, 23), pval2AFC(1, 23)] = corr([mfr1(:);mfr2(:)], val_C(:), 'type', typeCorr);
    [Corr2AFC(1, 24), pval2AFC(1, 24)] = corr([mfr2(:);mfr1(:)], val_C(:), 'type', typeCorr);
    
    % corr unchosen value
    [Corr2AFC(1, 25), pval2AFC(1, 25)] = corr([mfr1(:);mfr2(:)], val_U(:), 'type', typeCorr);
    [Corr2AFC(1, 26), pval2AFC(1, 26)] = corr([mfr2(:);mfr1(:)], val_U(:), 'type', typeCorr);
    
    % corr chosen salience
    [Corr2AFC(1, 27), pval2AFC(1, 27)] = corr([mfr1(:);mfr2(:)], abs(val_C(:)), 'type', typeCorr);
    [Corr2AFC(1, 28), pval2AFC(1, 28)] = corr([mfr2(:);mfr1(:)], abs(val_C(:)), 'type', typeCorr);
    
    % corr unchosen salience    
    [Corr2AFC(1, 29), pval2AFC(1, 29)] = corr([mfr1(:);mfr2(:)], abs(val_U(:)), 'type', typeCorr);
    [Corr2AFC(1, 30), pval2AFC(1, 30)] = corr([mfr2(:);mfr1(:)], abs(val_U(:)), 'type', typeCorr);
    



function [PanovaRat]= find_ancova_rating (mfrRat1, mfrRat2, val_rat_out, val_rat_in)

% val out is always the answer from the patient 
% val in is the input, so the expected value from the stimulus, and can be
% coded as as 4 types, ranting, ranking, AFC and total


PanovaRat = nan(1, 8);
TanovaRat = cell(1, 8);

    
    [PanovaRat(1, 1), TanovaRat{1, 1}] =  anovan([mfrRat1(:);mfrRat2(:)], val_rat_out(:), 'continuous',[1],'varnames',{'val_rat_out'},'display','off');
    [PanovaRat(1, 2), TanovaRat{1, 2}] =  anovan([mfrRat2(:);mfrRat1(:)], val_rat_out(:), 'continuous',[1],'varnames',{'val_rat_out'},'display','off');
    [PanovaRat(1, 3), TanovaRat{1, 3}] =  anovan([mfrRat1(:);mfrRat2(:)], val_rat_in(:), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat(1, 4), TanovaRat{1, 4}] =  anovan([mfrRat2(:);mfrRat1(:)], val_rat_in(:), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    
    [PanovaRat(1, 5), TanovaRat{1, 5}] =  anovan([mfrRat1(:);mfrRat2(:)], abs(val_rat_out(:)), 'continuous',[1],'varnames',{'sal_rat_out'},'display','off');
    [PanovaRat(1, 6), TanovaRat{1, 6}] =  anovan([mfrRat2(:);mfrRat1(:)], abs(val_rat_out(:)), 'continuous',[1],'varnames',{'sal_rat_out'},'display','off');
    [PanovaRat(1, 7), TanovaRat{1, 7}] =  anovan([mfrRat1(:);mfrRat2(:)], abs(val_rat_in(:)), 'continuous',[1],'varnames',{'sal_rat_in'},'display','off');
    [PanovaRat(1, 8), TanovaRat{1, 8}] =  anovan([mfrRat2(:);mfrRat1(:)], abs(val_rat_in(:)), 'continuous',[1],'varnames',{'sal_rat_in'},'display','off');
    




function  [Panova2AFC]= find_ancova_2AFC (mfr1, mfr2, val_2AFCout, val_A, val_B, val_C, val_U)

% val_out is always the answer from the patient 
% val A and B is the input, so the expected value from the stimulus, and can be
% coded as as 4 types, ranting, ranking, AFC and total
% val C is the iexpected value (input) of the chosen stimulus, and U of the
% unchosen one

Panova2AFC = nan(1, 36);
Tanova2AFC = cell(1, 36);


    
    % corr output value, so the revealed preferences
    [Panova2AFC(1, 1), Tanova2AFC{1, 1}] = anovan([mfr1(:);mfr2(:)], val_2AFCout(:), 'continuous',[1],'varnames',{'val_2AFCout'},'display','off');
    [Panova2AFC(1, 2), Tanova2AFC{1, 2}] = anovan([mfr2(:);mfr1(:)], val_2AFCout(:),'continuous',[1],'varnames',{'val_2AFCout'},'display','off');
    
    % corr input value during A
    [Panova2AFC(1, 3), Tanova2AFC{1, 3}] = anovan([mfr1(:);mfr2(:)], val_A(:), 'continuous',[1],'varnames',{'val_A'},'display','off');
    [Panova2AFC(1, 4), Tanova2AFC{1, 4}] = anovan([mfr2(:);mfr1(:)], val_A(:), 'continuous',[1],'varnames',{'val_A'},'display','off');
    
    % corr input value during B
    [Panova2AFC(1, 5), Tanova2AFC{1, 5}] = anovan([mfr1(:);mfr2(:)], val_B(:), 'continuous',[1],'varnames',{'val_B'},'display','off');
    [Panova2AFC(1, 6), Tanova2AFC{1, 6}] = anovan([mfr2(:);mfr1(:)], val_B(:), 'continuous',[1],'varnames',{'val_B'},'display','off');
    
    % corr absolute preference
    [Panova2AFC(1, 7), Tanova2AFC{1, 7}] = anovan([mfr1(:);mfr2(:)], abs(val_2AFCout(:)), 'continuous',[1],'varnames',{'val_2AFCout'},'display','off');
    [Panova2AFC(1, 8), Tanova2AFC{1, 8}] = anovan([mfr2(:);mfr1(:)], abs(val_2AFCout(:)), 'continuous',[1],'varnames',{'val_2AFCout'},'display','off');
    
    % corr salience of A
    [Panova2AFC(1, 9), Tanova2AFC{1, 9}] = anovan([mfr1(:);mfr2(:)], abs(val_A(:)), 'continuous',[1],'varnames',{'sal_A'},'display','off');
    [Panova2AFC(1, 10), Tanova2AFC{1, 10}] = anovan([mfr2(:);mfr1(:)], abs(val_A(:)), 'continuous',[1],'varnames',{'sal_A'},'display','off');
    
    % corr salience of B
    [Panova2AFC(1, 11), Tanova2AFC{1, 11}] = anovan([mfr1(:);mfr2(:)], abs(val_B(:)), 'continuous',[1],'varnames',{'sal_B'},'display','off');
    [Panova2AFC(1, 12), Tanova2AFC{1, 12}] = anovan([mfr2(:);mfr1(:)], abs(val_B(:)), 'continuous',[1],'varnames',{'sal_B'},'display','off');
    
    % corr difference value A-B
    [Panova2AFC(1, 13), Tanova2AFC{1, 13}] = anovan([mfr1(:);mfr2(:)], val_A(:) - val_B(:), 'continuous',[1],'varnames',{'diff'},'display','off');
    [Panova2AFC(1, 14), Tanova2AFC{1, 14}] = anovan([mfr2(:);mfr1(:)], val_A(:) - val_B(:), 'continuous',[1],'varnames',{'diff'},'display','off');

    % corr absolute difference value A-B
    [Panova2AFC(1, 15), Tanova2AFC{1, 15}] = anovan([mfr1(:);mfr2(:)], abs(val_B(:) - val_A(:)), 'continuous',[1],'varnames',{'abs_diff'},'display','off');
    [Panova2AFC(1, 16), Tanova2AFC{1, 16}] = anovan([mfr2(:);mfr1(:)], abs(val_B(:) - val_A(:)), 'continuous',[1],'varnames',{'abs_diff'},'display','off');
    
    % corr sum A-B
    [Panova2AFC(1, 17), Tanova2AFC{1, 17}] = anovan([mfr1(:);mfr2(:)], (val_B(:) + val_A(:)), 'continuous',[1],'varnames',{'sum'},'display','off');
    [Panova2AFC(1, 18), Tanova2AFC{1, 18}] = anovan([mfr2(:);mfr1(:)], (val_B(:) + val_A(:)), 'continuous',[1],'varnames',{'sum'},'display','off');
    
    % corr sum of salience A-B
    [Panova2AFC(1, 19), Tanova2AFC{1, 19}] = anovan([mfr1(:);mfr2(:)], abs(val_B(:)) + abs(val_A(:)), 'continuous',[1],'varnames',{'abs_sum'},'display','off');
    [Panova2AFC(1, 20), Tanova2AFC{1, 20}] = anovan([mfr2(:);mfr1(:)], abs(val_B(:)) + abs(val_A(:)), 'continuous',[1],'varnames',{'abs_sum'},'display','off');
    
    % corr absolute of sum A-B
    [Panova2AFC(1, 21), Tanova2AFC{1, 21}] = anovan([mfr1(:);mfr2(:)], abs(val_B(:) + val_A(:)), 'continuous',[1],'varnames',{'abs_abs_sum'},'display','off');
    [Panova2AFC(1, 22), Tanova2AFC{1, 22}] = anovan([mfr2(:);mfr1(:)], abs(val_B(:) + val_A(:)), 'continuous',[1],'varnames',{'sal_rat_in'},'display','off');
    
    % corr chosen value 
    [Panova2AFC(1, 23), Tanova2AFC{1, 23}] = anovan([mfr1(:);mfr2(:)], val_C(:), 'continuous',[1],'varnames',{'val_C'},'display','off');
    [Panova2AFC(1, 24), Tanova2AFC{1, 24}] = anovan([mfr2(:);mfr1(:)], val_C(:), 'continuous',[1],'varnames',{'val_C'},'display','off');
    
    % corr unchosen value
    [Panova2AFC(1, 25), Tanova2AFC{1, 25}] = anovan([mfr1(:);mfr2(:)], val_U(:), 'continuous',[1],'varnames',{'val_U'},'display','off');
    [Panova2AFC(1, 26), Tanova2AFC{1, 26}] = anovan([mfr2(:);mfr1(:)], val_U(:), 'continuous',[1],'varnames',{'val_U'},'display','off');
    
    % corr chosen salience
    [Panova2AFC(1, 27), Tanova2AFC{1, 27}] = anovan([mfr1(:);mfr2(:)], abs(val_C(:)), 'continuous',[1],'varnames',{'sal_C'},'display','off');
    [Panova2AFC(1, 28), Tanova2AFC{1, 28}] = anovan([mfr2(:);mfr1(:)], abs(val_C(:)), 'continuous',[1],'varnames',{'sal_C'},'display','off');
    
    % corr unchosen salience    
    [Panova2AFC(1, 29), Tanova2AFC{1, 29}] = anovan([mfr1(:);mfr2(:)], abs(val_U(:)), 'continuous',[1],'varnames',{'sal_U'},'display','off');
    [Panova2AFC(1, 30), Tanova2AFC{1, 30}] = anovan([mfr2(:);mfr1(:)], abs(val_U(:)), 'continuous',[1],'varnames',{'sal_U'},'display','off');
    
    
   
function  [Panova2way]= find_2wayancova_2AFC (mfr1, mfr2, val_2AFCout, val_A, val_B) 
    
Panova2way = nan(1, 36);


%% 2 way ancovas

     % corr output value, so the revealed preferences Hill
     [Panova2way(1,  1:3)] = anovan([mfr1(:);mfr2(:)], {val_2AFCout(:), val_A(:), val_B(:)}, 'continuous',[1,2,3],'varnames',{'val_2AFCout', 'val_A', 'val_B'},'display','off');
     [Panova2way(1,  4:6)] = anovan([mfr2(:);mfr1(:)], {val_2AFCout(:), val_A(:), val_B(:)}, 'continuous',[1,2,3],'varnames',{'val_2AFCout','val_A', 'val_B'},'display','off');
    
     % corr output value, so the revealed preferences Hill
     [Panova2way(1,  7:8)] = anovan([mfr1(:);mfr2(:)], {val_A(:), val_B(:)}, 'continuous',[1,2],'varnames',{'val_A', 'val_B' },'display','off');
     [Panova2way(1, 9:10)] = anovan([mfr2(:);mfr1(:)], {val_A(:), val_B(:)},'continuous',[1,2],'varnames',{'val_A', 'val_B'},'display','off');
     
     
     % corr output value, so the revealed preferences Hill
     [Panova2way(1, 11:13)] = anovan([mfr1(:);mfr2(:)], {val_A(:), val_B(:)}, 'continuous',[1,2],'varnames',{'val_A', 'val_B'}, 'display','off', 'model', 'interaction');
     [Panova2way(1, 14:16)] = anovan([mfr2(:);mfr1(:)], {val_A(:), val_B(:)}, 'continuous',[1,2],'varnames',{'val_A', 'val_B'}, 'display','off', 'model', 'interaction');
     
     % corr absolute output value, so the absolute revealed preferences
     % form DeMartino analysis
     [Panova2way(1, 17:19)] = anovan([mfr1(:);mfr2(:)], {abs(val_2AFCout(:)), val_A(:), val_B(:)}, 'continuous',[1,2,3],'varnames',{'abs_val2AFCout', 'val_A', 'val_B'},'display','off');
     [Panova2way(1, 20:22)] = anovan([mfr2(:);mfr1(:)], {abs(val_2AFCout(:)), val_A(:), val_B(:)}, 'continuous',[1,2,3],'varnames',{'abs_val2AFCout','val_A', 'val_B'},'display','off');    
     
     % corr value difference, controlled by A 
     [Panova2way(1, 23:24)] = anovan([mfr1(:);mfr2(:)], {(val_A(:)-val_B(:)), val_A(:)}, 'continuous',[1,2],'varnames',{'val_diff', 'val_A'},'display','off');
     [Panova2way(1, 25:26)] = anovan([mfr2(:);mfr1(:)], {(val_A(:)-val_B(:)), val_A(:)}, 'continuous',[1,2],'varnames',{'val_diff','val_A'},'display','off');    
     
     % corr value difference, controlled by B 
     [Panova2way(1, 27:28)] = anovan([mfr1(:);mfr2(:)], {(val_A(:)-val_B(:)), val_B(:)}, 'continuous',[1,2],'varnames',{'val_diff','val_B'},'display','off');
     [Panova2way(1, 29:30)] = anovan([mfr2(:);mfr1(:)], {(val_A(:)-val_B(:)), val_B(:)}, 'continuous',[1,2],'varnames',{'val_diff','val_B'},'display','off');    
     
     % corr value difference, controlled by A and B 
     [Panova2way(1, 31:33)] = anovan([mfr1(:);mfr2(:)], {(val_A(:)-val_B(:)), val_A(:), val_B(:)}, 'continuous',[1,2,3],'varnames',{'val_diff', 'val_A','val_B'},'display','off');
     [Panova2way(1, 34:36)] = anovan([mfr2(:);mfr1(:)], {(val_A(:)-val_B(:)), val_A(:), val_B(:)}, 'continuous',[1,2,3],'varnames',{'val_diff', 'val_A','val_B'},'display','off');    
        
     
     
function [Panova, Pkruskal]= anova_kruskall_detection (stimuliRat, stimuliA, stimuliB, mfrRat1, mfrRat2, mfrlikRat1, mfrlikRat2, mfrA11, mfrA12, mfrB1, mfrB2, mfrA21, mfrA22, mfrlikAFC1, mfrlikAFC2)

[Panova(1,1:2),  Pkruskal(1,1:2)]=  anova_kruskall_rat (mfrRat1, mfrRat2, stimuliRat);
[Panova(1,3:4),  Pkruskal(1,3:4)]=  anova_kruskall_rat (mfrlikRat1, mfrlikRat2, stimuliRat);
[Panova(1,5:8),  Pkruskal(1,5:8)]=  anova_kruskall_AFC (mfrA11, mfrA12,  stimuliA, stimuliB);
[Panova(1,9:12), Pkruskal(1,9:12)]=  anova_kruskall_AFC (mfrB1, mfrB2,  stimuliA, stimuliB);
[Panova(1,13:16),Pkruskal(1,13:16)]= anova_kruskall_AFC (mfrA21, mfrA22,  stimuliA, stimuliB);
[Panova(1,17:20),Pkruskal(1,17:20)]= anova_kruskall_AFC (mfrlikAFC1, mfrlikAFC2, stimuliA, stimuliB);



function  [Panova, Pkruskal]= anova_kruskall_rat (mfr1, mfr2, stimuli)


Panova(1,1) = anova1 (mfr1(:), stimuli(:,1), 'off');
Panova(1,2) = anova1 (mfr2(:), stimuli(:,2), 'off');
Pkruskal(1,1)=kruskalwallis(mfr1(:), stimuli(:,1), 'off');
Pkruskal(1,2)=kruskalwallis(mfr2(:), stimuli(:,2), 'off');



function  [Panova, Pkruskal]= anova_kruskall_AFC (mfr1, mfr2, stimuliA, stimuliB)


Panova(1,1) = anova1 (mfr1(:), stimuliA(:,1), 'off');
Panova(1,2) = anova1 (mfr2(:), stimuliA(:,2), 'off');
Panova(1,3) = anova1 (mfr1(:), stimuliB(:,1), 'off');
Panova(1,4) = anova1 (mfr2(:), stimuliB(:,2), 'off');
Pkruskal(1,1) = kruskalwallis(mfr1(:), stimuliA(:,1), 'off');
Pkruskal(1,2) = kruskalwallis(mfr2(:), stimuliA(:,2), 'off');
Pkruskal(1,3) = kruskalwallis(mfr1(:), stimuliB(:,1), 'off');
Pkruskal(1,4) = kruskalwallis(mfr2(:), stimuliB(:,2), 'off');



function [CorrRat, pvalRat]=change_valueRatingCorr(mfr1, mfr2, val_rat_in, typeCorr)

CorrRat = nan(1, 4);
pvalRat = nan(1, 4);

    
% correlation change vs mfr during rating
[CorrRat(1, 1), pvalRat(1, 1)] = corr(mfr1(:), val_rat_in(:,1)- val_rat_in(:,2), 'type', typeCorr);
[CorrRat(1, 2), pvalRat(1, 2)] = corr(mfr2(:), val_rat_in(:,1)- val_rat_in(:,2), 'type', typeCorr);

% correlation absolute change vs mfr during rating
[CorrRat(1, 3), pvalRat(1, 3)] = corr(mfr1(:), abs(val_rat_in(:,1)- val_rat_in(:,2)), 'type', typeCorr);
[CorrRat(1, 4), pvalRat(1, 4)] = corr(mfr2(:), abs(val_rat_in(:,1)- val_rat_in(:,2)), 'type', typeCorr);




function [Corr2AFC, pval2AFC]= change_valueAFCCorr(mfr1, mfr2, val_A, val_B, typeCorr)

Corr2AFC = nan(1, 8);
pval2AFC = nan(1, 8);

    
% correlation change vs mfr during A
[Corr2AFC(1, 1), pval2AFC(1, 1)] = corr(mfr1(:), val_A(:,1)- val_A(:,2), 'type', typeCorr);
[Corr2AFC(1, 2), pval2AFC(1, 2)] = corr(mfr2(:), val_A(:,1)- val_A(:,2), 'type', typeCorr);

% correlation absolute change vs mfr during A
[Corr2AFC(1, 3), pval2AFC(1, 3)] = corr(mfr1(:), abs(val_A(:,1)- val_A(:,2)), 'type', typeCorr);
[Corr2AFC(1, 4), pval2AFC(1, 4)] = corr(mfr2(:), abs(val_A(:,1)- val_A(:,2)), 'type', typeCorr);


% correlation change vs mfr during B
[Corr2AFC(1, 5), pval2AFC(1, 5)] = corr(mfr1(:), val_B(:,1)- val_B(:,2), 'type', typeCorr);
[Corr2AFC(1, 6), pval2AFC(1, 6)] = corr(mfr2(:), val_B(:,1)- val_B(:,2), 'type', typeCorr);

% correlation absolute change vs mfr during B
[Corr2AFC(1, 7), pval2AFC(1, 7)] = corr(mfr1(:), abs(val_B(:,1)- val_B(:,2)), 'type', typeCorr);
[Corr2AFC(1, 8), pval2AFC(1, 8)] = corr(mfr2(:), abs(val_B(:,1)- val_B(:,2)), 'type', typeCorr);




function [PancovaRat]=change_valueRatingAncova(mfr1, mfr2, val_rat_in)

PancovaRat = nan(1, 4);


% correlation change vs mfr during rating
[PancovaRat(1, 1)] = anovan(mfr1(:), val_rat_in(:,1)- val_rat_in(:,2), 'continuous',[1],'display','off');
[PancovaRat(1, 2)] = anovan(mfr2(:), val_rat_in(:,1)- val_rat_in(:,2), 'continuous',[1],'display','off');

% correlation absolute change vs mfr during rating
[PancovaRat(1, 3)] = anovan(mfr1(:), abs(val_rat_in(:,1)- val_rat_in(:,2)), 'continuous',[1],'display','off');
[PancovaRat(1, 4)] = anovan(mfr2(:), abs(val_rat_in(:,1)- val_rat_in(:,2)), 'continuous',[1],'display','off');




function [Pancova2AFC]= change_valueAFCAncova(mfr1, mfr2, val_A, val_B)
 

Pancova2AFC = nan(1, 8);

    
% correlation change vs mfr during A
[Pancova2AFC(1, 1)] = anovan(mfr1(:), val_A(:,1)- val_A(:,2), 'continuous',[1],'display','off');
[Pancova2AFC(1, 2)] = anovan(mfr2(:), val_A(:,1)- val_A(:,2), 'continuous',[1],'display','off');

% correlation absolute change vs mfr during A
[Pancova2AFC(1, 3)] = anovan(mfr1(:), abs(val_A(:,1)- val_A(:,2)), 'continuous',[1],'display','off');
[Pancova2AFC(1, 4)] = anovan(mfr2(:), abs(val_A(:,1)- val_A(:,2)), 'continuous',[1],'display','off');


% correlation change vs mfr during B
[Pancova2AFC(1, 5)] = anovan(mfr1(:), val_B(:,1)- val_B(:,2), 'continuous',[1],'display','off');
[Pancova2AFC(1, 6)] = anovan(mfr2(:), val_B(:,1)- val_B(:,2), 'continuous',[1],'display','off');

% correlation absolute change vs mfr during B
[Pancova2AFC(1, 7)] = anovan(mfr1(:), abs(val_B(:,1)- val_B(:,2)), 'continuous',[1],'display','off');
[Pancova2AFC(1, 8)] = anovan(mfr2(:), abs(val_B(:,1)- val_B(:,2)), 'continuous',[1],'display','off');

