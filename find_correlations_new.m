function find_correlations_new

tic
dbstop if error

%% settings 


rsp_t_AFC= [0,1000]; %use at the beginning this time period to find the corrlations, this could change
rsp_t_likertAFC= [0,1000]; %use at the beginning this time period to find the corrlations, this could change
rsp_t_rating= [0,1000]; %use at the beginning this time period to find the corrlations, this could change
rsp_t_likertRat= [0,1000]; %use at the beginning this time period to find the corrlations, this could change
rsp_t_bsl= [-400,0];

value_signal='Rank'; % this is the value signal used in this analysis, thus, Tot, Rat, Rank or AFC.  If 'All', then it loop over the four types of signlas and saves all of them 
typeCorr= 'Spearman'; % here choose which correlation to use 



%% 

anatitle = inputdlg('Please enter a meaningful title (file name) for this analysis','Analysis title',1,{''});
ana_title= char(anatitle{1,1});
mkdir(sprintf('/media/Projects/Alex/New Analysis/%s',ana_title));
folder_to_save = sprintf('/media/Projects/Alex/New Analysis/%s',ana_title);

 

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



if isequal(value_signal, 'All') 
    value_signal=[{'Rank'}, {'Rat'}, {'Tot'}, {'AFC'}];
    
    for hh=1:length(value_signal)
        main(value_signal(hh), rasts, rsp_t_AFC, rsp_t_rating, rsp_t_likertRat, rsp_t_likertAFC, rsp_t_bsl, folder_to_save, typeCorr, namefiles, outfiles, behavDir)
    end
    
else 
    main(value_signal, rasts, rsp_t_AFC, rsp_t_rating, rsp_t_likertRat, rsp_t_likertAFC, rsp_t_bsl, folder_to_save, typeCorr, namefiles, outfiles, behavDir)

end





function main(value_signal, rasts, rsp_t_AFC, rsp_t_rating, rsp_t_likertRat, rsp_t_likertAFC, rsp_t_bsl, folder_to_save, typeCorr, namefiles, outfiles, behavDir)

if length(namefiles) == length(outfiles)
    
    n_ses_units=size(rasts,1)/23;
    CorrRat= nan(n_ses_units, 8*2);
    CorrLikRat= nan(n_ses_units, 8*2);
    CorrA1= nan(n_ses_units, 30*2);
    CorrB= nan(n_ses_units, 30*2);
    CorrA2= nan(n_ses_units, 30*2);
    CorrLikAFC= nan(n_ses_units, 30*2);
    CorrA1_B= nan(n_ses_units, 30*2);
    CorrB_A2= nan(n_ses_units, 30*2);
    CorrA1_Babs= nan(n_ses_units, 30*2);
    CorrB_A2abs= nan(n_ses_units, 30*2);
   
    pvalRat= nan(n_ses_units, 8*2);
    pvalLikRat= nan(n_ses_units, 8*2);
    pvalA1= nan(n_ses_units, 30*2);
    pvalB= nan(n_ses_units, 30*2);
    pvalA2= nan(n_ses_units, 30*2);
    pvalLikAFC= nan(n_ses_units, 30*2);
    pvalA1_B= nan(n_ses_units, 30*2);
    pvalB_A2= nan(n_ses_units, 30*2);
    pvalA1_Babs= nan(n_ses_units, 30*2);
    pvalB_A2abs= nan(n_ses_units, 30*2);
    
    CorrRatbsl= nan(n_ses_units, 8*2);
    CorrLikRatbsl= nan(n_ses_units, 8*2);
    CorrA1bsl= nan(n_ses_units, 30*2);
    CorrBbsl= nan(n_ses_units, 30*2);
    CorrA2bsl= nan(n_ses_units, 30*2);
    CorrLikAFCbsl= nan(n_ses_units, 30*2);
    CorrA1_Bbsl= nan(n_ses_units, 30*2);
    CorrB_A2bsl= nan(n_ses_units, 30*2);
    CorrA1_Babsbsl= nan(n_ses_units, 30*2);
    CorrB_A2absbsl= nan(n_ses_units, 30*2);
   
    pvalRatbsl= nan(n_ses_units, 8*2);
    pvalLikRatbsl= nan(n_ses_units, 8*2);
    pvalA1bsl= nan(n_ses_units, 30*2);
    pvalBbsl= nan(n_ses_units, 30*2);
    pvalA2bsl= nan(n_ses_units, 30*2);
    pvalLikAFCbsl= nan(n_ses_units, 30*2);
    pvalA1_Bbsl= nan(n_ses_units, 30*2);
    pvalB_A2bsl= nan(n_ses_units, 30*2);
    pvalA1_Babsbsl= nan(n_ses_units, 30*2);
    pvalB_A2absbsl= nan(n_ses_units, 30*2);
    
    

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
        unitind=ismember(cell2mat(rasts(:,2:5)),cell2mat(outfiles(sess,2:5)),'rows');


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
            mfrRat2bsl= mfrRat2 - bslRat2;
            mfrlikRat1bsl =mfrlikRat1 - bslRat1;
            mfrlikRat2bsl =mfrlikRat2 - bslRat2;
            
            
            

                
      %%
            
            
            
            % select the values to use
            
            eval(['val_A = val_A_', cell2mat(value_signal)]);
            eval(['val_B = val_B_', cell2mat(value_signal)]);
            eval(['val_C = val_C_', cell2mat(value_signal)]);
            eval(['val_U = val_U_', cell2mat(value_signal)]);
            eval(['val_rat_afc = val_rat_', cell2mat(value_signal)]);

            
            
            % make correlation matrices for each trial (190 for AFC, and 60 for rating)
            
            [CorrRat(unitnum, :), pvalRat(unitnum, :)] = find_corr_rating (mfrRat1, mfrRat2, val_rat, val_rat_afc, typeCorr);
            [CorrLikRat(unitnum, :), pvalLikRat(unitnum, :)] = find_corr_rating (mfrlikRat1, mfrlikRat2, val_rat, val_rat_afc, typeCorr);
            [CorrA1(unitnum, :), pvalA1(unitnum, :)] = find_corr_2AFC  (mfrA11, mfrA12, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrB(unitnum, :), pvalB(unitnum, :)]   = find_corr_2AFC  (mfrB1,  mfrB2, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrA2(unitnum, :), pvalA2(unitnum, :)]  = find_corr_2AFC  (mfrA21, mfrA22, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrLikAFC(unitnum, :), pvalLikAFC(unitnum, :)] = find_corr_2AFC  (mfrlikAFC1, mfrlikAFC2, val_2AFC, val_A, val_B, val_C_Rank, val_U, typeCorr);
            
            [CorrA1_B(unitnum, :), pvalA1_B(unitnum, :)] = find_corr_2AFC  (mfrA11-mfrB1, mfrA12-mfrB2, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrB_A2(unitnum, :), pvalB_A2(unitnum, :)] = find_corr_2AFC  (mfrB1-mfrA21, mfrB2-mfrA22, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrA1_Babs(unitnum, :), pvalA1_Babs(unitnum, :)]  = find_corr_2AFC  (abs(mfrA11-mfrB1), abs(mfrA12-mfrB2), val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrB_A2abs(unitnum, :), pvalB_A2abs(unitnum, :) ] = find_corr_2AFC  (abs(mfrB1-mfrA21), abs(mfrB2-mfrA22), val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            
            
            
           
            % baseline subtraction
            [CorrRatbsl(unitnum, :), pvalRatbsl(unitnum, :)] = find_corr_rating (mfrRat1bsl, mfrRat2bsl, val_rat, val_rat_afc, typeCorr);
            [CorrLikRatbsl(unitnum, :), pvalLikRatbsl(unitnum, :)] = find_corr_rating (mfrlikRat1bsl, mfrlikRat2bsl, val_rat, val_rat_afc, typeCorr);
            [CorrA1bsl(unitnum, :), pvalA1bsl(unitnum, :)] = find_corr_2AFC  (mfrA11bsl, mfrA12bsl, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrBbsl(unitnum, :), pvalBbsl(unitnum, :)]   = find_corr_2AFC  (mfrB1bsl,  mfrB2bsl, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrA2bsl(unitnum, :), pvalA2bsl(unitnum, :)]  = find_corr_2AFC  (mfrA21bsl, mfrA22bsl, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrLikAFCbsl(unitnum, :), pvalLikAFCbsl(unitnum, :)] = find_corr_2AFC  (mfrlikAFC1bsl, mfrlikAFC2bsl, val_2AFC, val_A, val_B, val_C_Rank, val_U, typeCorr);
            
            [CorrA1_Bbsl(unitnum, :), pvalA1_Bbsl(unitnum, :)] = find_corr_2AFC  (mfrA11bsl-mfrB1bsl, mfrA12bsl-mfrB2bsl, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrB_A2bsl(unitnum, :), pvalB_A2bsl(unitnum, :)] = find_corr_2AFC  (mfrB1bsl-mfrA21bsl, mfrB2bsl-mfrA22bsl, val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrA1_Babsbsl(unitnum, :), pvalA1_Babsbsl(unitnum, :)]  = find_corr_2AFC  (abs(mfrA11bsl-mfrB1bsl), abs(mfrA12bsl-mfrB2bsl), val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            [CorrB_A2absbsl(unitnum, :), pvalB_A2absbsl(unitnum, :) ] = find_corr_2AFC  (abs(mfrB1bsl-mfrA21bsl), abs(mfrB2bsl-mfrA22bsl), val_2AFC, val_A, val_B, val_C, val_U, typeCorr);
            

            
            %% ancova analysis
            
            
            [PanovaRat{unitnum,1}, TanovaRat{unitnum,1}] = find_ancova_rating (mfrRat1, mfrRat2, val_rat, val_rat_afc);
            [PanovaLikRat{unitnum,1}, TanovaLikRat{unitnum,1}] = find_ancova_rating (mfrlikRat1, mfrlikRat2, val_rat, val_rat_afc);
            [PanovaA1{unitnum,1}, TanovaA1{unitnum,1}] = find_ancova_2AFC  (mfrA11, mfrA12, val_2AFC, val_A, val_B, val_C, val_U);
            [PanovaB{unitnum,1}, TanovaB{unitnum,1}]   = find_ancova_2AFC  (mfrB1,  mfrB2, val_2AFC, val_A, val_B, val_C, val_U);
            [PanovaA2{unitnum,1}, TanovaA2{unitnum,1}]  = find_ancova_2AFC  (mfrA21, mfrA22, val_2AFC, val_A, val_B, val_C, val_U);
            [PanovaLikAFC{unitnum,1}, TanovaLikAFC{unitnum,1}] = find_ancova_2AFC  (mfrlikAFC1, mfrlikAFC2, val_2AFC, val_A, val_B, val_C_Rank, val_U);
            
            [PanovaA1_B{unitnum,1}, TanovaA1_B{unitnum,1}] = find_ancova_2AFC  (mfrA11-mfrB1, mfrA12-mfrB2, val_2AFC, val_A, val_B, val_C, val_U);
            [PanovaB_A2{unitnum,1}, TanovaB_A2{unitnum,1}] = find_ancova_2AFC  (mfrB1-mfrA21, mfrB2-mfrA22, val_2AFC, val_A, val_B, val_C, val_U);
            [PanovaA1_Babs{unitnum,1}, TanovaA1_Babs{unitnum,1}]  = find_ancova_2AFC  (abs(mfrA11-mfrB1), abs(mfrA12-mfrB2), val_2AFC, val_A, val_B, val_C, val_U);
            [PanovaB_A2abs{unitnum,1}, TanovaB_A2abs{unitnum,1}] = find_ancova_2AFC  (abs(mfrB1-mfrA21), abs(mfrB2-mfrA22), val_2AFC, val_A, val_B, val_C, val_U);
            
            
            
           
            % baseline subtraction
            [PanovaRatbsl{unitnum,1}, TanovaRatbsl{unitnum,1}] = find_ancova_rating (mfrRat1bsl, mfrRat2bsl, val_rat, val_rat_afc);
            [PanovaLikRatbsl{unitnum,1}, TanovaLikRatbsl{unitnum,1}] = find_ancova_rating (mfrlikRat1bsl, mfrlikRat2bsl, val_rat, val_rat_afc);
            [PanovaA1bsl{unitnum,1}, TanovaA1bsl{unitnum,1}] = find_ancova_2AFC  (mfrA11bsl, mfrA12bsl, val_2AFC, val_A, val_B, val_C, val_U);
            [PanovaBbsl{unitnum,1}, TanovaBbsl{unitnum,1}]   = find_ancova_2AFC  (mfrB1bsl,  mfrB2bsl, val_2AFC, val_A, val_B, val_C, val_U);
            [PanovaA2bsl{unitnum,1}, TanovaA2bsl{unitnum,1}]  = find_ancova_2AFC  (mfrA21bsl, mfrA22bsl, val_2AFC, val_A, val_B, val_C, val_U);
            [PanovaLikAFCbsl{unitnum,1}, TanovaLikAFCbsl{unitnum,1}] = find_ancova_2AFC  (mfrlikAFC1bsl, mfrlikAFC2bsl, val_2AFC, val_A, val_B, val_C_Rank, val_U);
            
            [PanovaA1_Bbsl{unitnum,1}, TanovaA1_Bbsl{unitnum,1}] = find_ancova_2AFC  (mfrA11bsl-mfrB1bsl, mfrA12bsl-mfrB2bsl, val_2AFC, val_A, val_B, val_C, val_U);
            [PanovaB_A2bsl{unitnum,1}, TanovaB_A2bsl{unitnum,1}] = find_ancova_2AFC  (mfrB1bsl-mfrA21bsl, mfrB2bsl-mfrA22bsl, val_2AFC, val_A, val_B, val_C, val_U);
            [PanovaA1_Babsbsl{unitnum,1}, TanovaA1_Babsbsl{unitnum,1}]  = find_ancova_2AFC  (abs(mfrA11bsl-mfrB1bsl), abs(mfrA12bsl-mfrB2bsl), val_2AFC, val_A, val_B, val_C, val_U);
            [PanovaB_A2absbsl{unitnum,1}, TanovaB_A2absbsl{unitnum,1} ] = find_ancova_2AFC  (abs(mfrB1bsl-mfrA21bsl), abs(mfrB2bsl-mfrA22bsl), val_2AFC, val_A, val_B, val_C, val_U);
            

            unitnum = unitnum+1;
        end
        
        
        
        
    end

cd (folder_to_save)

name_to_saveCorr=[cell2mat(value_signal), '_' typeCorr '_' num2str(rsp_t_AFC(1)) '_' num2str(rsp_t_AFC(2)) '_' num2str(rsp_t_likertAFC(1)) '_' num2str(rsp_t_likertAFC(2)) '_' num2str(rsp_t_rating(1)) '_' num2str(rsp_t_rating(2)), '_' num2str(rsp_t_likertRat(1)) '_' num2str(rsp_t_likertRat(2))];
name_to_saveCorrbsl=[name_to_saveCorr '_' num2str(rsp_t_bsl(1)) '_' num2str(rsp_t_bsl(2)), 'bsl'];
name_to_saveAncova=[cell2mat(value_signal), '_ancova_' num2str(rsp_t_AFC(1)) '_' num2str(rsp_t_AFC(2)) '_' num2str(rsp_t_likertAFC(1)) '_' num2str(rsp_t_likertAFC(2)) '_' num2str(rsp_t_rating(1)) '_' num2str(rsp_t_rating(2)), '_' num2str(rsp_t_likertRat(1)) '_' num2str(rsp_t_likertRat(2))];
name_to_saveAncovabsl=[name_to_saveAncova '_' num2str(rsp_t_bsl(1)) '_' num2str(rsp_t_bsl(2)), 'bsl'];


save (name_to_saveCorr, 'CorrRat', 'CorrLikRat', 'CorrA1', 'CorrB', 'CorrA2', 'CorrLikAFC', 'CorrA1_B', 'CorrB_A2', 'CorrA1_Babs', 'CorrB_A2abs',  'pvalRat', 'pvalLikRat', 'pvalA1', 'pvalB', 'pvalA2', 'pvalLikAFC', 'pvalA1_B', 'pvalB_A2', 'pvalA1_Babs', 'pvalB_A2abs');
save (name_to_saveCorrbsl, 'CorrRatbsl', 'CorrLikRatbsl', 'CorrA1bsl', 'CorrBbsl', 'CorrA2bsl', 'CorrLikAFCbsl', 'CorrA1_Bbsl', 'CorrB_A2bsl', 'CorrA1_Babsbsl', 'CorrB_A2absbsl', 'pvalRatbsl', 'pvalLikRatbsl', 'pvalA1bsl', 'pvalBbsl', 'pvalA2bsl', 'pvalLikAFCbsl', 'pvalA1_Bbsl', 'pvalB_A2bsl', 'pvalA1_Babsbsl', 'pvalB_A2absbsl');
save (name_to_saveAncova, 'PanovaRat', 'PanovaLikRat', 'PanovaA1', 'PanovaB', 'PanovaA2', 'PanovaLikAFC', 'PanovaA1_B', 'PanovaB_A2', 'PanovaA1_Babs', 'PanovaB_A2abs',  'TanovaRat', 'TanovaLikRat', 'TanovaA1', 'TanovaB', 'TanovaA2', 'TanovaLikAFC', 'TanovaA1_B', 'TanovaB_A2', 'TanovaA1_Babs', 'TanovaB_A2abs');
save (name_to_saveAncovabsl, 'PanovaRatbsl', 'PanovaLikRatbsl', 'PanovaA1bsl', 'PanovaBbsl', 'PanovaA2bsl', 'PanovaLikAFCbsl', 'PanovaA1_Bbsl', 'PanovaB_A2bsl', 'PanovaA1_Babsbsl', 'PanovaB_A2absbsl', 'TanovaRatbsl', 'TanovaLikRatbsl', 'TanovaA1bsl', 'TanovaBbsl', 'TanovaA2bsl', 'TanovaLikAFCbsl', 'TanovaA1_Bbsl', 'TanovaB_A2bsl', 'TanovaA1_Babsbsl', 'TanovaB_A2absbsl');

cd ..
toc      
    
else
    error('Not the same length of behavioural and neural data')
end


function [CorrRat, pvalRat]= find_corr_rating (mfrRat1, mfrRat2, val_rat_out, val_rat_in, typeCorr)

% val out is always the answer from the patient 
% val in is the input, so the expected value from the stimulus, and can be
% coded as as 4 types, ranting, ranking, AFC and total

CorrRat=nan(1, 8*2);
pvalRat=nan(1, 8*2);

for k=1:2
    
    [CorrRat(1, 1*2 - mod(k, 2)), pvalRat(1, 1*2 - mod(k, 2))] =  corr(mfrRat1, val_rat_out(:, k), 'type', typeCorr);
    [CorrRat(1, 2*2 - mod(k, 2)), pvalRat(1, 2*2 - mod(k, 2))] =  corr(mfrRat2, val_rat_out(:, k), 'type', typeCorr);
    [CorrRat(1, 3*2 - mod(k, 2)), pvalRat(1, 3*2 - mod(k, 2))] =  corr(mfrRat1, val_rat_in(:, k), 'type', typeCorr);
    [CorrRat(1, 4*2 - mod(k, 2)), pvalRat(1, 4*2 - mod(k, 2))] =  corr(mfrRat2, val_rat_in(:, k), 'type', typeCorr);
    
    [CorrRat(1, 5*2 - mod(k, 2)), pvalRat(1, 5*2 - mod(k, 2))] =  corr(mfrRat1, abs(val_rat_out(:, k)), 'type', typeCorr);
    [CorrRat(1, 6*2 - mod(k, 2)), pvalRat(1, 6*2 - mod(k, 2))] =  corr(mfrRat2, abs(val_rat_out(:, k)), 'type', typeCorr);
    [CorrRat(1, 7*2 - mod(k, 2)), pvalRat(1, 7*2 - mod(k, 2))] =  corr(mfrRat1, abs(val_rat_in(:, k)), 'type', typeCorr);
    [CorrRat(1, 8*2 - mod(k, 2)), pvalRat(1, 8*2 - mod(k, 2))] =  corr(mfrRat2, abs(val_rat_in(:, k)), 'type', typeCorr);
    
end



function  [Corr2AFC, pval2AFC]= find_corr_2AFC (mfr1, mfr2, val_2AFCout, val_A, val_B, val_C, val_U, typeCorr)

% val_out is always the answer from the patient 
% val A and B is the input, so the expected value from the stimulus, and can be
% coded as as 4 types, ranting, ranking, AFC and total
% val C is the iexpected value (input) of the chosen stimulus, and U of the
% unchosen one


Corr2AFC = nan(1, 30*2);
pval2AFC = nan(1, 30*2);

for k=1:2
    
    
    % corr output value, so the revealed preferences
    [Corr2AFC(1, 1*2 - mod(k, 2)), pval2AFC(1, 1*2 - mod(k, 2))] = corr(mfr1, val_2AFCout(:, k), 'type', typeCorr);
    [Corr2AFC(1, 2*2 - mod(k, 2)), pval2AFC(1, 2*2 - mod(k, 2))] = corr(mfr2, val_2AFCout(:, k), 'type', typeCorr);
    
    % corr input value during A
    [Corr2AFC(1, 3*2 - mod(k, 2)), pval2AFC(1, 3*2 - mod(k, 2))] = corr(mfr1, val_A(:, k), 'type', typeCorr);
    [Corr2AFC(1, 4*2 - mod(k, 2)), pval2AFC(1, 4*2 - mod(k, 2))] = corr(mfr2, val_A(:, k), 'type', typeCorr);
    
    % corr input value during B
    [Corr2AFC(1, 5*2 - mod(k, 2)), pval2AFC(1, 5*2 - mod(k, 2))] = corr(mfr1, val_B(:, k), 'type', typeCorr);
    [Corr2AFC(1, 6*2 - mod(k, 2)), pval2AFC(1, 6*2 - mod(k, 2))] = corr(mfr2, val_B(:, k), 'type', typeCorr);
    
    % corr absolute preference
    [Corr2AFC(1, 7*2 - mod(k, 2)), pval2AFC(1, 7*2 - mod(k, 2))] = corr(mfr1, abs(val_2AFCout(:, k)), 'type', typeCorr);
    [Corr2AFC(1, 8*2 - mod(k, 2)), pval2AFC(1, 8*2 - mod(k, 2))] = corr(mfr2, abs(val_2AFCout(:, k)), 'type', typeCorr);
    
    % corr salience of A
    [Corr2AFC(1, 9*2 - mod(k, 2)), pval2AFC(1, 9*2 - mod(k, 2))] = corr(mfr1, abs(val_A(:, k)), 'type', typeCorr);
    [Corr2AFC(1, 10*2 - mod(k, 2)), pval2AFC(1, 10*2 - mod(k, 2))] = corr(mfr2, abs(val_A(:, k)), 'type', typeCorr);
    
    % corr salience of B
    [Corr2AFC(1, 11*2 - mod(k, 2)), pval2AFC(1, 11*2 - mod(k, 2))] = corr(mfr1, abs(val_B(:, k)), 'type', typeCorr);
    [Corr2AFC(1, 12*2 - mod(k, 2)), pval2AFC(1, 12*2 - mod(k, 2))] = corr(mfr2, abs(val_B(:, k)), 'type', typeCorr);
    
    % corr difference value A-B
    [Corr2AFC(1, 13*2 - mod(k, 2)), pval2AFC(1, 13*2 - mod(k, 2))] = corr(mfr1, val_A(:, k) - val_B(:, k), 'type', typeCorr);
    [Corr2AFC(1, 14*2 - mod(k, 2)), pval2AFC(1, 14*2 - mod(k, 2))] = corr(mfr2, val_A(:, k) - val_B(:, k), 'type', typeCorr);

    % corr absolute difference value A-B
    [Corr2AFC(1, 15*2 - mod(k, 2)), pval2AFC(1, 15*2 - mod(k, 2))] = corr(mfr1, abs(val_B(:, k) - val_A(:, k)), 'type', typeCorr);
    [Corr2AFC(1, 16*2 - mod(k, 2)), pval2AFC(1, 16*2 - mod(k, 2))] = corr(mfr2, abs(val_B(:, k) - val_A(:, k)), 'type', typeCorr);
    
    % corr sum A-B
    [Corr2AFC(1, 17*2 - mod(k, 2)), pval2AFC(1, 17*2 - mod(k, 2))] = corr(mfr1, (val_B(:, k) + val_A(:, k)), 'type', typeCorr);
    [Corr2AFC(1, 18*2 - mod(k, 2)), pval2AFC(1, 18*2 - mod(k, 2))] = corr(mfr2, (val_B(:, k) + val_A(:, k)), 'type', typeCorr);
    
    % corr sum of salience A-B
    [Corr2AFC(1, 19*2 - mod(k, 2)), pval2AFC(1, 19*2 - mod(k, 2))] = corr(mfr1, abs(val_B(:, k)) + abs(val_A(:, k)), 'type', typeCorr);
    [Corr2AFC(1, 20*2 - mod(k, 2)), pval2AFC(1, 20*2 - mod(k, 2))] = corr(mfr2, abs(val_B(:, k)) + abs(val_A(:, k)), 'type', typeCorr);
    
    % corr absolute of sum A-B
    [Corr2AFC(1, 21*2 - mod(k, 2)), pval2AFC(1, 21*2 - mod(k, 2))] = corr(mfr1, abs(val_B(:, k) + val_A(:, k)), 'type', typeCorr);
    [Corr2AFC(1, 22*2 - mod(k, 2)), pval2AFC(1, 22*2 - mod(k, 2))] = corr(mfr2, abs(val_B(:, k) + val_A(:, k)), 'type', typeCorr);
    
    % corr chosen value 
    [Corr2AFC(1, 23*2 - mod(k, 2)), pval2AFC(1, 23*2 - mod(k, 2))] = corr(mfr1, val_C(:, k), 'type', typeCorr);
    [Corr2AFC(1, 24*2 - mod(k, 2)), pval2AFC(1, 24*2 - mod(k, 2))] = corr(mfr2, val_C(:, k), 'type', typeCorr);
    
    % corr unchosen value
    [Corr2AFC(1, 25*2 - mod(k, 2)), pval2AFC(1, 25*2 - mod(k, 2))] = corr(mfr1, val_U(:, k), 'type', typeCorr);
    [Corr2AFC(1, 26*2 - mod(k, 2)), pval2AFC(1, 26*2 - mod(k, 2))] = corr(mfr2, val_U(:, k), 'type', typeCorr);
    
    % corr chosen salience
    [Corr2AFC(1, 27*2 - mod(k, 2)), pval2AFC(1, 27*2 - mod(k, 2))] = corr(mfr1, abs(val_C(:, k)), 'type', typeCorr);
    [Corr2AFC(1, 28*2 - mod(k, 2)), pval2AFC(1, 28*2 - mod(k, 2))] = corr(mfr2, abs(val_C(:, k)), 'type', typeCorr);
    
    % corr unchosen salience    
    [Corr2AFC(1, 29*2 - mod(k, 2)), pval2AFC(1, 29*2 - mod(k, 2))] = corr(mfr1, abs(val_U(:, k)), 'type', typeCorr);
    [Corr2AFC(1, 30*2 - mod(k, 2)), pval2AFC(1, 30*2 - mod(k, 2))] = corr(mfr2, abs(val_U(:, k)), 'type', typeCorr);
    
end


function [PanovaRat, TanovaRat]= find_ancova_rating (mfrRat1, mfrRat2, val_rat_out, val_rat_in)

% val out is always the answer from the patient 
% val in is the input, so the expected value from the stimulus, and can be
% coded as as 4 types, ranting, ranking, AFC and total


PanovaRat=cell(1, 8*2);
TanovaRat=cell(1, 8*2);

for k=1:2
    
    [PanovaRat{1, 1*2 - mod(k, 2)}, TanovaRat{1, 1*2 - mod(k, 2)}] =  anovan(mfrRat1, val_rat_out(:, k), 'continuous',[1],'varnames',{'val_rat_out'},'display','off');
    [PanovaRat{1, 2*2 - mod(k, 2)}, TanovaRat{1, 2*2 - mod(k, 2)}] =  anovan(mfrRat2, val_rat_out(:, k), 'continuous',[1],'varnames',{'val_rat_out'},'display','off');
    [PanovaRat{1, 3*2 - mod(k, 2)}, TanovaRat{1, 3*2 - mod(k, 2)}] =  anovan(mfrRat1, val_rat_in(:, k), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 4*2 - mod(k, 2)}, TanovaRat{1, 4*2 - mod(k, 2)}] =  anovan(mfrRat2, val_rat_in(:, k), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    
    [PanovaRat{1, 5*2 - mod(k, 2)}, TanovaRat{1, 5*2 - mod(k, 2)}] =  anovan(mfrRat1, abs(val_rat_out(:, k)), 'continuous',[1],'varnames',{'sal_rat_out'},'display','off');
    [PanovaRat{1, 6*2 - mod(k, 2)}, TanovaRat{1, 6*2 - mod(k, 2)}] =  anovan(mfrRat2, abs(val_rat_out(:, k)), 'continuous',[1],'varnames',{'sal_rat_out'},'display','off');
    [PanovaRat{1, 7*2 - mod(k, 2)}, TanovaRat{1, 7*2 - mod(k, 2)}] =  anovan(mfrRat1, abs(val_rat_in(:, k)), 'continuous',[1],'varnames',{'sal_rat_in'},'display','off');
    [PanovaRat{1, 8*2 - mod(k, 2)}, TanovaRat{1, 8*2 - mod(k, 2)}] =  anovan(mfrRat2, abs(val_rat_in(:, k)), 'continuous',[1],'varnames',{'sal_rat_in'},'display','off');
    
end



function  [Panova2AFC, TanovaAFC]= find_ancova_2AFC (mfr1, mfr2, val_2AFCout, val_A, val_B, val_C, val_U)

% val_out is always the answer from the patient 
% val A and B is the input, so the expected value from the stimulus, and can be
% coded as as 4 types, ranting, ranking, AFC and total
% val C is the iexpected value (input) of the chosen stimulus, and U of the
% unchosen one

Panova2AFC = cell(1, 36*2);
TanovaAFC = cell(1, 36*2);

for k=1:2
    
    
    % corr output value, so the revealed preferences
    [Panova2AFC{1, 1*2 - mod(k, 2)}, TanovaAFC{1, 1*2 - mod(k, 2)}] = anovan(mfr1, val_2AFCout(:, k), 'continuous',[1],'varnames',{'val_2AFCout'},'display','off');
    [Panova2AFC{1, 2*2 - mod(k, 2)}, TanovaAFC{1, 2*2 - mod(k, 2)}] = anovan(mfr2, val_2AFCout(:, k),'continuous',[1],'varnames',{'val_2AFCout'},'display','off');
    
    % corr input value during A
    [Panova2AFC{1, 3*2 - mod(k, 2)}, TanovaAFC{1, 3*2 - mod(k, 2)}] = anovan(mfr1, val_A(:, k), 'continuous',[1],'varnames',{'val_A'},'display','off');
    [Panova2AFC{1, 4*2 - mod(k, 2)}, TanovaAFC{1, 4*2 - mod(k, 2)}] = anovan(mfr2, val_A(:, k), 'continuous',[1],'varnames',{'val_A'},'display','off');
    
    % corr input value during B
    [Panova2AFC{1, 5*2 - mod(k, 2)}, TanovaAFC{1, 5*2 - mod(k, 2)}] = anovan(mfr1, val_B(:, k), 'continuous',[1],'varnames',{'val_B'},'display','off');
    [Panova2AFC{1, 6*2 - mod(k, 2)}, TanovaAFC{1, 6*2 - mod(k, 2)}] = anovan(mfr2, val_B(:, k), 'continuous',[1],'varnames',{'val_B'},'display','off');
    
    % corr absolute preference
    [Panova2AFC{1, 7*2 - mod(k, 2)}, TanovaAFC{1, 7*2 - mod(k, 2)}] = anovan(mfr1, abs(val_2AFCout(:, k)), 'continuous',[1],'varnames',{'val_2AFCout'},'display','off');
    [Panova2AFC{1, 8*2 - mod(k, 2)}, TanovaAFC{1, 8*2 - mod(k, 2)}] = anovan(mfr2, abs(val_2AFCout(:, k)), 'continuous',[1],'varnames',{'val_2AFCout'},'display','off');
    
    % corr salience of A
    [Panova2AFC{1, 9*2 - mod(k, 2)}, TanovaAFC{1, 9*2 - mod(k, 2)}] = anovan(mfr1, abs(val_A(:, k)), 'continuous',[1],'varnames',{'sal_A'},'display','off');
    [Panova2AFC{1, 10*2 - mod(k, 2)}, TanovaAFC{1, 10*2 - mod(k, 2)}] = anovan(mfr2, abs(val_A(:, k)), 'continuous',[1],'varnames',{'sal_A'},'display','off');
    
    % corr salience of B
    [Panova2AFC{1, 11*2 - mod(k, 2)}, TanovaAFC{1, 11*2 - mod(k, 2)}] = anovan(mfr1, abs(val_B(:, k)), 'continuous',[1],'varnames',{'sal_B'},'display','off');
    [Panova2AFC{1, 12*2 - mod(k, 2)}, TanovaAFC{1, 12*2 - mod(k, 2)}] = anovan(mfr2, abs(val_B(:, k)), 'continuous',[1],'varnames',{'sal_B'},'display','off');
    
    % corr difference value A-B
    [Panova2AFC{1, 13*2 - mod(k, 2)}, TanovaAFC{1, 13*2 - mod(k, 2)}] = anovan(mfr1, val_A(:, k) - val_B(:, k), 'continuous',[1],'varnames',{'diff'},'display','off');
    [Panova2AFC{1, 14*2 - mod(k, 2)}, TanovaAFC{1, 14*2 - mod(k, 2)}] = anovan(mfr2, val_A(:, k) - val_B(:, k), 'continuous',[1],'varnames',{'diff'},'display','off');

    % corr absolute difference value A-B
    [Panova2AFC{1, 15*2 - mod(k, 2)}, TanovaAFC{1, 15*2 - mod(k, 2)}] = anovan(mfr1, abs(val_B(:, k) - val_A(:, k)), 'continuous',[1],'varnames',{'abs_diff'},'display','off');
    [Panova2AFC{1, 16*2 - mod(k, 2)}, TanovaAFC{1, 16*2 - mod(k, 2)}] = anovan(mfr2, abs(val_B(:, k) - val_A(:, k)), 'continuous',[1],'varnames',{'abs_diff'},'display','off');
    
    % corr sum A-B
    [Panova2AFC{1, 17*2 - mod(k, 2)}, TanovaAFC{1, 17*2 - mod(k, 2)}] = anovan(mfr1, (val_B(:, k) + val_A(:, k)), 'continuous',[1],'varnames',{'sum'},'display','off');
    [Panova2AFC{1, 18*2 - mod(k, 2)}, TanovaAFC{1, 18*2 - mod(k, 2)}] = anovan(mfr2, (val_B(:, k) + val_A(:, k)), 'continuous',[1],'varnames',{'sum'},'display','off');
    
    % corr sum of salience A-B
    [Panova2AFC{1, 19*2 - mod(k, 2)}, TanovaAFC{1, 19*2 - mod(k, 2)}] = anovan(mfr1, abs(val_B(:, k)) + abs(val_A(:, k)), 'continuous',[1],'varnames',{'abs_sum'},'display','off');
    [Panova2AFC{1, 20*2 - mod(k, 2)}, TanovaAFC{1, 20*2 - mod(k, 2)}] = anovan(mfr2, abs(val_B(:, k)) + abs(val_A(:, k)), 'continuous',[1],'varnames',{'abs_sum'},'display','off');
    
    % corr absolute of sum A-B
    [Panova2AFC{1, 21*2 - mod(k, 2)}, TanovaAFC{1, 21*2 - mod(k, 2)}] = anovan(mfr1, abs(val_B(:, k) + val_A(:, k)), 'continuous',[1],'varnames',{'abs_abs_sum'},'display','off');
    [Panova2AFC{1, 22*2 - mod(k, 2)}, TanovaAFC{1, 22*2 - mod(k, 2)}] = anovan(mfr2, abs(val_B(:, k) + val_A(:, k)), 'continuous',[1],'varnames',{'sal_rat_in'},'display','off');
    
    % corr chosen value 
    [Panova2AFC{1, 23*2 - mod(k, 2)}, TanovaAFC{1, 23*2 - mod(k, 2)}] = anovan(mfr1, val_C(:, k), 'continuous',[1],'varnames',{'val_C'},'display','off');
    [Panova2AFC{1, 24*2 - mod(k, 2)}, TanovaAFC{1, 24*2 - mod(k, 2)}] = anovan(mfr2, val_C(:, k), 'continuous',[1],'varnames',{'val_C'},'display','off');
    
    % corr unchosen value
    [Panova2AFC{1, 25*2 - mod(k, 2)}, TanovaAFC{1, 25*2 - mod(k, 2)}] = anovan(mfr1, val_U(:, k), 'continuous',[1],'varnames',{'val_U'},'display','off');
    [Panova2AFC{1, 26*2 - mod(k, 2)}, TanovaAFC{1, 26*2 - mod(k, 2)}] = anovan(mfr2, val_U(:, k), 'continuous',[1],'varnames',{'val_U'},'display','off');
    
    % corr chosen salience
    [Panova2AFC{1, 27*2 - mod(k, 2)}, TanovaAFC{1, 27*2 - mod(k, 2)}] = anovan(mfr1, abs(val_C(:, k)), 'continuous',[1],'varnames',{'sal_C'},'display','off');
    [Panova2AFC{1, 28*2 - mod(k, 2)}, TanovaAFC{1, 28*2 - mod(k, 2)}] = anovan(mfr2, abs(val_C(:, k)), 'continuous',[1],'varnames',{'sal_C'},'display','off');
    
    % corr unchosen salience    
    [Panova2AFC{1, 29*2 - mod(k, 2)}, TanovaAFC{1, 29*2 - mod(k, 2)}] = anovan(mfr1, abs(val_U(:, k)), 'continuous',[1],'varnames',{'sal_U'},'display','off');
    [Panova2AFC{1, 30*2 - mod(k, 2)}, TanovaAFC{1, 30*2 - mod(k, 2)}] = anovan(mfr2, abs(val_U(:, k)), 'continuous',[1],'varnames',{'sal_U'},'display','off');
    
    
    %% 2 way ancovas
    
    % corr output value, so the revealed preferences
    [Panova2AFC{1, 31*2 - mod(k, 2)}, TanovaAFC{1, 31*2 - mod(k, 2)}] = anovan(mfr1, {val_2AFCout(:, k), val_A(:, k), val_B(:, k)}, 'continuous',[1,2,3],'varnames',{'val_2AFCout', 'val_A', 'val_B'},'display','off');
    [Panova2AFC{1, 32*2 - mod(k, 2)}, TanovaAFC{1, 32*2 - mod(k, 2)}] = anovan(mfr2, {val_2AFCout(:, k), val_A(:, k), val_B(:, k)},'continuous',[1,2,3],'varnames',{'val_2AFCout','val_A', 'val_B'},'display','off');
         
    % corr output value, so the revealed preferences
    [Panova2AFC{1, 33*2 - mod(k, 2)}, TanovaAFC{1, 33*2 - mod(k, 2)}] = anovan(mfr1, {val_A(:, k), val_B(:, k)}, 'continuous',[1,2],'varnames',{'val_A', 'val_B' },'display','off');
    [Panova2AFC{1, 34*2 - mod(k, 2)}, TanovaAFC{1, 34*2 - mod(k, 2)}] = anovan(mfr2, {val_A(:, k), val_B(:, k)},'continuous',[1,2],'varnames',{'val_A', 'val_B'},'display','off');
        
    % corr output value, so the revealed preferences
    [Panova2AFC{1, 35*2 - mod(k, 2)}, TanovaAFC{1, 35*2 - mod(k, 2)}] = anovan(mfr1, {val_B(:, k),  val_A(:, k)}, 'continuous',[1,2],'varnames',{'val_B', 'val_A'},'display','off');
    [Panova2AFC{1, 36*2 - mod(k, 2)}, TanovaAFC{1, 36*2 - mod(k, 2)}] = anovan(mfr2, {val_B(:, k),  val_A(:, k)}, 'continuous',[1,2],'varnames',{'val_B', 'val_A'},'display','off');
        
end

