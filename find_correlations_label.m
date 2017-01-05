function find_correlations_label

tic
dbstop if error

%% settings 


rsp_t_AFC= [0,1000]; %use at the beginning this time period to find the corrlations, this could change
rsp_t_likertAFC= [0,1000]; %use at the beginning this time period to find the corrlations, this could change
rsp_t_rating= [0,1000]; %use at the beginning this time period to find the corrlations, this could change
rsp_t_likertRat= [0,1000]; %use at the beginning this time period to find the corrlations, this could change
rsp_t_bsl= [-400,0];

value_signal='All'; % this is the value signal used in this analysis, thus, Tot, Rat, Rank or AFC.  If 'All', then it loop over the four types of signlas and saves all of them 
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

end



function main(value_signal, rasts, rsp_t_AFC, rsp_t_rating, rsp_t_likertRat, rsp_t_likertAFC, rsp_t_bsl, folder_to_save, typeCorr, namefiles, outfiles, behavDir)

if length(namefiles) == length(outfiles)
    

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
            
            
            
            
          %% convert MFR to mean of MFR, so the mean of every time it was
            % shown. 
            
            MFRB1= nan(20,10); MFRA11= nan(20,10); MFRA21= nan(20,10); MFRB2= nan(20,10); MFRA12= nan(20,10); MFRA22= nan(20,10);
            MFRB1bsl= nan(20,10); MFRB2bsl= nan(20,10); MFRA11bsl= nan(20,10); MFRA12bsl= nan(20,10); MFRA21bsl= nan(20,10); MFRA22bsl= nan(20,10);
            MFRRat1= nan(20,3); MFRlikRat1= nan(20,3);  MFRRat2= nan(20,3); MFRlikRat2= nan(20,3);
            MFRRat1bsl= nan(20,3); MFRlikRat1bsl= nan(20,3);  MFRRat2bsl= nan(20,3); MFRlikRat2bsl= nan(20,3);
            
            for jj=[3,6]
                for i=1:20
                    
                    
                    if jj==3
                        
                        ind3=find(cell2mat(satcue{jj-1}(:,2))==i);
                        MFRRat1(i,:)   = (mfrRat1(ind3,1));
                        MFRlikRat1(i,:)   = (mfrlikRat1(ind3,1));
                        MFRRat1bsl(i,:)    = ((mfrRat1(ind3,1) - bslRat1(ind3,1)));
                        MFRlikRat1bsl(i,:) = ((mfrlikRat1(ind3,1) - bslRat1(ind3,1)));
                        
                    elseif jj==6
                        
                        ind3=find(cell2mat(satcue{jj-1}(:,2))==i);
                        MFRRat2(i,:)   = (mfrRat2(ind3,1));
                        MFRlikRat2(i,:)   = (mfrlikRat2(ind3,1));
                        MFRRat2bsl(i,:)    = ((mfrRat2(ind3,1) - bslRat2(ind3,1)));
                        MFRlikRat2bsl(i,:) = ((mfrlikRat2(ind3,1) - bslRat2(ind3,1)));
                    end
                    
                    
                    
                    indl = find(cell2mat(satcue{jj}(:,2))==i);            % when this stimulus was on left
                    indr = find(cell2mat(satcue{jj}(:,3))==i);           % when this stimulus was on right
                    lr_to_ab = cell2mat(satcue{jj}(indl,8));              % conversion factor from left_rigth to A_B, this factor saved on the column 8
                    lr_to_ab2 = cell2mat(satcue{jj}(indr,8));            % conversion factor from left_rigth to A_B
                    lr_to_cu = cell2mat(satcue{jj}(:,14))>0;              % conversion factor from left_rigth to A_B, just if positive or negative the output, 1 if right chosen, zero if left chosen
                    indA = cat(1,((lr_to_ab) .* indl), (~lr_to_ab2).* indr);  indA(indA==0)= [];  %indices for A
                    indB = cat(1,((~lr_to_ab) .* indl), (lr_to_ab2).* indr);  indB(indB==0)= [];  %indices for B
                    indC = cat(1,indl(~lr_to_cu(indl)), indr(lr_to_cu(indr)));   % indeces for chosen, 1 if right chosen, so 1*indr, and zero if left chosen, so ~0 * indl
                    indU = cat(1,indl(lr_to_cu(indl)), indr(~lr_to_cu(indr)));
                    
                    
                    if jj==3
                        
                        MFRB1(i, 1:length(indB))   = mfrB1(indB,1);
                        MFRA11(i, 1:length(indA))   = mfrA11(indA,1);
                        MFRA21(i, 1:length(indA))   = mfrA21(indA,1);
%                         MFRlikAFC1(i,1) = mfrlikAFC1(ind3,1);
                       
                        MFRB1bsl(i, 1:length(indB))   = (mfrB1(indB,1)- MFRbsl1(indB,1));
                        MFRA11bsl(i, 1:length(indA))   = (mfrA11(indA,1)- MFRbsl1(indA,1));
                        MFRA21bsl(i, 1:length(indA))   = (mfrA21(indA,1)- MFRbsl1(indA,1));
%                         MFRlikAFC1bsl(i,1) = (mfrlikAFC1(ind3,1) - bslRat1(ind3,1));
                    
                    elseif jj==6
                        
                        MFRB2(i, 1:length(indB))   = mfrB2(indB,1);
                        MFRA12(i, 1:length(indA))   = mfrA12(indA,1);
                        MFRA22(i, 1:length(indA))   = mfrA22(indA,1);
%                         MFRlikAFC2(i,1) = mfrlikAFC2(ind3,1);
                                                                  
                        
                        MFRB2bsl(i, 1:length(indB))   = (mfrB2(indB,1)- MFRbsl2(indB,1));
                        MFRA12bsl(i, 1:length(indA))   = (mfrA12(indA,1)- MFRbsl2(indA,1));
                        MFRA22bsl(i, 1:length(indA))   = (mfrA22(indA,1)- MFRbsl2(indA,1));
%                         MFRlikAFC2bsl(i,1) = (mfrlikAFC2(ind3,1) - bslRat2(ind3,1));
                        
                    end
                    
                end
                
           % compute value of each stimulus in order 
           
           lab_valRat(:,jj/3)=  ranking{jj-1}(:,4); % ranking
           lab_valAFC(:,jj/3)= ranking{jj}(:,4);
           lab_valRank(:,jj/3)= ranking{jj}(:,3);
           lab_valTot(:,jj/3)= ranking{jj}(:,5);
           end
           

                
      %%
            
            
            
            % select the values to use

            eval(['label_val = lab_val', cell2mat(value_signal)]);
            

            % make correlation matrices for each trial (190 for AFC, and 60 for rating)

            % 20 labels analysus of rating, pooling together 3 times shown each
            [CorrLabelVal(unitnum, :), pvalLabelVal(unitnum, :)] = find_corr_rating_label (MFRRat1, MFRRat2, MFRB1, MFRB2, MFRA11, MFRA12, MFRA21, MFRA22, MFRlikRat1, MFRlikRat2,  label_val, typeCorr, 0);
            [CorrLabelSal(unitnum, :), pvalLabelSal(unitnum, :)] = find_corr_rating_label (MFRRat1, MFRRat2, MFRB1, MFRB2, MFRA11, MFRA12, MFRA21, MFRA22, MFRlikRat1, MFRlikRat2,  label_val, typeCorr, 1);
            [CorrLabelBslVal(unitnum, :), pvalLabelBslVal(unitnum, :)] = find_corr_rating_label (MFRRat1bsl, MFRRat2bsl,MFRB1bsl, MFRB2bsl, MFRA11bsl, MFRA12bsl, MFRA21bsl, MFRA22bsl, MFRlikRat1bsl, MFRlikRat2bsl,  label_val, typeCorr,0);
            [CorrLabelBslSal(unitnum, :), pvalLabelBslSal(unitnum, :)] = find_corr_rating_label (MFRRat1bsl, MFRRat2bsl,MFRB1bsl, MFRB2bsl, MFRA11bsl, MFRA12bsl, MFRA21bsl, MFRA22bsl, MFRlikRat1bsl, MFRlikRat2bsl,  label_val, typeCorr, 1);

            
            
            
            %% ancova analysis
            
            

            % 20 labels analysus of rating, pooling together 3 times shown each
            [PanovaLabelVal{unitnum,1}, TanovaLabelVal{unitnum,1}] = find_ancova_rating_label (MFRRat1, MFRRat2, MFRB1, MFRB2, MFRA11, MFRA12, MFRA21, MFRA22, MFRlikRat1, MFRlikRat2, label_val, 0);
            [PanovaLabelSal{unitnum,1}, TanovaLabelSal{unitnum,1}] = find_ancova_rating_label (MFRRat1, MFRRat2, MFRB1, MFRB2, MFRA11, MFRA12, MFRA21, MFRA22, MFRlikRat1, MFRlikRat2, label_val, 1);
            [PanovaLabelBslVal{unitnum,1}, TanovaLabelBslVal{unitnum,1}] = find_ancova_rating_label (MFRRat1bsl, MFRRat2bsl,MFRB1bsl, MFRB2bsl, MFRA11bsl, MFRA12bsl, MFRA21bsl, MFRA22bsl, MFRlikRat1bsl, MFRlikRat2bsl,  label_val, 0);
            [PanovaLabelBslSal{unitnum,1}, TanovaLabelBslSal{unitnum,1}] = find_ancova_rating_label (MFRRat1bsl, MFRRat2bsl,MFRB1bsl, MFRB2bsl, MFRA11bsl, MFRA12bsl, MFRA21bsl, MFRA22bsl, MFRlikRat1bsl, MFRlikRat2bsl,  label_val, 1);

            
            
            unitnum = unitnum+1;
        end
        
        
        
        
    end

cd (folder_to_save)

name_to_saveCorr=[cell2mat(value_signal), '_label_' typeCorr '_' num2str(rsp_t_AFC(1)) '_' num2str(rsp_t_AFC(2)) '_' num2str(rsp_t_likertAFC(1)) '_' num2str(rsp_t_likertAFC(2)) '_' num2str(rsp_t_rating(1)) '_' num2str(rsp_t_rating(2)), '_' num2str(rsp_t_likertRat(1)) '_' num2str(rsp_t_likertRat(2))];
name_to_saveCorrbsl=[name_to_saveCorr '_' num2str(rsp_t_bsl(1)) '_' num2str(rsp_t_bsl(2)), 'bsl'];
name_to_saveAncova=[cell2mat(value_signal), '_label_ancova_' num2str(rsp_t_AFC(1)) '_' num2str(rsp_t_AFC(2)) '_' num2str(rsp_t_likertAFC(1)) '_' num2str(rsp_t_likertAFC(2)) '_' num2str(rsp_t_rating(1)) '_' num2str(rsp_t_rating(2)), '_' num2str(rsp_t_likertRat(1)) '_' num2str(rsp_t_likertRat(2))];
name_to_saveAncovabsl=[name_to_saveAncova '_' num2str(rsp_t_bsl(1)) '_' num2str(rsp_t_bsl(2)), 'bsl'];


save (name_to_saveCorr,  'CorrLabelVal', 'CorrLabelSal',  'pvalLabelVal', 'pvalLabelSal');
save (name_to_saveCorrbsl,'CorrLabelBslVal', 'CorrLabelBslSal', 'pvalLabelBslVal', 'pvalLabelBslSal' );
save (name_to_saveAncova,  'PanovaLabelVal', 'PanovaLabelSal', 'TanovaLabelVal', 'TanovaLabelSal');
save (name_to_saveAncovabsl, 'PanovaLabelBslVal' ,'PanovaLabelBslSal','TanovaLabelBslVal','TanovaLabelBslSal');

cd ..
toc      
    
else
    error('Not the same length of behavioural and neural data')
end
end 



function [CorrRat, pvalRat]= find_corr_rating_label (MFRrat1, MFRrat2, MFRB1, MFRB2, MFRA11, MFRA12, MFRA21, MFRA22,  MFRlikRat1, MFRlikRat2,  val_rat_mean, typeCorr, absolute)

% val out is always the answer from the patient 
% val in is the input, so the expected value from the stimulus, and can be
% coded as as 4 types, ranting, ranking, AFC and total

CorrRat=nan(1, 24);
pvalRat=nan(1, 24);

    
if absolute==1
    val_rat_mean= abs(val_rat_mean);
end

    % R (Rating pic)
    [CorrRat(1, 1), pvalRat(1, 1)] =  corr(mean(MFRrat1,2), val_rat_mean(:, 1), 'type', typeCorr);
    [CorrRat(1, 2), pvalRat(1, 2)] =  corr(mean(MFRrat2,2), val_rat_mean(:, 2), 'type', typeCorr);
    % B
    [CorrRat(1, 3), pvalRat(1, 3)] =  corr(nanmean(MFRB1,2), val_rat_mean(:, 1), 'type', typeCorr);
    [CorrRat(1, 4), pvalRat(1, 4)] =  corr(nanmean(MFRB2,2), val_rat_mean(:, 2), 'type', typeCorr);
    % A1
    [CorrRat(1, 5), pvalRat(1, 5)] =  corr(nanmean(MFRA11,2), val_rat_mean(:, 1), 'type', typeCorr);
    [CorrRat(1, 6), pvalRat(1, 6)] =  corr(nanmean(MFRA12,2), val_rat_mean(:, 2), 'type', typeCorr);
    % A2
    [CorrRat(1, 7), pvalRat(1, 7)] =  corr(nanmean(MFRA21,2), val_rat_mean(:, 1), 'type', typeCorr);
    [CorrRat(1, 8), pvalRat(1, 8)] =  corr(nanmean(MFRA22,2), val_rat_mean(:, 2), 'type', typeCorr);
    % LikRat
    [CorrRat(1, 9), pvalRat(1, 9)] =  corr(mean(MFRlikRat1,2), val_rat_mean(:, 1), 'type', typeCorr);
    [CorrRat(1, 10), pvalRat(1,10)] =  corr(mean(MFRlikRat2,2), val_rat_mean(:, 2), 'type', typeCorr);
    
    
    % pooling combinations
    
    % RA
    [CorrRat(1, 11), pvalRat(1, 11)]  =  corr(nanmean([MFRrat1, MFRA11],2), val_rat_mean(:, 1), 'type', typeCorr);
    [CorrRat(1, 12), pvalRat(1, 12)] =  corr(nanmean([MFRrat2, MFRA12],2), val_rat_mean(:, 2), 'type', typeCorr);
    % RAB
    [CorrRat(1, 13), pvalRat(1, 13)] =  corr(nanmean([MFRrat1, MFRA11, MFRB1],2), val_rat_mean(:, 1), 'type', typeCorr);
    [CorrRat(1, 14), pvalRat(1, 14)] =  corr(nanmean([MFRrat2, MFRA12, MFRB2],2), val_rat_mean(:, 2), 'type', typeCorr);
    % RABA
    [CorrRat(1, 15), pvalRat(1, 15)] =  corr(nanmean([MFRrat1, MFRA11, MFRB1, MFRA21],2), val_rat_mean(:, 1), 'type', typeCorr);
    [CorrRat(1, 16), pvalRat(1, 16)] =  corr(nanmean([MFRrat2, MFRA12, MFRB2, MFRA22],2), val_rat_mean(:, 2), 'type', typeCorr);
    % R-Lik
    [CorrRat(1, 17), pvalRat(1, 17)] =  corr(nanmean([MFRrat1, MFRlikRat1],2), val_rat_mean(:, 1), 'type', typeCorr);
    [CorrRat(1, 18), pvalRat(1, 18)] =  corr(nanmean([MFRrat2, MFRlikRat2],2), val_rat_mean(:, 2), 'type', typeCorr);
    % AB
    [CorrRat(1, 19), pvalRat(1, 19)] =  corr(nanmean([MFRA11, MFRB1],2), val_rat_mean(:, 1), 'type', typeCorr);
    [CorrRat(1, 20), pvalRat(1, 20)] =  corr(nanmean([MFRA12, MFRB2],2), val_rat_mean(:, 2), 'type', typeCorr);
    % ABA
    [CorrRat(1, 21), pvalRat(1, 21)] =  corr(nanmean([MFRA11, MFRB1, MFRA21],2), val_rat_mean(:, 1), 'type', typeCorr);
    [CorrRat(1, 22), pvalRat(1, 22)] =  corr(nanmean([MFRA12, MFRB2, MFRA22],2), val_rat_mean(:, 2), 'type', typeCorr);
    % AA
    [CorrRat(1, 23), pvalRat(1, 23)] =  corr(nanmean([MFRA11,  MFRA21],2), val_rat_mean(:, 1), 'type', typeCorr);
    [CorrRat(1, 24), pvalRat(1, 24)] =  corr(nanmean([MFRA12,  MFRA22],2), val_rat_mean(:, 2), 'type', typeCorr);

    % now both both before and after lunch, I pool them the same way HIll
    % did it, so, doing the mean of both ratings before and after and
    % correlating with the pooled data
    
    % R, A1, B, A2, LikRat
    [CorrRat(1, 25), pvalRat(1, 25)]  =  corr(nanmean([MFRrat1,MFRrat2],2), mean(val_rat_mean,2), 'type', typeCorr);
    [CorrRat(1, 26), pvalRat(1, 26)] =  corr(nanmean([MFRA11, MFRA12],2), mean(val_rat_mean,2), 'type', typeCorr);
    [CorrRat(1, 27), pvalRat(1, 27)] =  corr(nanmean([MFRB1, MFRB2],2), mean(val_rat_mean,2), 'type', typeCorr);
    [CorrRat(1, 28), pvalRat(1, 28)] =  corr(nanmean([MFRA21, MFRA22],2), mean(val_rat_mean,2), 'type', typeCorr);
    [CorrRat(1, 29), pvalRat(1, 29)] =  corr(nanmean([MFRlikRat1, MFRlikRat2],2), mean(val_rat_mean,2), 'type', typeCorr);
    
    % RA, RAB, RABA, R-Lik
    [CorrRat(1, 30), pvalRat(1, 30)] =  corr(nanmean([MFRrat1, MFRrat2, MFRA11, MFRA12],2), mean(val_rat_mean,2), 'type', typeCorr);
    [CorrRat(1, 31), pvalRat(1, 31)] =  corr(nanmean([MFRrat1, MFRrat2, MFRA11, MFRA12, MFRB1, MFRB2],2), mean(val_rat_mean,2), 'type', typeCorr);  
    [CorrRat(1, 32), pvalRat(1, 32)] =  corr(nanmean([MFRrat1, MFRrat2, MFRA11, MFRA12, MFRB1, MFRB2, MFRA21, MFRA22],2), mean(val_rat_mean,2), 'type', typeCorr);  
    [CorrRat(1, 33), pvalRat(1, 33)] =  corr(nanmean([MFRrat1, MFRrat2, MFRlikRat1, MFRlikRat2],2), mean(val_rat_mean,2), 'type', typeCorr);
    
    % AB, ABA, AA
    [CorrRat(1, 34), pvalRat(1, 34)] =  corr(nanmean([MFRA11, MFRA12, MFRB1, MFRB2],2), mean(val_rat_mean,2), 'type', typeCorr);
    [CorrRat(1, 35), pvalRat(1, 35)] =  corr(nanmean([MFRA11, MFRA12, MFRB1, MFRB2, MFRA21, MFRA22],2), mean(val_rat_mean,2), 'type', typeCorr);
    [CorrRat(1, 36), pvalRat(1, 36)] =  corr(nanmean([MFRA11, MFRA12, MFRA21, MFRA22],2), mean(val_rat_mean,2), 'type', typeCorr);

    
end 



function [PanovaRat, TanovaRat]= find_ancova_rating_label (MFRrat1, MFRrat2, MFRB1, MFRB2, MFRA11, MFRA12, MFRA21, MFRA22,  MFRlikRat1, MFRlikRat2,  val_rat_mean, absolute)

% val out is always the answer from the patient
% val in is the input, so the expected value from the stimulus, and can be
% coded as as 4 types, ranting, ranking, AFC and total


PanovaRat=cell(1, 24);
TanovaRat=cell(1, 24);

if absolute==1
    val_rat_mean= abs(val_rat_mean);
end

    
    % R
    [PanovaRat{1, 1}, TanovaRat{1, 1}] =  anovan(mean(MFRrat1,2), val_rat_mean(:, 1), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 2}, TanovaRat{1, 2}] =  anovan(mean(MFRrat2,2), val_rat_mean(:, 2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    % B
    [PanovaRat{1, 3}, TanovaRat{1, 3}] =  anovan(nanmean(MFRB1,2), val_rat_mean(:, 1), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 4}, TanovaRat{1, 4}] =  anovan(nanmean(MFRB2,2), val_rat_mean(:, 2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    % A1
    [PanovaRat{1, 5}, TanovaRat{1, 5}] =  anovan(nanmean(MFRA11,2), val_rat_mean(:, 1), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 6}, TanovaRat{1, 6}] =  anovan(nanmean(MFRA12,2), val_rat_mean(:, 2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    % A2
    [PanovaRat{1, 7}, TanovaRat{1, 7}] =  anovan(nanmean(MFRA21,2), val_rat_mean(:, 1), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 8}, TanovaRat{1, 8}] =  anovan(nanmean(MFRA22,2), val_rat_mean(:, 2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    % LikRat
    [PanovaRat{1, 9}, TanovaRat{1, 9}] =  anovan(mean(MFRlikRat1,2), val_rat_mean(:, 1), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 10}, TanovaRat{1,10}] =  anovan(mean(MFRlikRat2,2), val_rat_mean(:, 2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    
    % poolig combinations
    
    % RA
    [PanovaRat{1, 11}, TanovaRat{1, 11}] =  anovan(nanmean([MFRrat1, MFRA11],2), val_rat_mean(:, 1), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 12}, TanovaRat{1, 12}] =  anovan(nanmean([MFRrat2, MFRA12],2), val_rat_mean(:, 2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    % RAB
    [PanovaRat{1, 13}, TanovaRat{1, 13}] =  anovan(nanmean([MFRrat1, MFRA11, MFRB1],2), val_rat_mean(:, 1), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1,14}, TanovaRat{1, 14}] =  anovan(nanmean([MFRrat2, MFRA12, MFRB2],2), val_rat_mean(:, 2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    % RABA
    [PanovaRat{1, 15}, TanovaRat{1, 15}] =  anovan(nanmean([MFRrat1, MFRA11, MFRB1, MFRA21],2), val_rat_mean(:, 1), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 16}, TanovaRat{1, 16}] =  anovan(nanmean([MFRrat2, MFRA12, MFRB2, MFRA22],2), val_rat_mean(:, 2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    % R-Lik
    [PanovaRat{1, 17}, TanovaRat{1, 17}] =  anovan(nanmean([MFRrat1, MFRlikRat1],2), val_rat_mean(:, 1), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 18}, TanovaRat{1, 18}] =  anovan(nanmean([MFRrat2, MFRlikRat2],2), val_rat_mean(:, 2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    % AB
    [PanovaRat{1, 19}, TanovaRat{1, 19}] =  anovan(nanmean([MFRA11, MFRB1],2), val_rat_mean(:, 1), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 20}, TanovaRat{1, 20}] =  anovan(nanmean([MFRA12, MFRB2],2), val_rat_mean(:, 2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    % ABA
    [PanovaRat{1, 21}, TanovaRat{1, 21}] =  anovan(nanmean([MFRA11, MFRB1, MFRA21],2), val_rat_mean(:, 1), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 22}, TanovaRat{1, 22}] =  anovan(nanmean([MFRA12, MFRB2, MFRA22],2), val_rat_mean(:, 2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');    
    % AA
    [PanovaRat{1, 23}, TanovaRat{1, 23}] =  anovan(nanmean([MFRA11, MFRA21],2), val_rat_mean(:, 1), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 24}, TanovaRat{1, 24}] =  anovan(nanmean([MFRA12, MFRA22],2), val_rat_mean(:, 2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');    
    
    
    % now both both before and after lunch, I pool them the same way HIll
    % did it, so, doing the mean of both ratings before and after and
    % correlating with the pooled data
    
    
    % R, A1, B, A2, LikRat
    [PanovaRat{1, 25}, TanovaRat{1, 25}] =  anovan(nanmean([MFRrat1,MFRrat2],2), mean(val_rat_mean,2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 26}, TanovaRat{1, 26}] =  anovan(nanmean([MFRA11,MFRA12],2), mean(val_rat_mean,2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');    
    [PanovaRat{1, 27}, TanovaRat{1, 27}] =  anovan(nanmean([MFRB1,MFRB2],2), mean(val_rat_mean,2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 28}, TanovaRat{1, 28}] =  anovan(nanmean([MFRA21,MFRA22],2), mean(val_rat_mean,2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 29}, TanovaRat{1, 29}] =  anovan(nanmean([MFRlikRat1,MFRlikRat2],2), mean(val_rat_mean,2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
      
    % RA, RAB, RABA, R-Lik
    [PanovaRat{1, 30}, TanovaRat{1, 30}] =  anovan(nanmean([MFRrat1, MFRrat2, MFRA11, MFRA12],2), mean(val_rat_mean,2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 31}, TanovaRat{1, 31}] =  anovan(nanmean([MFRrat1, MFRrat2, MFRA11, MFRA12, MFRB1, MFRB2],2), mean(val_rat_mean,2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 32}, TanovaRat{1, 32}] =  anovan(nanmean([MFRrat1, MFRrat2, MFRA11, MFRA12, MFRB1, MFRB2, MFRA21, MFRA22],2), mean(val_rat_mean,2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 33}, TanovaRat{1, 33}] =  anovan(nanmean([MFRrat1,MFRrat2, MFRlikRat1,MFRlikRat2],2), mean(val_rat_mean,2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
   
    % AB, ABA, AA
    [PanovaRat{1, 34}, TanovaRat{1,34}]  =  anovan(nanmean([MFRA11, MFRA12, MFRB1, MFRB2],2), mean(val_rat_mean,2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 35}, TanovaRat{1, 35}] =  anovan(nanmean([MFRA11, MFRA12, MFRB1, MFRB2, MFRA21, MFRA22],2), mean(val_rat_mean,2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    [PanovaRat{1, 36}, TanovaRat{1, 36}] =  anovan(nanmean([MFRA11, MFRA12, MFRA21, MFRA22],2), mean(val_rat_mean,2), 'continuous',[1],'varnames',{'val_rat_in'},'display','off');
    
end



