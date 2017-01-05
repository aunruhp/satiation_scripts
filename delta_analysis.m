function delta_analysis(paradigmBlock, shiftNum, posShift)


% delta coefficient: from Padoa-Schioppa 2008, Range-adapting
% representation in OFC. Journal of Neuroscience

% In this paper, they analyse what they call partial adaptation. This is
% the FAST trial by trial adaptation. They define DELTA as the proportion
% of neural activity change from one trial to another. For examople, take 
% trials with chosen value =4. Look at previous trial's chosen value. And
% then if the previous is higher, it shoud be now lower. Calculate and 
% normalize the difference according to mean firing rate of this neuron.
% Finally, define success as delta>0, and binomial test. Also, calculate
% mean delta and SD and do one-way t-test. 
% paradigmBlock is a number equals 2,3,5 or 6 depending on which block you want
% to test


if nargin<3
    posShift=0;
end

if paradigmBlock~=6
    loc_var=floor(paradigmBlock/2);
else 
    loc_var= 2;
end
   
%%

[satcue, ranking] = load_behavioural_data;

%%



% extract information about covariates, ie, mean values of stimulus over time to correlate.



val_rat_rat= nan(length(satcue{2}(:,14)),2);
val_rat_2AFC=nan(length(satcue{2}(:,14)),2);
val_rat_rank=nan(length(satcue{2}(:,14)),2);

for i=1:20
    ind=find(cell2mat(satcue{2}(:,2))==i);
    valRat=  ranking{2}(i,4); % ranking
    val2AFC= ranking{3}(i,4);
    valRank= ranking{3}(i,3);
    
    val_rat_rat(ind,1)= valRat;
    val_rat_2AFC(ind,1)= val2AFC;
    val_rat_rank(ind,1)= valRank;
    
    ind2=find(cell2mat(satcue{5}(:,2))==i);
    valRat=  ranking{5}(i,4); % ranking
    val2AFC= ranking{6}(i,4);
    valRank= ranking{6}(i,3);
    val_rat_rat(ind2,2)= valRat;
    val_rat_2AFC(ind2,2)= val2AFC;
    val_rat_rank(ind2,2)= valRank;
    
    
end

%scaling values
val_rat_rat= val_rat_rat./3;             % so it goes from -100 to +100
val_rat_2AFC= (val_rat_2AFC ./9.5) -100; % so it goes from -100 to +100





classif=  [(cell2mat(satcue{paradigmBlock}(:,14))), val_rat_rank(:,loc_var), val_rat_rat(:,loc_var), val_rat_2AFC(:,loc_var)];  


for kk=[-posShift:-1, 1:shiftNum]
    
    for hh = 1:4



    small=cell(20,1);
    big=cell(20,1);
    equal=cell(20,1);

    for i=1:20
        index= find(cell2mat(satcue{paradigmBlock}(:,2))==i);

        for j=1:length(index)
     
            if index(j)- kk <= 0 || index(j)- kk > length(cell2mat(satcue{paradigmBlock}(:,14)))  % if at the beginning or the end, imposible to subtract, so do nothing!
                     
            
            elseif classif(index(j)-kk,hh) > classif(index(j),hh)
                small{i}(end+1)=  cell2mat(satcue{paradigmBlock}(index(j),14));
            elseif classif(index(j)-kk,hh) < classif(index(j),hh)
                big{i}(end+1)  =  cell2mat(satcue{paradigmBlock}(index(j),14));
            elseif classif(index(j)-kk,hh) == classif(index(j),hh)
                equal{i}(end+1)=  cell2mat(satcue{paradigmBlock}(index(j),14));
            end

        end
    end


    ind=[];
    for i=1:20
        if length(small{i})>0 && length(small{i})<3
           if length(big{i})>0 && length(big{i})<3
              ind(end+1)=i;
           end
        end
    end
    kkk= kk+posShift +1; %shift so that no negatice nor zero index is included
    eval(['delta' num2str(kkk), num2str(hh) '=[];'])

    
    h=1;
    for i=ind
        normalize_val = abs(ranking{3}(i,4) ./9.5 -100);  % other possibility   abs(ranking{2}(i,4) ./3);
        eval(['delta' num2str(kkk), num2str(hh) '(h,1:2)= (small{i} - big{i}) / normalize_val ;'] )  %difference betwwen when previous bigger vs smaller, normalized by normal value of the object. 
        
        h=h+1;
    end



    eval(['delta =delta' num2str(kkk), num2str(hh) ';'])

    binotest(kkk,hh)= 1- binocdf(sum(sum(delta>0)), numel(delta), 0.05);  % if no influence whatsoever, results all the same. error of 0.05 accepted, but if higher, then reject stability and accept influence
    binoneg (kkk,hh) = 1- binocdf(sum(sum(delta<0)), numel(delta), 0.05);
    
    comptest(kkk,hh)= 1- binocdf(sum(sum(delta>0)), numel(delta), 0.1);  % if no influence whatsoever, results all the same. error of 0.5 accepted, but if higher, then reject stability and accept influence
    compneg (kkk,hh) = 1- binocdf(sum(sum(delta<0)), numel(delta), 0.1); 
%     small 
%     big
%     equal
  


   
    end
end
 binotest
 comptest
 binoneg
 compneg

end 
