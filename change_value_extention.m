function change_value_extention

% load form change value
% the question here is, how can quentufy how much the sellective satiation
% decreases the value of only the object selected. Is there any other type
% of correlation? For example, there could be a correlation between the 
% higher values decrease more. Or the products with same taste type
% decrease more as well. For example, how is the ranking change for sweet
% and salty?  

% Moreover 2 patients did not decrease the value of the specific selective
% satiation, even increased it. However they could have a 2 peak in ranking
% and they ate a different one, check this out. Moroever, this could also be



% This from here first analyses if any correlation exists correlating the
% value previous to the food and the change in value, expectinng to see a
% positive correlation, but this could be a spurious correlation and mean
% that they are bored of  the same taste. Check the nex analysis out bro!.

try
    cd('/media/Projects/Alex/Satiation Analysis');      %/Volumes/MTL/CS4 output')   %orginal form michael '/Volumes/Barnacle/Data/Satiation_data/Satiation_session_data'
catch
    cd('/media/Projects/Alex/');      %/Volumes/MTL/CS4 output')
end

DirName=dir;
h=1;
for a=1:length(DirName)
   
    if regexp(DirName(a).name, 'Sat_[0-9][0-9][0-9]')
        if isempty(regexp(DirName(a).name, 'Sat_039')) % don't analyze pat 39, he did not eat.
            
            [satcue{h},ranking{h}]=load_behavioural_data(DirName(a).name);
            rankingMatRat{h}=ranking{h}{2}(:,4)-ranking{h}{5}(:,4);
            rankingMatRank{h}=ranking{h}{3}(:,3)-ranking{h}{6}(:,3);
            rankingMatAFC{h}=ranking{h}{3}(:,4)-ranking{h}{6}(:,4);
            rankingMatTot{h}=ranking{h}{3}(:,5)-ranking{h}{6}(:,5);
            
            
            % do higher values decrease more or less than lower values?
%             figure;
%             plot(ranking{h}{2}(:,4),rankingMatRat{h}(:), 'o'); %plot 
%             title(['rating Pre vs change ration patient ', DirName(a).name])

            [spco1{h}, pco1{h}]=corr(ranking{h}{2}(:,4),rankingMatRat{h}(:), 'type', 'Spearman');  % 6 out of 9 sig positive correlation coeff, 3 after Boferonni
            [spco2{h}, pco2{h}]=corr(ranking{h}{3}(:,3),rankingMatRank{h}(:), 'type', 'Spearman'); % 3 out of 9 sig positive correlation coeff, 3 after Bonferoni
            [spco3{h}, pco3{h}]=corr(ranking{h}{3}(:,4),rankingMatAFC{h}(:), 'type', 'Spearman');  % 7 out of 9 sig positive correlation coeff, 5 after Bonferoni
            [spco4{h}, pco4{h}]=corr(ranking{h}{3}(:,5),rankingMatTot{h}(:), 'type', 'Spearman');  % 7 out of 9 sig positive correlation coeff, 5 after Bonferoni
            
            h=h+1;
        end
    end
end


% Here I extract the information of the change in value of the products
% that match the taste type of the Eaten product and the ones that not. 
% Seems to be that the diffenrence in rating of the same taste is higher
% than the ones that not. This makes sense if you think that they are
% bored of sweet or salty stuff. \

% Again be careful that the most rated one is the one they chose. Even
% though they use to give higher rates to things with the same rate, they
% and thus this would not make a sense in this case, this could not  be so.



% Finally  this is an analysis to check how many products get the highest
% ranking possible, this could be a problem to now what have they selected.



DirName=dir;
h=1;
for a=1:length(DirName)
   
    if regexp(DirName(a).name, 'Sat_[0-9][0-9][0-9]')
        if isempty(regexp(DirName(a).name, 'Sat_039')) % don't analyze pat 39, he did not eat.
            
            [satcue{h},ranking{h}]=load_behavioural_data(DirName(a).name);
            
           if  isempty(regexp(DirName(a).name, 'Sat_037')) %  pat 37, he ate the 2nd most prefered product, the first one was expired..... so the ERdnusse, number 16. 
               if  isempty(regexp(DirName(a).name, 'Sat_046')) %  pat 46, he also ate the second most prefered one, as Michael sorted by rank. 2 products have the highst rank in this patient, and while sorting he chose the wrong one, as the other one had higher values in all 3 other measures.
                   if  isempty(regexp(DirName(a).name, 'Sat_040')) %  pat 40, he also gave the highest rank to 2 products, but he luckily ate the most prefered one. 
                            [val,ind]= max(ranking{h}{3}(:,3));
                   else
                       ind=16; %erdnusee
                   end
               else
                   ind=17;
               end
           else
               ind=16; %erdnusee
           end   
           
           
           maxind{h}=find(ranking{h}{3}(:,3)==val);
           maxrank{h}= ranking{h}{3}(maxind{h},3);
           maxAFC{h} = ranking{h}{3}(maxind{h},4);
           maxTot{h}= ranking{h}{3}(maxind{h},5);
           maxRat{h} = ranking{h}{2}(maxind{h},4);
           
           
           if ind>10
               
               rankingMatRatSel{h}=ranking{h}{2}(11:20,4)-ranking{h}{5}(11:20,4);
               rankingMatRankSel{h}=ranking{h}{3}(11:20,3)-ranking{h}{6}(11:20,3);
               rankingMatAFCSel{h}=ranking{h}{3}(11:20,4)-ranking{h}{6}(11:20,4);
               rankingMatTotSel{h}=ranking{h}{3}(11:20,5)-ranking{h}{6}(11:20,5);
               
               
               rankingMatRatUn{h}=ranking{h}{2}(1:10,4)-ranking{h}{5}(1:10,4);
               rankingMatRankUn{h}=ranking{h}{3}(1:10,3)-ranking{h}{6}(1:10,3);
               rankingMatAFCUn{h}=ranking{h}{3}(1:10,4)-ranking{h}{6}(1:10,4);
               rankingMatTotUn{h}=ranking{h}{3}(1:10,5)-ranking{h}{6}(1:10,5);
               
           elseif ind<10
               
               rankingMatRatUn{h}=ranking{h}{2}(11:20,4)-ranking{h}{5}(11:20,4);
               rankingMatRankUn{h}=ranking{h}{3}(11:20,3)-ranking{h}{6}(11:20,3);
               rankingMatAFCUn{h}=ranking{h}{3}(11:20,4)-ranking{h}{6}(11:20,4);
               rankingMatTotUn{h}=ranking{h}{3}(11:20,5)-ranking{h}{6}(11:20,5);
               
               
               rankingMatRatSel{h}=ranking{h}{2}(1:10,4)-ranking{h}{5}(1:10,4);
               rankingMatRankSel{h}=ranking{h}{3}(1:10,3)-ranking{h}{6}(1:10,3);
               rankingMatAFCSel{h}=ranking{h}{3}(1:10,4)-ranking{h}{6}(1:10,4);
               rankingMatTotSel{h}=ranking{h}{3}(1:10,5)-ranking{h}{6}(1:10,5);
               
            end 
            
            
            h=h+1;
        end
    end
end

[~,p5{1}]=ttest2(cell2mat(rankingMatRankSel'), cell2mat(rankingMatRankUn'));
[~,p6{1}]=ttest2(cell2mat(rankingMatRatSel'), cell2mat(rankingMatRatUn'));
[~,p7{1}]=ttest2(cell2mat(rankingMatAFCSel'), cell2mat(rankingMatAFCUn'));
[~,p8{1}]=ttest2(cell2mat(rankingMatTotSel'), cell2mat(rankingMatTotUn'));

j=0;k=0;l=0;m=0;
for i=1:10

[~,p5{i+1}]=ttest2((rankingMatRankSel{i}), (rankingMatRankUn{i}));
[~,p6{i+1}]=ttest2((rankingMatRatSel{i}), (rankingMatRatUn{i}));
[~,p7{i+1}]=ttest2((rankingMatAFCSel{i}), (rankingMatAFCUn{i}));
[~,p8{i+1}]=ttest2((rankingMatTotSel{i}), (rankingMatTotUn{i}));
if mean(rankingMatRatSel{i})> mean(rankingMatRatUn{i})
    j=j+1;
end
if mean(rankingMatRankSel{i})> mean(rankingMatRankUn{i})
    k=k+1;
end
if mean(rankingMatAFCSel{i})> mean(rankingMatAFCUn{i})
    l=l+1;
end
if mean(rankingMatTotSel{i})> mean(rankingMatTotUn{i})
    m=m+1;
end
end




