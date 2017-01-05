function change_value

f={};
try
    cd('/media/Projects/Alex/Satiation Analysis');      %/Volumes/MTL/CS4 output')   %orginal form michael '/Volumes/Barnacle/Data/Satiation_data/Satiation_session_data'
catch
    cd('/media/Projects/Alex/');      %/Volumes/MTL/CS4 output')
end

DirName= dir;
j=1;
for a = 1:length(DirName)
    if regexp(DirName(a).name, 'Sat_+[0-9][0-9][0-9]')
        if isempty(regexp(DirName(a).name, 'Sat_039'))
            f{end+1,1} = [DirName(a).name];
            [satcue, ranking]=load_behavioural_data(DirName(a).name);
            
            
            if  isempty(regexp(DirName(a).name, 'Sat_037')) %  pat 37, he ate the 2nd most prefered product, the first one was expired..... so the ERdnusse, number 16.
                if  isempty(regexp(DirName(a).name, 'Sat_046')) %  pat 46, he also ate the second most prefered one, as Michael sorted by rank. 2 products have the highst rank in this patient, and while sorting he chose the wrong one, as the other one had higher values in all 3 other measures.
                    if  isempty(regexp(DirName(a).name, 'Sat_040')) %  pat 40, he also gave the highest rank to 2 products, but he luckily ate the most prefered one.
                        [val,indPrefered]= max(ranking{3}(:,3));
                    else
                        indPrefered=16; %erdnusee
                    end
                else
                    indPrefered=17;
                end
            else
                indPrefered=16; %erdnusee
            end
            
         
            ratPre(j)=  ranking{2}(indPrefered,4);
            ratPost(j)= ranking{5}(indPrefered,4);
            
            rankPre(j)= ranking{3}(indPrefered,3);
            rankPost(j)=ranking{6}(indPrefered,3);
            
            AFCPre(j)=  ranking{3}(indPrefered,4);
            AFCPost(j)= ranking{6}(indPrefered,4);
            
            TotPre(j)=  ranking{3}(indPrefered,5);
            TotPost(j)= ranking{6}(indPrefered,5);
            
            j=j+1;
        end
    end
end
[~,p1]=ttest(rankPre,rankPost);
[~,p2]=ttest(ratPre,ratPost);
[~,p3]=ttest(AFCPre,AFCPost);
[~,p4]=ttest(TotPre,TotPost);

% meanTot = mean((TotPost+1900)./(TotPre+1900));    stdTot = std((TotPost+1900)./(TotPre+1900));
% meanACF = mean((AFCPost+950)./(AFCPre+950));    stdACF = std((AFCPost+950)./(AFCPre+950)); 
% meanrat = mean((ratPost+300)./(ratPre+300)); stdrat = std((ratPost+300)./(ratPre+300));
% meanrank = mean(rankPost./rankPre);    stdrank = std(rankPost./rankPre);

meanTot = mean((TotPost)./(TotPre));    stdTot = std((TotPost)./(TotPre));
meanACF = mean((AFCPost)./(AFCPre));    stdACF = std((AFCPost)./(AFCPre)); 
meanrat = mean((ratPost)./(ratPre)); stdrat = std((ratPost)./(ratPre));
meanrank = mean(rankPost./rankPre);    stdrank = std(rankPost./rankPre);

% d=[(rankPre./rankPre)', (rankPost./rankPre)',(ratPre./ratPre)',(ratPost./ratPre)',(AFCPre./AFCPre)',(AFCPost./AFCPre)',(TotPre./TotPre)', (TotPost./TotPre)'];

meanMatrix=[1, meanrank; 1, meanrat; 1,meanACF; 1, meanTot];
stdMatrix= [0,  stdrank; 0,  stdrat; 0, stdACF; 0,  stdTot];

for k=1:4
h=0; 
eval(['p=p' num2str(k) ';']);
while p<0.05
    p=p*10;
    h=h+1;
end
    bridgeMarix(:,:,k)=ones(2)*(h+1);
end


barweb(meanMatrix, stdMatrix);

% 

%% By taste type 












% ,bridgeMarix
% barWhiskerBridge(meanMatrix,stdMatrix)

% 
% DirName= dir;
% j=1;
% for a = 1:length(DirName)
%     if regexp(DirName(a).name, 'Sat_+[0-9][0-9][0-9]')
%         if isempty(regexp(DirName(a).name, 'Sat_039'))
%             f{end+1,1} = [DirName(a).name];
%             [satcue, ranking]=load_behavioural_data(DirName(a).name);
%             [ ~,indPrefered]=max(ranking{3}(:,3));
%             
%             ratPre(j)=  ranking{2}(indPrefered,4); 
%             ratPost(j)= ranking{5}(indPrefered,4); 
%             
%             rankPre(j)= ranking{3}(indPrefered,3); 
%             rankPost(j)=ranking{6}(indPrefered,3);
%    
%             AFCPre(j)=  ranking{3}(indPrefered,4); 
%             AFCPost(j)= ranking{6}(indPrefered,4);
%             
%             TotPre(j)=  ranking{3}(indPrefered,5);
%             TotPost(j)= ranking{6}(indPrefered,5);
% 
%             j=j+1;
%         end 
%     end
% end
% 
% [~,p1]=ttest(TotPre,TotPost);
% [~,p2]=ttest(AFCPre,AFCPost);
% [~,p3]=ttest(rankPre,rankPost);
% [~,p4]=ttest(ratPre,ratPost);
% 
% meanTot = mean(TotPost./TotPre);
% meanACF = mean(AFCPost./AFCPre); 
% meanrank = mean(rankPost./rankPre);
% meanrat = mean(ratPost./ratPre);


% % 
% % d = randi(10, 20, 3);
% figure(1)
% boxplot(d)
% yt = get(gca, 'YTick');
% axis([xlim    0  ceil(max(yt)*1.2)])
% set(gca, 'xtick', 1:3);
% xt = get(gca, 'XTick');
% hold on
% plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k')
% hold off

