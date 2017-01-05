function [mama, mama2, papa, papa2, m1, m2, p1, p2, p_Val, p_Sal, lm_rat]=change_RT_salience


type=0;  % type 0 means Pearsons matrix of correlation, if 1 means Spearmans

f={};
% mama={}; mama2={}; papa={}; papa2={};
mama=[]; mama2=[]; papa=[]; papa2=[]; diffRTout=[]; diffSalout=[]; diffValout=[];diffValabsout=[];diffSalabsout=[]; lm_rat={};
try
    cd('/media/Projects/Alex/Satiation Analysis');      %/Volumes/MTL/CS4 output')   %orginal form michael '/Volumes/Barnacle/Data/Satiation_data/Satiation_session_data'
catch
    cd('/media/Projects/Alex/');      %/Volumes/MTL/CS4 output')
end

DirName= dir;
for a = 1:length(DirName)
    if regexp(DirName(a).name, 'Sat_+[0-9][0-9][0-9]')
        if isempty(regexp(DirName(a).name, 'Sat_039'))
            f{end+1,1} = [DirName(a).name];
            [satcue, ranking, events] = load_behavioural_data(DirName(a).name);
            [mama(end+1,:), papa(end+1,:), mama2(end+1,:), papa2(end+1,:), diffRT,  diffVal, diffValabs, diffSal, diffSalabs, lm_rat{end+1}]= findcorr_change_RT_salience(satcue, ranking,events);

            diffRTout=[diffRTout;diffRT(:)];
            diffSalout=[diffSalout;diffSal(:)];
            diffValout=[diffValout; diffVal(:)]; 
            diffValabsout=[diffValabsout; diffValabs(:)];
            diffSalabsout=[diffSalabsout; diffSalabs(:)];
        end
    end
end

[m1,m2,p1,p2]=total_corr(diffRTout,diffValout,diffValabsout,diffSalout, diffSalabsout);

for i=1:length(lm_rat)
gg(i)=lm_rat{i}.Coefficients.Estimate(3);
hh(i)=lm_rat{i}.Coefficients.Estimate(2);
end
[~,p_Val]=ttest(hh);
[~,p_Sal]=ttest(gg);


end



function [mama, papa, mama2, papa2, diffRT,  diffVal, diffValabs, diffSal, diffSalabs, lm_rat]= findcorr_change_RT_salience(satcue, ranking, events)

[RT] = extract_reaction_times(events);

for i=1:20
    
    ind3=find(cell2mat(satcue{2}(:,2))==i);
    ind4=find(cell2mat(satcue{5}(:,2))==i);
    
    diffRT(i)= mean(RT{2}(ind3)) -mean(RT{5}(ind4));          % RT for rating vector, 60 entries
    diffVal(i)= ranking{2}(i,4) - ranking{5}(i,4);
    diffValabs(i)= abs(ranking{2}(i,4) - ranking{5}(i,4));
    diffSal(i)= abs(ranking{2}(i,4)) - abs(ranking{5}(i,4));
    diffSalabs(i)= abs(abs(ranking{2}(i,4)) - abs(ranking{5}(i,4)));
end

[mama(1), papa(1)]=corr(diffRT', diffVal');
[mama(2), papa(2)]=corr(diffRT', diffValabs');
[mama(3), papa(3)]=corr(diffRT', diffSal');
[mama(4), papa(4)]=corr(diffRT', diffSalabs');

[mama2(1), papa2(1)]=corr(diffRT', diffVal', 'type', 'Spearman');
[mama2(2), papa2(2)]=corr(diffRT', diffValabs', 'type', 'Spearman');
[mama2(3), papa2(3)]=corr(diffRT', diffSal', 'type', 'Spearman');
[mama2(4), papa2(4)]=corr(diffRT', diffSalabs', 'type', 'Spearman');

 
tbl1=table(zscore(diffRT'), zscore(diffVal'), zscore(diffSal'),'VariableNames',{'RT','Val','Sal'});
lm_rat =fitlm(tbl1,'RT ~ Val+ Sal');

end

function [m1,m2,p1,p2]=total_corr(diffRT,diffVal,diffValabs,diffSal, diffSalabs)

[m1(1), p1(1)]=corr(diffRT, diffVal);
[m1(2), p1(2)]=corr(diffRT, diffValabs);
[m1(3), p1(3)]=corr(diffRT, diffSal);
[m1(4), p1(4)]=corr(diffRT, diffSalabs);

[m2(1), p2(1)]=corr(diffRT, diffVal, 'type', 'Spearman');
[m2(2), p2(2)]=corr(diffRT, diffValabs, 'type', 'Spearman');
[m2(3), p2(3)]=corr(diffRT, diffSal, 'type', 'Spearman');
[m2(4), p2(4)]=corr(diffRT, diffSalabs, 'type', 'Spearman');

    
end