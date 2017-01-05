function responses_to_sort_rating


cd '/media/Projects/Alex/New Analysis/saturday_afternoon/'
folderdir=dir;

numPos=0;
numNeg=0;
numunits=0;
units_region=zeros(length(folderdir),1);
units_region2=zeros(length(folderdir),1);

for l=1:length(folderdir)
    if ~isempty(regexp(folderdir(l).name, 'Hcoeff_responses_'))
        load (folderdir(l).name)

        % extract information from the matafiles cells
        units_region2(l)=length(response_judge1);
        posResponsesStim=[];negResponsesStim=[];
        for unit=1:length(response_judge1)
            
                posResponsesStim(unit)= response_judge1{unit}{2}(1);
                negResponsesStim(unit)= response_judge1{unit}{2}(2);
                numunits=numunits+1;
                units_region(l)=units_region(l)+1;
 
        
       
        end 
        numPos=numPos+sum(posResponsesStim>1);
        numNeg=numNeg+sum(negResponsesStim>1);
        
    end
end


%Binomial test

numPos= length(left_value);
if numPos <= 190/2
    p_bino_l = sum(binopdf([1:numPos], numunits, 0.5));
    if p_bino_l< 0.05
        disp('This patient has a significat bias due to the position on screen when shown at choice time (left vs right), two-tail binomial test')
    end
elseif numPos > 190/2
    p_bino_r = sum(binopdf([numPos:190], 190, 0.5));
    if p_bino_r<0.05
        disp('This patient has a significat bias due to the position on screen when shown at choice time (left vs right), two-tail binomial test')
    end
end



cd '/media/Projects/Alex/New Analysis/saturday_afternoon/'
folderdir=dir;

numPos=0;
numNeg=0;
numunits=0;
units_region=zeros(length(folderdir),1);
for l=1:length(folderdir)
    if ~isempty(regexp(folderdir(l).name, 'Hcoeff_responses_'))
        load (folderdir(l).name)
       
        % extract information from the matafiles cells

        for unit=1:length(response_rating5)
            for win=1:length(response_rating5{1})
            posResponsesStim(unit)= response_rating5{unit}{2}(1);
            negResponsesStim(unit)= response_rating5{unit}{2}(2);
            numunits=numunits+1;
            units_region(l)=units_region(l)+1;
        end
        
        numPos=numPos+sum(posResponsesStim>1);
        numNeg=numNeg+sum(negResponsesStim>1);
        
        end 
    end
end



%% 



cd '/media/Projects/Alex/New Analysis/saturday_afternoon/'
folderdir=dir;



posResponsesStim=[];negResponsesStim=[];

for l=1:length(folderdir)
    if ~isempty(regexp(folderdir(l).name, 'Hcoeff_responses_'))
        load (folderdir(l).name)

        % extract information from the matafiles cells
        
        toanalyze=response_rating5;
        
        for unit=1:length(toanalyze)
            
                posResponsesStim = [posResponsesStim, toanalyze{unit}{2}(1)];
                negResponsesStim = [negResponsesStim, toanalyze{unit}{2}(2)];

 
        
       
        end 

    end
end
numPos=sum(posResponsesStim>=1);
numNeg=sum(negResponsesStim>=1);
numunit=length(negResponsesStim);


