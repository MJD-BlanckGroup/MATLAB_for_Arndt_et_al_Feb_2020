%% discovery-based program correlating gene expression to KM survival distinction
% inputs: RNASeq is the exome file from cBioPortal, CBIO contains the read
% groups depicted by KM analysis, polif_apop_genes and immune_genes
% contains lists of well-researched effector and immune genes for sample
% comparison

function HNSC(RNASeq,CBIO,prolif_apop_genes, immune_genes)
for i=2
    gothroughallgenesincolumn(RNASeq,i,CBIO,prolif_apop_genes,immune_genes)
end
end

% need to automatically index through thousands of genes
function gothroughallgenesincolumn(RNASeq,i,CBIO, pa, immune)
% "i" refers to column
for r= 2:20532
    gene= RNASeq(r,1);
    sigtesting(gene,i,RNASeq,CBIO, pa, immune);
end
end

%determing signifance in gene expression between two groups
function sigtesting(gene,column, RNASeq, CBIO, pa, immune)
 
e_row=RNASeq(:,1);e_col=RNASeq(1,:);
indapop=strfind(e_row,gene);
 
for i =1:numel(indapop)
    if indapop{i,1}==1
        genelocation=i;
    end
    if exist('genelocation') == 1
        break
    end
end

hcolval=[]; lcolval=[];
z= CBIO(:,column); 
Upper_percentile=[];Lower_percentile=[];
 
if column==1
    for o=1:numel(z)
        if o<17
            hstr= z(o,1);
            fileind = strfind(e_col,hstr);
             for i=1:numel(fileind)
                if fileind{1,i}==1
                hcolval(end+1)=i;
                break
                end
                end
        end
                      
        if o>17
            lstr = z(o,1);
            fileind2=strfind(e_col,lstr);
            for i=1:numel(fileind2)
                if fileind2{1,i}==1
                    lcolval(end+1)= i;
                    break
                end
            end
        end
    end
    
        for u=1:numel(lcolval)
            Upper_percentile(end+1)= RNASeq(genelocation, lcolval(u));
        end
        for l=1:numel(hcolval)
            Lower_percentile(end+1)= RNASeq(genelocation, hcolval(l));
        end
    if isempty(Upper_percentile) || isempty(Lower_percentile)
        disp('no rows returned')
    else
    [~,p] = ttest2(Upper_percentile,Lower_percentile);
    end
end


if column==2 
    for o=1:numel(z)
        if o<17
            hstr= z(o,1);
            fileind= strfind(e_col,hstr);
            for i =1:numel(fileind)
                if fileind{1,i}==1
                  hcolval(end+1)= i;
                  break  
                end
            end
        end
        if o>17
            lstr = z(o,1);
            fileind2=strfind(e_col,lstr);
            for i=1:numel(fileind2)
                if fileind2{1,i}==1
                    lcolval(end+1)= i;
                    break
                end
            end
        end
    end
    for u=1:numel(hcolval)
        Upper_percentile(end+1) = RNASeq(genelocation, hcolval(u));
    end
    for l=1:numel(lcolval)
        Lower_percentile(end+1) = RNASeq(genelocation, lcolval(l));
       
    end
    if isempty(Upper_percentile) || isempty(Lower_percentile)
        disp('no rows returned')
    else
    [~,p] = ttest2(Upper_percentile,Lower_percentile);
    end
end

if p<0.05 %alpha level of 0.05 used

if column==1
    Attribute = 'J1-3 01 & DPB1*03:01';
elseif column==2 
    Attribute = 'J2-1 01 & DPB1*02:01';
end

%checking for proliferation genes
proliffiltered=pa(:,1);
apoptofiltered=pa(:,2);
cp=contains(proliffiltered,gene);ca=contains(apoptofiltered,gene);ci=contains(immune,gene);

if (sum(cp) == 1 && (p<0.05))
     PROLIFERATION='Y'; APOPTOSIS_EFF='N';IMMUNE='N';
%checking for apoptosis effector genes
elseif (sum(ca) == 1 && (p<0.05))
     APOPTOSIS_EFF='Y'; PROLIFERATION='N'; IMMUNE='N';
elseif (sum(ci) == 1 && (p<0.05))
    IMMUNE='Y'; APOPTOSIS_EFF='N'; PROLIFERATION='N';
elseif p<0.05
     PROLIFERATION='N'; APOPTOSIS_EFF='N'; IMMUNE='N';
end

%calculating top&bottom half averages
Top_50_Avg= mean(Upper_percentile, 'all');
Bot_50_Avg = mean(Lower_percentile, 'all');

results = [gene, Attribute, p, Top_50_Avg, Bot_50_Avg, PROLIFERATION, APOPTOSIS_EFF, IMMUNE];
fileID = fopen('HN_SCC.txt','a+');
fprintf(fileID,'%s %s %f %f %f %s %s %s\n', results);
fclose(fileID);

end

end

