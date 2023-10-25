clear all;

% Load in sequences
[h,seqs]=fastaread('221025_Whi3_RNA_Sequences_V2.txt')

%pairs={'cln3_bni1.txt','spa2_bni1.txt','cln3_bni1mut.txt','cln3_cln3.txt','bni1_bni1.txt','bni1mut_bni1mut.txt'};
pairs={'cln3_cln3.txt','cln3_bni1mut.txt','cln3_bni1.txt','bni1mut_bni1mut.txt','bni1_bni1.txt','spa2_bni1.txt'};

for p=1:length(pairs)
    m = 1;
    d = fopen(pairs{p});
    % Skip 7 headerlines
    for n=1:5
        tline = fgetl(d);
    end
    while ischar(tline)
        tline = fgetl(d);
        matches{m} = tline;
        m = m+1;
        % Skip next two lines
        for n=1:3
            tline = fgetl(d);
        end
    end
    fclose(d);
    %matches
    for i=1:length(matches)
        tmp=strsplit(matches{i},' ');
        rnast(i,1)=str2num(tmp{5});
        rnast(i,2)=str2num(tmp{9});
        rnast(i,3)=str2num(tmp{2});
        %return
        clear tmp; 
    end

    % For the ones that are the same RNA
    if p>3
        for i=1:length(rnast)
            tmp1=rnast(i,1);
            tmp2=rnast(i,2);
            if tmp1>tmp2
               rnast(i,1)=tmp2;
               rnast(i,2)=tmp1;
            end
        end
        urnast=unique(rnast,'rows');
    else
        urnast=unique(rnast,'rows');
    end
    numcompsites(p)=length(urnast)

    % Plot sites
    tmp=strsplit(pairs{p},'_')
    rna1=tmp{1};
    tmp2=strsplit(tmp{2},'.');
    rna2=tmp2{1};
    pos1=find(strcmp(h,rna1)==1);
    pos2=find(strcmp(h,rna2)==1);
    figure;
    plot([1 length(seqs{pos1})],[2 2],'-k'); hold on; 
    plot([1 length(seqs{pos2})],[1 1],'-k');
    for i=1:length(urnast)
        rectangle('Position',[urnast(i,1) 1.75 urnast(i,3) 0.5]); 
        rectangle('Position',[urnast(i,2) 0.75 urnast(i,3) 0.5]); 
    end
    %return
    clear rnast; clear sitelength; clear urnast; clear matches; 
end

figure;
bar(numcompsites)
ylabel('Number of Complementary Sites');
set(gca,'xticklabel',{'cln3-cln3','cln3-bni1mut','cln3-bni1','bni1mut-bni1mut','bni1-bni1','bni1-spa2'});
xtickangle(45);



