clear all;

% Load in sequences
[h,seqs]=fastaread('221025_Whi3_RNA_Sequences_V2.txt')

%rna1={'cln3','spa2','cln3','bni1'};
%rna2={'bni1','bni1','cln3','bni1'};
rna1={'spa2','bni1','cln3','cln3'};
rna2={'bni1','bni1','bni1','cln3'};
mycolor=[1 0 0; 0 1 0; 0 0 1; 0 0 0];

for p=1:length(rna1)
    m = 1;
    d = fopen([rna1{p} '_' rna2{p} '.txt']);
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

    rna1shapeinfo=load([rna1{p} 'shapeinfo.mat']);
    rna2shapeinfo=load([rna2{p} 'shapeinfo.mat']);
    rna1seq=rna1shapeinfo.finalseq;
    rna1shape=rna1shapeinfo.shapevals;
    rna2seq=rna2shapeinfo.finalseq;
    rna2shape=rna2shapeinfo.shapevals;
    pos=find(rna1shape==-999);
    rna1shape(pos)=NaN;
    pos=find(rna2shape==-999);
    rna2shape(pos)=NaN;

    for i=1:length(urnast)
        % If have shape data for the sequence then calculate the mean shape values for this complementary pairing
        if urnast(i,1)+urnast(i,3)-1<=length(rna1seq) & urnast(i,2)+urnast(i,3)-1<=length(rna2seq)
           rna1val=rna1shape(urnast(i,1):urnast(i,1)+urnast(i,3)-1);
           rna2val=rna2shape(urnast(i,2):urnast(i,2)+urnast(i,3)-1);
           mshapeval=mean([rna1val rna2val]);
           urnast(i,4)=mshapeval;
           %return
           clear rna1val; clear rna2val; clear mshapeval; 
        else
           urnast(i,4)=NaN;
	end
    end
    figure(29); 
    histogram(urnast(:,4),[-0.5:0.1:2.5],'facecolor',mycolor(p,:),'facealpha',0.5); hold on; 
    ylabel('Number of Complementary Sites');
    xlabel('Mean Shape');
    meanshape(p)=nanmean(urnast(:,4));
    legend('bni1-spa2','bni1-bni1','cln3-bni1','cln3-cln3');
    urnast(:,4)

    %% Plot sites
    %pos1=find(strcmp(h,rna1)==1);
    %pos2=find(strcmp(h,rna2)==1);
    %figure;
    %plot([1 length(seqs{pos1})],[2 2],'-k'); hold on; 
    %plot([1 length(seqs{pos2})],[1 1],'-k');
    %for i=1:length(urnast)
    %    rectangle('Position',[urnast(i,1) 1.75 urnast(i,3) 0.5]); 
    %    rectangle('Position',[urnast(i,2) 0.75 urnast(i,3) 0.5]); 
    %end

    %return
    clear rnast; clear sitelength; clear urnast; clear matches; 
end
