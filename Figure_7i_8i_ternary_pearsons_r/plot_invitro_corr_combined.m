clear all;

[da,txt]=xlsread('220831_combined_invitro_correlation_data.xlsx','Sheet1');

systems={'BC Del','BC Sim','BS Del','BS Sim','Bmut C Del','Bmut C Sim'};
%systems={'C3C5 Del','C3C5 Sim','B3B5 Del','B3B5 Sim','Bmut3Bmut5 Del','Bmut3Bmut5 Sim'};


systemlist=txt(2:end,1);
pearsonr=da(:,3);

f=figure; 
%f.Position = [100 100 600 150];
f.Position = [100 100 800 150];

mycolor = [202 33 39; 128 58 150; 202 33 39; 128 58 150; 202 33 39; 128 58 150]./[255 255 255];
mycolor2 = [220 125 132; 166 123 182; 220 125 132; 166 123 182; 220 125 132; 166 123 182]./[255 255 255];

count=0;
for s=1:length(systems)
    count=count+1;
    pos=find(strcmp(systemlist,systems{s}));
    pearsonr(pos)
    apearsonr(count)=mean(pearsonr(pos));
    sem(count)=std(pearsonr(pos))/sqrt(length(pos));

    for p=1:length(pos)
        plot(count+randsample([-0.1:0.01:0.1],1),pearsonr(pos(p)),'o','MarkerEdgeColor',mycolor(s,:),'MarkerFaceColor',[1 1 1]); hold on; 
    end
    %return
    clear pos;
    
    if rem(s,2)==0
        count=count+1;
    end
end

%f=figure; 
%f.Position = [100 100 600 150];
%mycolor = [202 33 39; 128 58 150; 202 33 39; 128 58 150; 202 33 39; 128 58 150]./[255 255 255];
%mycolor2 = [220 125 132; 166 123 182; 220 125 132; 166 123 182; 220 125 132; 166 123 182]./[255 255 255];
%%bp=bar(avgrgroup2); hold on;

count=0;
for b = 1:length(apearsonr)
    h=bar(b,apearsonr(b)); hold on;
    if apearsonr(b)~=0
        count=count+1;
        set(h,'EdgeColor',mycolor(count,:));
        set(h,'FaceColor',mycolor2(count,:));
    end
end

er = errorbar(1:1:length(apearsonr),apearsonr,sem);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

ylabel('Pearsons r');

%print -painters -depsc 'invitro_pearsonsr_combined_data_with_points.eps'

%% Add in simultaneous time titration
timvec={'0 min','5 min','15 min','30 min','1 hr','2 hr','3 hr','4 hr'}
rep(1,:)=[0.94	0.95	0.96	0.96	0.94	0.95	0.93	0.91];
rep(2,:)=[0.9	0.92	0.94	0.96	0.94	0.96	0.93	0.94];
rep(3,:)=[0.9	0.94	0.95	0.94	0.95	0.93	0.95	0.93];

apearsonrtime=mean(rep);
semtime=std(rep)/sqrt(size(rep,1));

f2=figure; 
f2.Position = [100 100 800 150];

for i=1:size(rep,2)
    for p=1:size(rep,1)
        plot(i+randsample([-0.1:0.01:0.1],1),rep(p,i),'o','MarkerEdgeColor',mycolor(2,:),'MarkerFaceColor',[1 1 1]); hold on; 
    end
end
%return

for b = 1:length(apearsonrtime)
    h=bar(b,apearsonrtime(b)); hold on;
    if apearsonrtime(b)~=0
        set(h,'EdgeColor',mycolor(2,:));
        set(h,'FaceColor',mycolor2(2,:));
    end
end

er = errorbar(1:1:length(apearsonrtime),apearsonrtime,semtime);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

ylabel('Pearsons r');

set(gca,'xtick',1:1:length(timvec))
set(gca,'XTickLabel',timvec)
