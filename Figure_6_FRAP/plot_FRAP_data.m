clear all;

%% Whi3
system={'Whi3_Trace_Cln3_221011','Whi3_Trace_BNI1','Whi3_Trace_Spa2','Whi3_Trace_BNI1_Mut'};
mycolor=[144 95 167; 60 130 63; 47 71 153; 154 216 154]/255;

f=figure;
f.Position=[100 100 800 600]
for s=1:length(system)
    [da,txt]=xlsread('220926_Compiled_FRAP_Traces_Avg_SEM_Max_update.xlsx',system{s});
    norm_net_max(s)=da(1,7);
    norm_net_max_sem(s)=da(1,9);
    time=da(:,1);
    norm_mean=da(:,2);
    norm_sem=da(:,3);
    subplot(2,2,s)
    myeplot=errorbar(time,norm_mean,norm_sem,'-o','color',mycolor(s,:),'markerfacecolor',[1 1 1],'markersize',4,'LineWidth',1.5); hold on;
    myeplot.CapSize=0;
    xlim([-5, 185])
    ylim([-0.05,1.1])
    xlabel('Time (s)')
    ylabel('Intensity (AU)')
    clear time; clear norm_mean; clear norm_sem; clear myeplot; 
end
%print -painters -depsc 'FRAP_traces_Whi3_Protein.eps'
%norm_net_max
%return

f=figure; 
f.Position = [100 100 250 150];
mycolor=[144 95 167; 60 130 63; 47 71 153; 154 216 154]/255;

count=0;
for b = 1:length(norm_net_max)
    h=bar(b,norm_net_max(b)); hold on;
    if norm_net_max(b)~=0
        count=count+1;
        set(h,'EdgeColor',mycolor(count,:));
        set(h,'FaceColor',mycolor(count,:));
    end
end

er = errorbar(1:1:length(norm_net_max),norm_net_max,norm_net_max_sem);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
clear f; 

%return
%% RNAs
system={'Whi3_RNA_Cln3_v2','Whi3_RNA_BNI1','Whi3_RNA_Spa2','Whi3_RNA_BNI1Mut'};

f=figure;
f.Position=[100 100 800 600]
for s=1:length(system)
    [da,txt]=xlsread('220926_Compiled_FRAP_Traces_Avg_SEM_Max_update.xlsx',system{s});
    norm_net_max(s)=da(1,7);
    norm_net_max_sem(s)=da(1,9);
    time=da(:,1);
    norm_mean=da(:,2);
    norm_sem=da(:,3);
    subplot(2,2,s)
    myeplot=errorbar(time,norm_mean,norm_sem,'-o','color',mycolor(s,:),'markerfacecolor',[1 1 1],'markersize',4,'LineWidth',1.5); hold on;
    myeplot.CapSize=0;
    xlim([-5, 185])
    ylim([-0.05,1.1])
    xlabel('Time (s)')
    ylabel('Intensity (AU)')
    clear time; clear norm_mean; clear norm_sem; clear myeplot;
    %return
end
%print -painters -depsc 'FRAP_traces_RNAs.eps'

f=figure; 
f.Position = [100 100 250 150];
mycolor=[144 95 167; 60 130 63; 47 71 153; 154 216 154]/255;

count=0;
for b = 1:length(norm_net_max)
    h=bar(b,norm_net_max(b)); hold on;
    if norm_net_max(b)~=0
        count=count+1;
        set(h,'EdgeColor',mycolor(count,:));
        set(h,'FaceColor',mycolor(count,:));
    end
end

er = errorbar(1:1:length(norm_net_max),norm_net_max,norm_net_max_sem);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
clear f; 

%% Whi3 & RNAs
system={'Whi3_Trace_Cln3_221011','Whi3_Trace_BNI1','Whi3_Trace_Spa2','Whi3_Trace_BNI1_Mut','Whi3_RNA_Cln3_v2','Whi3_RNA_BNI1','Whi3_RNA_Spa2','Whi3_RNA_BNI1Mut'};

f=figure; 
f.Position = [100 100 600 150];
mycolor=[144 95 167; 60 130 63; 47 71 153; 154 216 154; 144 95 167; 60 130 63; 47 71 153; 154 216 154]/255;

net_max_trials=nan(length(system),4);
for s=1:length(system)
    [da,txt]=xlsread('220926_Compiled_FRAP_Traces_Avg_SEM_Max_update.xlsx',system{s});
    norm_net_max(s)=da(1,7);
    norm_net_max_sem(s)=da(1,9);
    if s<3 | s>3
        net_max_trials(s,1:3)=da(1:3,10);
    else
        net_max_trials(s,1:4)=da(1:4,10);
    end
    %return
end

count=0;
for b = 1:length(norm_net_max)
    h=bar(b,norm_net_max(b)); hold on;
    if norm_net_max(b)~=0
        count=count+1;
        set(h,'EdgeColor',mycolor(count,:));
        set(h,'FaceColor',mycolor(count,:));
    end
end

er = errorbar(1:1:length(norm_net_max),norm_net_max,norm_net_max_sem);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel('Net Max FRAP')

for s=1:length(system)
    for p=1:4
        plot(s+randsample([-0.1:0.01:0.1],1),net_max_trials(s,p),'o','MarkerEdgeColor',mycolor(s,:),'MarkerFaceColor',[1 1 1]); hold on; 
    end
end
%print -painters -depsc 'Net_Max_FRAP_bars.eps'
