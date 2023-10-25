clear all;

constructs={'Control','Mutant'};
reps=3;

bpdata=nan(50,8);
for c=1:length(constructs)
    fccln3nuc=[];
    fcbni1nuc=[]; 
    fccln3nonnuc=[];
    fcbni1nonnuc=[]; 
    for r=1:reps
        da=load([constructs{c} '_Rep_' num2str(r) '_colocalization_nuc_full_data.mat']);
        fccln3nuc=[fccln3nuc; da.fraccoloccln3nuc];
        fcbni1nuc=[fcbni1nuc; da.fraccolocbni1nuc];
        fccln3nonnuc=[fccln3nonnuc; da.fraccoloccln3nonnuc];
        fcbni1nonnuc=[fcbni1nonnuc; da.fraccolocbni1nonnuc];
        if r==1 
           pixelshift=da.pixelshift;
           radius=da.radius; 
        end
        clear da; 
    end
    bpdata(1:size(fccln3nuc,1),2*(c-1)+1)=fccln3nuc(:,1); % CLN3 fraction colocalized
    bpdata(1:size(fccln3nonnuc,1),2*(c-1)+2)=fccln3nonnuc(:,1); % CLN3 fraction colocalized with pixel shift
    bpdata(1:size(fcbni1nuc,1),2*(c-1)+5)=fcbni1nuc(:,1); % BNI1 fraction colocalized
    bpdata(1:size(fcbni1nonnuc,1),2*(c-1)+6)=fcbni1nonnuc(:,1); % BNI1 fraction colocalized with pixel shift
    %return
end

f=figure; 
f.Position=[100 100 500 250];

mycolor = [202 33 39; 202 33 39; 128 58 150; 128 58 150; 202 33 39; 202 33 39; 128 58 150; 128 58 150]./[255 255 255];
for i=1:8
    boxchart(i*ones(size(bpdata(:,i))), bpdata(:,i), 'BoxFaceColor', mycolor(i,:),'MarkerStyle','none'); hold on; 
    xc(:,i)=i*ones(1,size(bpdata,1));
end
swarmchart(xc,bpdata,20,mycolor,'filled','MarkerFaceAlpha',1.0,'MarkerEdgeAlpha',1.0);
ylabel('Fraction Colocalized')
set(gca,'Xtick',[1:1:8]);
set(gca,'Xticklabel',{'WT nuc','WT nonnuc','CLN3+/BNI1+ nuc','CLN3+/BNI1+ nonnuc','WT nuc','WT nonnuc','CLN3+/BNI1+ nuc','CLN3+/BNI1+ shift nonnuc'})
%xtickangle(45)

x=bpdata(:,1); y=bpdata(:,3); 
[p,b]=ranksum(x(~isnan(x)), y(~isnan(y)))
p
return
%print -painters -depsc 'invivo_colocalization_w_2rad_shift.eps'

%% No pixel shifts
clear xc; 

bpnsdata=bpdata(:,[1 3 5 7])

f=figure; 
f.Position=[100 100 500 250];

mycolor = [202 33 39; 128 58 150; 202 33 39; 128 58 150]./[255 255 255];
for i=1:4
    boxchart(i*ones(size(bpnsdata(:,i))), bpnsdata(:,i), 'BoxFaceColor', mycolor(i,:),'MarkerStyle','none'); hold on; 
    xc(:,i)=i*ones(1,size(bpnsdata,1));
end
swarmchart(xc,bpnsdata,40,mycolor,'filled','MarkerFaceAlpha',1.0,'MarkerEdgeAlpha',1.0);
ylabel('Fraction Colocalized')
set(gca,'Xtick',[1:1:4]);
set(gca,'Xticklabel',{'WT','CLN3+/BNI1+','WT','CLN3+/BNI1+'})
xtickangle(45)

%print -painters -depsc 'invivo_colocalization.eps'

%kstest: test for normality
% 1 implies not normal
% If data is normal use two-sample t-test (ttest2) rather than wilcoxon rank sum (ranksum)

%[pcln3,hcln3]=ranksum(bpnsdata(:,1),bpnsdata(:,2))
x=bpnsdata(:,1); y=bpnsdata(:,2); 
[pcln3,hcln3]=ranksum(x(~isnan(x)), y(~isnan(y)))
%[h,p]=kstest2(bpnsdata(:,1),bpnsdata(:,2))

%[pbni1,hbni1]=ranksum(bpnsdata(:,3),bpnsdata(:,4))
x=bpnsdata(:,3); y=bpnsdata(:,4); 
[pbni1,hbni1]=ranksum(x(~isnan(x)), y(~isnan(y)))
%[h,p]=kstest2(bpnsdata(:,3),bpnsdata(:,4))

pbni1-pcln3


