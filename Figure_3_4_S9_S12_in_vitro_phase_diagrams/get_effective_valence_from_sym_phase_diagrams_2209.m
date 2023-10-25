clear all;

constructs={'Cln3','Bni1','Spa2'};
stval=[1 50 200];
enval=[50 200 900];
skip=[1 1 1];
%reval=NaN(4,100);
reval=NaN(3,150);
%mycolor=[106 59 155; 170 103 43; 41 56 153]/255;
mycolor=[144 95 167; 60 130 63; 47 71 153]/255;

sheetnames={'Figure 3 CLN3', 'Figure 3 BNI1', 'Figure 3 SPA2'}


for c=1:length(constructs)
    %% Load in phase separation data
    da=load(['220927_Andrew_' constructs{c} '_LLPS_uM_nM.txt']);
    da2=load(['220927_Andrew_' constructs{c} '_noLLPS_uM_nM.txt']);
    
    %% Find edges of 2-phase and 1-phase regime as these are the values that will be assessed to determine how symmetric the data is
    
    % First find lower boundary of phase separated points
    uwhi3con=unique(da2(:,1));
    for i=1:length(uwhi3con)
        tmp=find(da(:,1)==uwhi3con(i));
        if isempty(tmp)==0
            tmp2(i)=min(da(tmp,2));
        else
            tmp2(i)=0;
        end
        clear tmp; 
    end
    
    xPS=[uwhi3con];
    yPS=[tmp2'];

    clear uwhi3con; clear tmp2;
    
    uwhi3con=unique(da2(:,2));
    for i=1:length(uwhi3con)
        tmp=find(da(:,2)==uwhi3con(i));
        if isempty(tmp)==0
            tmp2(i)=min(da(tmp,1));
            %return
        else
            tmp2(i)=0;
        end
        clear tmp; 
    end

    xPS=[xPS; tmp2'];
    yPS=[yPS; uwhi3con];
    clear tmp2; clear uwhi3con; 
    
    % Next find boundary of non phase separated points
    xnoPS=[];
    ynoPS=[];
    pos=find(yPS>=xPS);
    for i=1:length(pos)
        pos2=find(da2(:,2)==yPS(pos(i)) & da2(:,1)<xPS(pos(i)));
        ynoPS=[ynoPS da2(pos2,2)'];
        xnoPS=[xnoPS da2(pos2,1)'];
        clear pos2; 
    end  
    clear pos; 
    
    pos=find(xPS>=yPS);
    for i=1:length(pos)
        pos2=find(da2(:,2)<yPS(pos(i)) & da2(:,1)==xPS(pos(i)));
        ynoPS=[ynoPS da2(pos2,2)'];
        xnoPS=[xnoPS da2(pos2,1)'];
        clear pos2; 
    end
    clear pos; 
  
    % Test plot
    %figure; plot(log10(xPS), log10(yPS),'o');
    %hold on; plot(log10(xnoPS), log10(ynoPS),'o');
    %hold on; plot(log10(da(:,1)),log10(da(:,2)),'*');
    %hold on; plot(log10(da2(:,1)),log10(da2(:,2)),'*');
    %set(gca,'Yscale','log');
    %set(gca,'Xscale','log');
    %return
    
    %% Now rescale the y-axis and check symmetry by seeing if rescaling leads to a low overlap between the area of the 2-phase and 1-phase regimes
    count=0;
    for s=[stval(c):skip(c):enval(c)]
        count=count+1;
        
        % Test plot
        %figure; plot(log10(xPS), log10(s*0.001.*yPS),'o');
        %hold on; plot(log10(xnoPS), log10(s*0.001.*ynoPS),'o');
        %hold on; plot(log10(da(:,1)),log10(s*0.001.*da(:,2)),'*');
        %hold on; plot(log10(da2(:,1)),log10(s*0.001.*da2(:,2)),'*');
        
        % Get boundary of 1-phase regime
        nx=log10(xnoPS);
        ny=log10(s*0.001.*ynoPS);
        k=boundary(nx',ny');
        %plot((nx(k)),(ny(k)),'-k'); hold on;
        polynoPS=polyshape((nx(k)),(ny(k)));
        %plot(polynoPS);
        clear k; 

        % Get overlap between this boundary and the boundary of symmetric fit
        x2=log10([xPS; s*0.001.*yPS]);
        y2=log10([s*0.001.*yPS; xPS]);
        pos1=find(x2==-Inf);
        pos2=find(y2==-Inf);
        x2([pos1; pos2])=[];
        y2([pos1; pos2])=[];
        k=boundary(x2,y2);
        %plot(x2(k),y2(k),'-r'); hold on;
        polyPS=polyshape(x2(k),y2(k));
        %plot(polyPS);
        %plot([-6 4],[-6 4],'k');
        clear k; 
        
        polyout=intersect(polynoPS,polyPS);
        myarea(count)=area(polyout);
        %myarea(count)=area(polyout)/area(polynoPS); % would go between 0
        %and 1 but the problem is this depends on how many points in the
        %1-phase regime so still can't compare across phase diagrams
        psarea(count)=area(polyPS);
        nopsarea(count)=area(polynoPS);
        %return
        clear nx; clear ny; clear polynoPS; 
        clear x2; clear y2; clear pos1; clear pos2; clear polyPS; 
        clear polyout; 
    end
    
    %% Plot the best symmetrized data
    srange=[stval(c):skip(c):enval(c)];
    
    mycolors=[0 0 1; 0.5 0.5 0; 0.5 0 0.5; 1 0 0];
    
    %figure(29)
    %subplot(1,3,c)
    f3=figure;
    %f3.Position = [100 100 250 300];
    bar(myarea,'facecolor',mycolor(c,:));
    title([constructs{c}],'fontsize',22);
    if c<3
        set(gca,'XTick',[1:5:length(srange)]);
        set(gca,'XTickLabel',srange(1:5:end));
    else
        set(gca,'XTick',[1:20:length(srange)]);
        set(gca,'XTickLabel',srange(1:20:end));
    end
    xtickangle(90)
    ylabel('Overlap Area','fontsize',20)
    xlabel('Effective Valence','fontsize',20)
    %myarea'
    %return
    
    m=min(myarea);
    pos=find(m==myarea);
    s=srange(pos);
    
    perchange=0.05;
    fmin=m-m*perchange;
    fmax=m+m*perchange;
    pos2=find(myarea>=fmin & myarea<=fmax);
    
    effval(c)=s;
    reval(c,1:length(pos2))=srange(pos2);
    
    
    % Plot non-scaled
    f1=figure; 
    %f1.Position = [100 100 250 300];
    f1.Position = [100 100 250 250];
    plot(log10(da(:,1)),log10(0.001.*da(:,2)),'o','markerfacecolor',[183 34 37]/255,'markeredgecolor',[183 34 37]/255,'markersize',4,'linewidth',0.1); hold on;
    plot(log10(da2(:,1)), log10(0.001.*da2(:,2)),'ok','markerfacecolor','k','markersize',4,'linewidth',0.1);
    plot([-7 4],[-7 4],'k');
    xlim([-7 4]); ylim([-7 4])
    xlabel('[Whi3 Molecules]');
    ylabel('[RNA Molecules]');
    set(gca,'YTick',[-6:2:2]);
    axis square
    %return
    
    % Plot scaled no fits
    f2=figure(200+c); 
    %f2.Position = [100 100 250 300];
    f2.Position = [100 100 250 250];
    plot(log10(da(:,1)),log10(s*0.001.*da(:,2)),'o','markerfacecolor',[183 34 37]/255,'markeredgecolor',[183 34 37]/255,'markersize',4,'linewidth',0.1); hold on;
    plot(log10(da2(:,1)), log10(s*0.001.*da2(:,2)),'ok','markerfacecolor','k','markersize',4,'linewidth',0.1);
    plot([-7 4],[-7 4],'k');
    xlim([-7 4]); ylim([-7 4])
    xlabel('[Whi3 Stickers]');
    ylabel('[RNA Stickers]');

    %% Create data for text file
    mycell=cell(max(length(da(:,1)),length(da2(:,1))),8);
    mycell(1:length(da(:,1)),1)=num2cell(da(:,1)) % [Whi3] phase separation
    mycell(1:length(da(:,1)),2)=num2cell(0.001.*da(:,2)) % [RNA] phase separation
    mycell(1:length(da2(:,1)),3)=num2cell(da2(:,1)) % [Whi3] no phase separation
    mycell(1:length(da2(:,1)),4)=num2cell(0.001.*da2(:,2)) % [RNA] no phase separation
    mycell(1:length(da(:,1)),5)=num2cell(da(:,1)) % [Whi3] phase separation
    mycell(1:length(da(:,1)),6)=num2cell(s*0.001.*da(:,2)) % [RNA] phase separation
    mycell(1:length(da2(:,1)),7)=num2cell(da2(:,1)) % [Whi3] no phase separation
    mycell(1:length(da2(:,1)),8)=num2cell(s*0.001.*da2(:,2)) % [RNA] no phase separation

    % Create table to add to excel
    mytable=cell2table(mycell,'VariableNames',{'[Whi3 phase separation molecule]','[RNA phase separation molecule]','[Whi3 no phase separation molecule]','[RNA no phase separation molecule]','[Whi3 phase separation sticker]','[RNA phase separation sticker]','[Whi3 no phase separation sticker]','[RNA no phase separation sticker]'});
    
    % Write to excel
    %writetable(mytable,'../../../2023_10/Source_Data/Lin_etal_Source_Data_File.xlsx','Sheet',sheetnames{c})
    %return
    
    figure;
    plot(log10(xPS), log10(s*0.001.*yPS),'ob'); hold on;
    plot(log10(s*0.001.*yPS),log10(xPS),'*b'); hold on;
    plot(log10(xnoPS), log10(s*0.001.*ynoPS),'or');
    plot(log10(da(:,1)),log10(s*0.001.*da(:,2)),'.b');
    %plot(log10(da2(:,1)),log10(s*0.001.*da2(:,2)),'*');
    
    nx=log10(xnoPS);
    ny=log10(s*0.001.*ynoPS);
    k=boundary(nx',ny');
    polynoPS=polyshape((nx(k)),(ny(k)));
    plot(polynoPS,'FaceColor','red','FaceAlpha',0.1);
    clear k; 
    
    x2=log10([xPS; s*0.001.*yPS]);
    y2=log10([s*0.001.*yPS; xPS]);
    pos1=find(x2==-Inf);
    pos2=find(y2==-Inf);
    x2([pos1; pos2])=[];
    y2([pos1; pos2])=[];
    k=boundary(x2,y2);
    polyPS=polyshape(x2(k),y2(k));
    plot(polyPS,'FaceColor','blue','FaceAlpha',0.1);
    plot([-6 4],[-6 4],'k');
    
    title([constructs{c} ' Apparent Valence: ' num2str(s)]);
    xlabel('[Whi3 Stickers]');
    ylabel('[RNA Stickers]');
    
    % Get goodness of symmetry by getting R^2 of points fit to ellipse
    data_points=[x2 y2]';
    [theta_guaranteed] = compute_guaranteedellipse_estimates(data_points);
    
    % Plot data points
    x = data_points';
    n = length(x);
    figure('Color',[1 1 1])
    plot(x(:,1),x(:,2),'b.');
    hold on
    
    % determine data range
    minX = min(min(x(:,1))) - 1;
    minY = min(min(x(:,2))) - 1;
    maxX = max(max(x(:,1))) + 1;
    maxY = max(max(x(:,2)))  + 1;
    
    % plot the guaranteed ellipse fit
    a1 = theta_guaranteed(1); b1 = theta_guaranteed(2); c1 = theta_guaranteed(3);
    d1 = theta_guaranteed(4); e1 = theta_guaranteed(5); f1 = theta_guaranteed(6);
    fh = @(x,y) (a1*x.^2 + b1*x.*y + c1*y.^2 + d1*x + e1*y + f1);
    h = ezplot(fh,[minX maxX minY maxY]);
    
    myval=x;
    for i=1:size(myval,1)
        % Solve for y
        myfun=@(y) (a1*myval(i,1)^2 + b1*myval(i,1)*y + c1*y^2 + d1*myval(i,1) + e1*y + f1);
        y0=myval(i,2);
        yp(i) = fsolve(myfun,y0);
        plot(myval(i,1),yp(i),'*'); 
        clear x0; clear myfun; 
    end
    
    figure(200+c)
    ezplot(fh,[minX maxX minY maxY]);
    xlabel('[Whi3 Stickers]');
    ylabel('[RNA Stickers]');
    title([constructs{c} ' Apparent Valence: ' num2str(s)]);
    set(gca,'YTick',[-6:2:2])
    axis square
    
    SSR=sum((myval(:,2)-yp').^2);
    SST=sum((myval(:,2)-mean(myval(:,2))).^2);
    r2(c)=1-(SSR/SST)
    
    f10=figure(200)
    f10.Position=[100 100 250 250]
    h=ezplot(fh,[minX maxX minY maxY]); hold on; 
    set(h, 'Color', mycolor(c,:),'linewidth',2);
    plot(log10(xPS), log10(s*0.001.*yPS),'*','markeredgecolor',mycolor(c,:)); hold on;
    plot(log10(da(:,1)),log10(s*0.001.*da(:,2)),'.','color',mycolor(c,:));
    xlabel('[Whi3 Stickers]');
    ylabel('[RNA Stickers]');
    title([])

    % Find concentration of RNA sticker at 5uM Whi3
    tmppos=find(xPS==5);
    rnaconc(c)=s*0.001.*yPS(tmppos)
    %tmpnostick(c)=yPS(tmppos)
    
    whi3conclist=[0.05 0.1 0.5 1 5 10 100];
    for w=1:length(whi3conclist)
        tmppos=find(xPS==whi3conclist(w));
        if length(tmppos)==1
            rnaconcmat(w,c)=s*0.001.*yPS(tmppos);
        else
            tmpval=min(yPS(tmppos));
            rnaconcmat(w,c)=s*0.001.*tmpval; 
            clear tmpval; 
        end
    end
    
    %return
    clear srange; clear m; clear pos; clear s; 
    clear nx; clear ny; clear polynoPS; clear myarea; 
    clear fmin; clear fmax; clear pos2; 
    clear x2; clear y2; clear pos1; clear pos2; clear k; clear polyPS; 
    clear da; clear da2; clear xPS; clear yPS; clear xnoPS; clear ynoPS; 
    clear data_points; clear theta_guaranteed; clear x; clear n; 
    clear minX; clear minY; clear maxX; clear maxY; 
    clear a1; clear b1; clear c1; clear d1; clear e1; clear f1; 
    clear yp; clear SSR; clear SST; 
end

%return
% Three RNAs
nt=[1594 6689 10908];

%mycolors=[1 0 0; 0.5 0 0.5; 0.5 0.5 0; 0 0 1];
%mycolors=[41 56 153; 170 103 43; 106 59 155]/255;
mycolors=[47 71 153; 60 130 63; 144 95 167]/255;

f4=figure;
f4.Position = [100 100 250 300];
h=boxplot(reval', 'positions', nt, 'labels', nt); hold on;
bx = findobj(gca,'Tag','Box');
for b = 1:length(bx)
    bx(b).Color = mycolors(b,:);
end
bx = findobj(gca,'Tag','Median');
for b = 1:length(bx)
    bx(b).Color = mycolors(b,:);
end
set(h,'LineWidth',2);

%% Create data for text file
mycell2=cell(size(reval,1),3);
mycell2(1:length(reval(1,:)),1)=num2cell(reval(1,:));
mycell2(1:length(reval(2,:)),2)=num2cell(reval(2,:)');
mycell2(1:length(reval(3,:)),3)=num2cell(reval(3,:)');

% Create table to add to excel
mytable2=cell2table(mycell2,'VariableNames',{'CLN3 5% Overlap Change','BNI1 5% Overlap Change','SPA2 5% Overlap Change'});

% Write to excel
%writetable(mytable2,'../../../2023_10/Source_Data/Lin_etal_Source_Data_File.xlsx','Sheet','Figure 4b')
%return

mycolors2=[106 59 155; 170 103 43; 41 56 153]/255;
for i=1:length(nt)
    plot(nt(i),effval(i),'o','markeredgecolor',mycolors2(i,:));
end
ylabel('Apparent Valency');
xlabel('RNA Length');
ylim([0 800])

% Plot csats by Whi3 concentration and RNA

yvalues = whi3conclist;
xvalues = {'CLN3','BNI1','SPA2'};
figure;
h = heatmap(xvalues,yvalues,1000*rnaconcmat);
%h.XLabel = '';
h.YLabel = 'Whi3 Concentration (\muM)';

% Plot log-log of RNA length and apparent valence
ntall=[];
effvalall=[];
for i=1:length(nt)
    tmp=reval(i,:);
    tmp2 = tmp(~isnan(tmp))
    ntall=[ntall nt(i)*ones(1,length(tmp2))];
    effvalall=[effvalall tmp2];
end

figure; plot(log10(ntall),log10(effvalall),'o'); hold on;
p=polyfit(log10(ntall),log10(effvalall),1)
x=3:0.1:4.5;
y=p(1).*x+p(2);
plot(x,y,'-')

f5=figure;
f5.Position = [100 100 250 300]
mycolors2=[106 59 155; 170 103 43; 41 56 153]/255;
for i=1:length(nt)
    plot(log10(nt(i)),log10(effval(i)),'o','markerfacecolor',mycolors2(i,:),'markeredgecolor',mycolors2(i,:),'markersize',12); hold on;
end
p=polyfit(log10(nt),log10(effval),1)
x=2.95:0.1:4.5;
y=p(1).*x+p(2);
plot(x,y,'-k');
xlabel('Nucleotides')
ylabel('Apparent Valence')
xlim([2.95 4.2])
set(gca,'xtick',[3 4])
set(gca,'XTickLabel',{'10^3','10^4'})
set(gca,'ytick',[1,2,3])
set(gca,'YTickLabel',{'10^1','10^2','10^3'})
text(3.1,3.4,['Slope = ' num2str(p(1))])


