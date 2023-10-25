clear all;

clear all; 

constructs={'230705_Andrew_Bni1mut','220927_Andrew_Bni1'};
%nt=[6689 6689];
nt=[1 1];
mycolor=[60 130 63; 60 130 63]/255;

%figure;
for c=1:length(constructs)
    %% Load in phase separation data
    da=load([constructs{c} '_LLPS_uM_nM.txt']);
    da2=load([constructs{c} '_noLLPS_uM_nM.txt']);


    f1=figure(2090);
    f1.Position = [100 100 250 250];
    if c==2
        plot((da(:,1)),(0.001.*da(:,2)),'o','markerfacecolor',[204 113 118]/255,'markeredgecolor',[204 113 118]/255,'markersize',8,'linewidth',0.1); hold on;
        plot((da2(:,1)), (0.001.*da2(:,2)),'o','markerfacecolor',[0.4 0.4 0.4],'markeredgecolor',[0.4 0.4 0.4],'markersize',8,'linewidth',0.1);
    elseif c==1
        plot((da(:,1)),(0.001.*da(:,2)),'o','markerfacecolor',[183 34 37]/255,'markeredgecolor',[183 34 37]/255,'markersize',8,'linewidth',0.1); hold on;5
        plot((da2(:,1)), (0.001.*da2(:,2)),'ok','markerfacecolor','k','markersize',8,'linewidth',0.1);
    end
    %plot([-7 4],[-7 4],'k');
    %xlim([-7 4]); ylim([-7 4])
    xlim([0 1]); ylim([0 10^-2])
    xlabel('[Whi3 Molecules]');
    ylabel('[RNA Molecules]');
    %set(gca,'YTick',[-6:2:2]);
    axis square
    
    f1=figure(2091);
    f1.Position = [100 100 250 250];
    if c==2
        plot(log10(da(:,1)),log10(0.001.*da(:,2)),'.','markerfacecolor',[204 113 118]/255,'markeredgecolor',[204 113 118]/255,'markersize',5,'linewidth',0.1); hold on;
        plot(log10(da2(:,1)), log10(0.001.*da2(:,2)),'.','markerfacecolor',[0.4 0.4 0.4],'markeredgecolor',[0.4 0.4 0.4],'markersize',5,'linewidth',0.1);
    elseif c==1
        plot(log10(da(:,1)),log10(0.001.*da(:,2)),'o','markerfacecolor',[183 34 37]/255,'markeredgecolor',[183 34 37]/255,'markersize',4,'linewidth',0.1); hold on;
        plot(log10(da2(:,1)), log10(0.001.*da2(:,2)),'ok','markerfacecolor','k','markersize',4,'linewidth',0.1);
    end
    %plot([-7 4],[-7 4],'k');
    %xlim([-7 4]); ylim([-7 4])
    xlim([-3 3]); ylim([-6 0])
    xlabel('[Whi3 Molecules]');
    ylabel('[RNA Molecules]');
    set(gca,'YTick',[-6:2:2]);
    axis square
    
    if c==2
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
        
        % Get goodness of symmetry by getting R^2 of points fit to ellipse
        x2=log10(xPS);
        y2=log10(yPS*0.001*nt(c));
        pos1=find(x2==-Inf);
        pos2=find(y2==-Inf);
        x2([pos1; pos2])=[];
        y2([pos1; pos2])=[];
        data_points=[x2 y2]';
        [theta_guaranteed] = compute_guaranteedellipse_estimates(data_points);

        % Plot data points
        x = data_points';
        n = length(x);
        %figure(11)
        %plot(x(:,1),x(:,2),'o','color',mycolor(c,:));
        hold on
        %return

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
        set(h, 'Color', mycolor(c,:),'linewidth',2);
        title([])
        xlabel('[Whi3 Molecules] (\muM)')
        ylabel('[RNA Molecules] (\muM)')
        
    elseif c==1
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

        x2=log10(xPS);
        y2=log10(yPS*0.001*nt(c));
        pos1=find(x2==-Inf);
        pos2=find(y2==-Inf);
        x2([pos1; pos2])=[];
        y2([pos1; pos2])=[];
        k=boundary(x2,y2,0.5);        
        %plot(x2(k),y2(k),'-r'); hold on;
        polyData=polyshape(x2(k),y2(k));
        plot(polyData)

        %% Create data for text file
        mycell=cell(max(length(da(:,1)),length(da2(:,1))),4);
        mycell(1:length(da(:,1)),1)=num2cell(da(:,1)) % [Whi3] phase separation
        mycell(1:length(da(:,1)),2)=num2cell(0.001.*da(:,2)) % [RNA] phase separation
        mycell(1:length(da2(:,1)),3)=num2cell(da2(:,1)) % [Whi3] no phase separation
        mycell(1:length(da2(:,1)),4)=num2cell(0.001.*da2(:,2)) % [RNA] no phase separation


        % Create table to add to excel
        mytable=cell2table(mycell,'VariableNames',{'[Whi3 phase separation molecule]','[RNA phase separation molecule]','[Whi3 no phase separation molecule]','[RNA no phase separation molecule]'});
    
        % Write to excel
        writetable(mytable,'../../../2023_10/Source_Data/Lin_etal_Source_Data_File.xlsx','Sheet','4d BNI1mut')
        %return
    end
        
    
    %plot(da(:,1),da(:,2)*0.001*nt(c),'o','markerfacecolor',mycolor(c,:),'markeredgecolor','k','markersize',20); hold on;
    %plot(da2(:,1),da2(:,2)*0.001*nt(c),'s','markerfacecolor',mycolor(c,:),'markeredgecolor','k','markersize',12); hold on;
end
%set(gca,'yscale','log');
%set(gca,'xscale','log');
