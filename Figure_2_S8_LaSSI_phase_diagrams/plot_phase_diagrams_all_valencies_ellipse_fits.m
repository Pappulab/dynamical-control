clear all;

mol1=8;

figure;
for t=1:3
    for m2=2:8
        mol2=m2;
        if t==1
            da=importdata(['HetOnly_' num2str(mol1) '_' num2str(mol2) '.txt']);
        elseif t==2
            da=importdata(['HomA_' num2str(mol1) '_' num2str(mol2) '.txt']);
        elseif t==3
            da=importdata([num2str(mol1) '_' num2str(mol2) '_HetAB_HomA_HomB.txt']);
        end

        co1=da.data(:,1);
        co2=da.data(:,2);
        cdil1=da.data(:,7);
        cdil2=da.data(:,8);
        cden1=da.data(:,3);
        cden2=da.data(:,4);

        % Remove data that has dense phase concentration less than dilute phase
        % Or dilute phase concentration greater than starting concentration
        % Also add the dilute and dense must be over an order of magnitude
        % difference and starting and dilute must be 1.5 fold different
        pos=find(cden1-cdil1<0 | cden2-cdil2<0 | cdil1>co1 | cdil2>co2 | cden2./cdil2<10 | cden1./cdil1<10 | co2./cdil2<1.5 | co1./cdil1<1.5 | cdil2==0 | cdil1 ==0 | cden1./co1<3.5 | cden2./co2<3.5);

        cdil1(pos)=[];
        cdil2(pos)=[];
        cden1(pos)=[];
        cden2(pos)=[];
        co1(pos)=[];
        co2(pos)=[];

        %% Create data for excel file
        mycell=cell(length(cdil1),4);
        mycell(:,1)=num2cell(cdil2) % Whi3 dilute sticker conc
        mycell(:,2)=num2cell(cdil1) % RNA dilute sticker conc
        mycell(:,3)=num2cell(cden2) % Whi3 dense sticker conc
        mycell(:,4)=num2cell(cden1) % RNA dilute sticker conc

        % Create table to add to excel
        mytable=cell2table(mycell,'VariableNames',{'[B dilute sticker]','[A dilute sticker]','[B dense sticker]','[A dense sticker]'});
        
        % Write to excel
        if t==1
            sheetnames=['Supp Figure 8 HetAB_' num2str(mol1) '_' num2str(mol2)]
        elseif t==2
            sheetnames=['Supp Figure 8 HetAB_HomoA_' num2str(mol1) '_' num2str(mol2)]
        elseif t==3 
            sheetnames=['Supp Figure 8 HetAB_HomoAB_' num2str(mol1) '_' num2str(mol2)]
        end
        %writetable(mytable,'../../../2023_10/Source_Data/Lin_etal_Source_Data_File.xlsx','Sheet',sheetnames)
        %return

        % Potential outlier removal
        %xy=[cdil2, cdil1];
        %[result,rmv] = rmoutliers(xy,"movmedian",5);
        %pos=find(rmv==1);
        %cdil1(pos)=[];
        %cdil2(pos)=[];
        %cden1(pos)=[];
        %cden2(pos)=[];
        %co1(pos)=[];
        %co2(pos)=[];
        %return

        subplot(3,7,7*(t-1)+(m2-1))
        plot(log10(cdil2),log10(cdil1),'o','markerfacecolor',[183 211 168]/255,'markeredgecolor',[183 211 168]/255,'markersize',4,'linewidth',0.1); hold on;
        %plot(log10(cdil1),log10(cdil2),'*b'); hold on;
        plot(log10(cden2),log10(cden1),'o','markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0],'markersize',4,'linewidth',0.1); hold on;
        plot([-8 0],[-8 0],'-k');
        ylim([-8 0]); xlim([-8 0]);

        %% Get shape of 2-phase regime
        % Actual data from simulations
        x2=[log10(cdil2); log10(cden2)];
        y2=[log10(cdil1); log10(cden1)];
        pos1=find(x2==-Inf);
        pos2=find(y2==-Inf);
        x2([pos1; pos2])=[];
        y2([pos1; pos2])=[];

        % For homotypic phase separation csat - create limits
        % For homotypic interactions should still phase separate if the
        % other molecule goes to zero so using the minimum cdil value to
        % create this limit
        if t>1
            idx=find(cdil2==min(cdil2));
            tmpxdil=[-10:0.5:log10(cdil2(idx))-0.5];
            tmpydil=log10(cdil1(idx))*ones(1,length(tmpxdil));
            plot(tmpxdil,tmpydil,'*','markerfacecolor',[183 211 168]/255,'markeredgecolor',[183 211 168]/255,'markersize',4,'linewidth',0.1); hold on;

            idx=find(cden2==min(cden2));
            tmpxden=[-10:0.5:log10(cden2(idx))-0.5];
            tmpyden=log10(cden1(idx))*ones(1,length(tmpxden));
            plot(tmpxden,tmpyden,'*','markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0],'markersize',4,'linewidth',0.1); hold on;

            y2=[y2; tmpydil'; tmpyden'];
            x2=[x2; tmpxdil'; tmpxden'];
            %return
        end
        if t==3 & mol2>2
            idx=find(cdil1==min(cdil1));
            tmpydil=[-10:0.5:log10(cdil1(idx))-0.5];
            tmpxdil=log10(cdil2(idx))*ones(1,length(tmpydil));
            plot(tmpxdil,tmpydil,'*','markerfacecolor',[183 211 168]/255,'markeredgecolor',[183 211 168]/255,'markersize',4,'linewidth',0.1); hold on;

            idx=find(cden1==min(cden1));
            tmpyden=[-10:0.5:log10(cden1(idx))-0.5];
            tmpxden=log10(cden2(idx))*ones(1,length(tmpyden));
            plot(tmpxden,tmpyden,'*','markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0],'markersize',4,'linewidth',0.1); hold on;

            y2=[y2; tmpydil'; tmpyden'];
            x2=[x2; tmpxdil'; tmpxden'];
        end

        % Convert points to clockwise
        cx = mean(x2);
        cy = mean(y2);
        a = atan2(y2 - cy, x2 - cx);
        [~, sortIdx] = sort(a);

        sx2=x2(sortIdx); sy2=y2(sortIdx);
        k=boundary(sx2,sy2,0.5);        
        %plot(x2(k),y2(k),'-r'); hold on;
        polyData=polyshape(sx2(k),sy2(k));
        plot(polyData)
        clear x2; clear y2; 

        %% Fit dilute arm of phase diagram to ellipse
        % If scaled by valence of stickers
        x2=[log10(cdil2); log10(cdil1)];
        y2=[log10(cdil1); log10(cdil2)];
        pos1=find(x2==-Inf);
        pos2=find(y2==-Inf);
        x2([pos1; pos2])=[];
        y2([pos1; pos2])=[];
        %return
        % if not rescaled by valence
        %x2=[log10(cdil2./mol2); log10(cdil1./mol1)];
        %y2=[log10(cdil1./mol1); log10(cdil2./mol2)];

        % For homotypic phase separation csat - create limits
        if t>1
            idx=find(cdil2==min(cdil2));
            tmpxdil=[-10:0.5:log10(cdil2(idx))-0.5];
            tmpydil=log10(cdil1(idx))*ones(1,length(tmpxdil));
            plot(tmpxdil,tmpydil,'*','markerfacecolor',[183 211 168]/255,'markeredgecolor',[183 211 168]/255,'markersize',4,'linewidth',0.1); hold on;

            y2=[y2; tmpydil'; tmpxdil'];
            x2=[x2; tmpxdil'; tmpydil'];
            %return
        end
        if t==3 & mol2>2
            idx=find(cdil1==min(cdil1));
            tmpydil=[-10:0.5:log10(cdil1(idx))-0.5];
            tmpxdil=log10(cdil2(idx))*ones(1,length(tmpydil));
            plot(tmpxdil,tmpydil,'*','markerfacecolor',[183 211 168]/255,'markeredgecolor',[183 211 168]/255,'markersize',4,'linewidth',0.1); hold on;

            y2=[y2; tmpydil'; tmpxdil'];
            x2=[x2; tmpxdil'; tmpydil'];
        end


        % Get goodness of symmetry by getting R^2 of points fit to ellipse
        data_points=[x2 y2]';
        [theta_guaranteed] = compute_guaranteedellipse_estimates(data_points);

        % Plot data points
        x = data_points';
        n = length(x);
        %figure('Color',[1 1 1])
        %plot(x(:,1),x(:,2),'b.');
        %hold on

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
            %return
            clear y0; clear myfun; 
        end
        SSR=sum((myval(:,2)-yp').^2);
        SST=sum((myval(:,2)-mean(myval(:,2))).^2);
        r2=1-(SSR/SST);
        if t==1
            text(-7.5,-0.5,['R^2 = ' num2str(round(r2,2))]);
        else
            text(-7.5,-7.5,['R^2 = ' num2str(round(r2,2))]);
        end


        ezplot(fh,[minX maxX minY maxY]);
        set(h, 'Color', 'k','linewidth',2); hold on; 
        if m2==2
            ylabel('[Molecule A]')
        else
            ylabel('')
        end
        if t==3
            xlabel('[Molecule B]')
        else
            xlabel('')
        end
        title('')


        %return
        clear da; clear co1; clear co2; clear cdil1; clear cdil2; 
        clear cden1; clear cden2; clear x2; clear y2; 
        clear SSR; clear SST; clear r2; 
        clear data_points; clear theta_guaranteed; clear x; clear n; 
        clear minX; clear minY; clear maxX; clear maxY; 
        clear a1; clear b1; clear c1; clear d1; clear e1; clear f1; 
        clear fh; clear h; 
        clear myval; clear yp; 
    end
end



