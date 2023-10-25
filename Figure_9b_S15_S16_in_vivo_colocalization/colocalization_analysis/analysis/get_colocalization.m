clear all;

constructs={'Control_Rep_1','Control_Rep_2','Control_Rep_3','Mutant_Rep_1','Mutant_Rep_2','Mutant_Rep_3'};
constructs2={'416_ctrl_60x','416_ctrl','ctrl_pfa','416_bidi_60x','416_bidi','bidi_pfa'};
basedir=pwd
radius=1.3888888888888888;
pixelshift=radius*[0:0.5:10];

for c=1:length(constructs)
    currdir=['../results/' constructs{c} '/100'];
    cd(['../results/' constructs{c} '/100'])
    fileList = dir('*puncta.txt');

    for f=1:length(fileList)
        tmp=fileList(f).name;
        tmp2=strsplit(tmp,'_');
        if c==1 | c==4
            imagenum{f}=tmp2{5};
        else
            imagenum{f}=tmp2{4};
        end
        clear tmp; clear tmp2; 
    end

    uimagenum=unique(imagenum);
    cd([basedir]);

    % Create matrices for colocalization 
    fraccoloccln3=nan(length(uimagenum),length(pixelshift));
    fraccolocbni1=nan(length(uimagenum),length(pixelshift));
    
    for i=1:length(uimagenum)
        % Load in CLN3 data
        [cln3y, cln3x]=textread(['../results/' constructs{c} '/100/CLN3-561nm_' constructs2{c} '_' uimagenum{i} '_puncta.txt'],'%d %d','headerlines',7);
        numcln3spots(i)=length(cln3x);

        % Load in BNI1 data
        [bni1y, bni1x]=textread(['../results/' constructs{c} '/100/BNI1-640nm_' constructs2{c} '_' uimagenum{i} '_puncta.txt'],'%d %d','headerlines',7);
        numbni1spots(i)=length(bni1x);   

        % only consider images in which the difference in number of spots is less than an order of magnitude
        if abs(log10(numcln3spots(i)/numbni1spots(i)))<1
           
            for p=1:length(pixelshift)
                % CLN3 colocalization
                numcoloc=0;
		for j=1:length(cln3x)
		    currx=cln3x(j)+pixelshift(p);
		    curry=cln3y(j)+pixelshift(p);
		    mydist=sqrt((bni1x-currx).^2+(bni1y-curry).^2);
		    coloc=find(mydist<2*radius);
		    if ~isempty(coloc)
		       numcoloc=numcoloc+1;
                       %xy1(numcoloc,1)=currx;
                       %xy1(numcoloc,2)=curry;
                       %xy2(numcoloc,1)=bni1x(coloc(1));
                       %xy2(numcoloc,2)=bni1y(coloc(1));
                       %if p==1
                       %   figure(290+c);
                       %   %subplot(5,5,i)
                       %   viscircles([currx curry],radius,'color','b'); hold on;
                       %   viscircles([bni1x(coloc(1)) bni1y(coloc(1))],radius,'color','r'); hold on;
                       %   ylim([0 2250]); xlim([0 2250]);
		       %end
		       %return
		    end
		    clear currx; clear curry; clear mydist; clear coloc;
		 end
		 fraccoloccln3(i,p)=numcoloc/length(cln3x);
                 numcoloccln3(i,p)=numcoloc;
                 %return

                 % BNI1 colocalization
                 numcoloc=0;
		 for j=1:length(bni1x)
		     currx=bni1x(j)+pixelshift(p);
		     curry=bni1y(j)+pixelshift(p);
		     mydist=sqrt((cln3x-currx).^2+(cln3y-curry).^2);
		     coloc=find(mydist<2*radius);
		     if ~isempty(coloc)
		        numcoloc=numcoloc+1;
		        %return
		     end
		     clear currx; clear curry; clear mydist; clear coloc; 
		  end
                  fraccolocbni1(i,p)=numcoloc/length(bni1x);
                  numcolocbni1(i,p)=numcoloc;
            end

        end

        clear da; clear da2; clear cln3x; clear cln3y; clear bni1x; clear bni1y;
    end

    figure;
    subplot(1,2,1)
    boxplot(fraccoloccln3);

    subplot(1,2,2)
    boxplot(fraccolocbni1);

    return

    % Save cln3 colocalization data
    fccln3(1,:)=pixelshift;
    fccln3(2:size(fraccoloccln3,1)+1,:)=fraccoloccln3;
    csvwrite(['../results/' constructs{c} '/100/' constructs{c} '_CLN3_colocalization.csv'],fccln3);

    fcbni1(1,:)=pixelshift;
    fcbni1(2:size(fraccolocbni1,1)+1,:)=fraccolocbni1;
    csvwrite(['../results/' constructs{c} '/100/' constructs{c} '_BNI1_colocalization.csv'],fcbni1);

    save(['../results/' constructs{c} '/100/' constructs{c} '_colocalization_full_data.mat'],'fraccoloccln3','fraccolocbni1','pixelshift','radius','numcln3spots','numbni1spots','uimagenum','numcoloccln3','numcolocbni1');

    %return
    clear currdir; clear fileList; clear imagenum; clear uimagenum; clear numcln3spots; clear numbni1spots; 
    clear fraccoloccln3; clear fraccolocbni1; clear fccln3; clear fcbni1; clear numcoloccln3; clear numcolocbni1; 
end
