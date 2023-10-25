clear all;

constructs={'Control_Rep_1','Control_Rep_2','Control_Rep_3','Mutant_Rep_1','Mutant_Rep_2','Mutant_Rep_3'};
constructs2={'416_ctrl_60x','416_ctrl','ctrl_pfa','416_bidi_60x','416_bidi','bidi_pfa'};
basedir=pwd
radius=1.3888888888888888;
pixelshift=0; %radius*[0:0.5:10];

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
    fraccoloccln3nuc=nan(length(uimagenum),length(pixelshift));
    fraccoloccln3nonnuc=nan(length(uimagenum),length(pixelshift));
    numcln3nuc=nan(length(uimagenum),length(pixelshift));
    numcln3nonnuc=nan(length(uimagenum),length(pixelshift));

    fraccolocbni1nuc=nan(length(uimagenum),length(pixelshift));
    fraccolocbni1nonnuc=nan(length(uimagenum),length(pixelshift));
    numbni1nuc=nan(length(uimagenum),length(pixelshift));
    numbni1nonnuc=nan(length(uimagenum),length(pixelshift));
    
    for i=1:length(uimagenum)
        % Load in CLN3 data
        [cln3y, cln3x]=textread(['../results/' constructs{c} '/100/CLN3-561nm_' constructs2{c} '_' uimagenum{i} '_puncta.txt'],'%d %d','headerlines',7);
        numcln3spots(i)=length(cln3x);

        % Load in BNI1 data
        [bni1y, bni1x]=textread(['../results/' constructs{c} '/100/BNI1-640nm_' constructs2{c} '_' uimagenum{i} '_puncta.txt'],'%d %d','headerlines',7);
        numbni1spots(i)=length(bni1x); 

        % Load in Nuclei data
        [nucy, nucx, nuca]=textread(['../nuclei_results/' constructs{c} '/Nuclei_' constructs2{c} '_' uimagenum{i} '_cellpose.txt'],'%f %f %f','headerlines',5);
        numnuc(i)=length(nucy);
        nucrad=sqrt(nuca/pi);

        % only consider images in which the difference in number of spots is less than an order of magnitude
        if abs(log10(numcln3spots(i)/numbni1spots(i)))<1
           
            for p=1:length(pixelshift)
                % CLN3 colocalization
                numcolocnuc=0;
                numcolocnonnuc=0;
                numnoncolocnuc=0;
                numnoncolocnonnuc=0;

                %figure;
                %for n=1:length(nucx)
                %    viscircles([nucx(n), nucy(n)],nucrad(n),'color',[0 0 0.5]); hold on; 
                %    ylim([0 2250]); xlim([0 2250]);
                %end

		for j=1:length(cln3x)
                    % Get each cln3 spots
		    currx=cln3x(j)+pixelshift(p);
		    curry=cln3y(j)+pixelshift(p);

		    % Distance between cln3 and all bni1 spots
		    mydist=sqrt((bni1x-currx).^2+(bni1y-curry).^2);
		    coloc=find(mydist<2*radius);

		    % Distance between cln3 and nucleus
                    mydist2=sqrt((nucx-currx).^2+(nucy-curry).^2);
                    nuccoloc=find(mydist2<radius+nucrad);

		    if ~isempty(coloc) & ~isempty(nuccoloc)
		        numcolocnuc=numcolocnuc+1;
                        %viscircles([currx curry],radius,'color','g'); hold on;
                        %viscircles([bni1x(coloc(1)) bni1y(coloc(1))],radius,'color','r'); hold on;
                        %viscircles([nucx(nuccoloc(1)) nucy(nuccoloc(1))],nucrad(nuccoloc(1)),'color','b'); hold on;
                    elseif isempty(coloc) & ~isempty(nuccoloc)
                        numnoncolocnuc=numnoncolocnuc+1;
                    elseif ~isempty(coloc) & isempty(nuccoloc)
                        numcolocnonnuc=numcolocnonnuc+1;
                        %viscircles([currx curry],radius,'color',[0 0.5 0]); hold on;
                        %viscircles([bni1x(coloc(1)) bni1y(coloc(1))],radius,'color',[0.5 0 0]); hold on;
                    else
                        numnoncolocnonnuc=numnoncolocnonnuc+1;
                    end

		    clear currx; clear curry; clear mydist; clear coloc; clear mydist2; clear nuccoloc; 
		 end
		 fraccoloccln3nuc(i,p)=numcolocnuc/(numcolocnuc+numnoncolocnuc);
                 fraccoloccln3nonnuc(i,p)=numcolocnonnuc/(numcolocnonnuc+numnoncolocnonnuc);
                 numcln3nuc(i,p)=numcolocnuc+numnoncolocnuc;
                 numcln3nonnuc(i,p)=numcolocnonnuc+numnoncolocnonnuc;
                 %return

	         % BNI1 colocalization
                 numcolocnuc=0;
                 numcolocnonnuc=0;
                 numnoncolocnuc=0;
                 numnoncolocnonnuc=0;

		 for j=1:length(bni1x)
                     % Get each bni1 spots
		     currx=bni1x(j)+pixelshift(p);
		     curry=bni1y(j)+pixelshift(p);

                     % Distance between bni1 and all cln3 spots
		     mydist=sqrt((cln3x-currx).^2+(cln3y-curry).^2);
		     coloc=find(mydist<2*radius);

	             % Distance between bni1 and nucleus
                     mydist2=sqrt((nucx-currx).^2+(nucy-curry).^2);
                     nuccoloc=find(mydist2<radius+nucrad);

		     if ~isempty(coloc) & ~isempty(nuccoloc)
		         numcolocnuc=numcolocnuc+1;
                         %viscircles([currx curry],radius,'color','g'); hold on;
                         %viscircles([cln3x(coloc(1)) cln3y(coloc(1))],radius,'color','r'); hold on;
                         %viscircles([nucx(nuccoloc(1)) nucy(nuccoloc(1))],nucrad(nuccoloc(1)),'color','b'); hold on;
                     elseif isempty(coloc) & ~isempty(nuccoloc)
                         numnoncolocnuc=numnoncolocnuc+1;
                     elseif ~isempty(coloc) & isempty(nuccoloc)
                         numcolocnonnuc=numcolocnonnuc+1;
                         %viscircles([currx curry],radius,'color',[0 0.5 0]); hold on;
                         %viscircles([cln3x(coloc(1)) cln3y(coloc(1))],radius,'color',[0.5 0 0]); hold on;
                     else
                         numnoncolocnonnuc=numnoncolocnonnuc+1;
                     end

		     clear currx; clear curry; clear mydist; clear coloc; clear mydist2; clear nuccoloc; 
		  end
		  fraccolocbni1nuc(i,p)=numcolocnuc/(numcolocnuc+numnoncolocnuc);
                  fraccolocbni1nonnuc(i,p)=numcolocnonnuc/(numcolocnonnuc+numnoncolocnonnuc);
                  numbni1nuc(i,p)=numcolocnuc+numnoncolocnuc;
                  numbni1nonnuc(i,p)=numcolocnonnuc+numnoncolocnonnuc;
            end

        end
        %return

        clear da; clear da2; clear cln3x; clear cln3y; clear bni1x; clear bni1y;
        clear nucy; clear nucx; clear nuca; clear nucrad; 
    end

    %numcln3nuc    
    %fraccoloccln3nuc
    %numcln3nonnuc
    %fraccoloccln3nonnuc
    numbni1nuc    
    fraccolocbni1nuc
    numbni1nonnuc
    fraccolocbni1nonnuc

    save(['../nuclei_results/' constructs{c} '/' constructs{c} '_colocalization_nuc_full_data.mat'],'pixelshift','radius','numcln3spots','numbni1spots','uimagenum','numcln3nuc','fraccoloccln3nuc','numcln3nonnuc','fraccoloccln3nonnuc','numbni1nuc','fraccolocbni1nuc','numbni1nonnuc','fraccolocbni1nonnuc');

    %return
    clear currdir; clear fileList; clear imagenum; clear uimagenum; clear numcln3spots; clear numbni1spots; 
    clear fraccoloccln3; clear fraccolocbni1; clear fccln3; clear fcbni1; clear numcoloccln3; clear numcolocbni1; 
    clear fraccolocbni1nuc; clear fraccolocbni1nonnuc; clear numbni1nuc; clear numbni1nonnuc; 
    clear fraccoloccln3nuc; clear fraccoloccln3nonnuc; clear numcln3nuc; clear numcln3nonnuc; 
end 
        
