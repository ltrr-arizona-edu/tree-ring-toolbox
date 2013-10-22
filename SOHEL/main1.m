%
% USAGE : MAIN1
%   A script .m file for time series analysis of Tree Rings
%
% S. Anwar     Last Revised :  8/17/94
%
%
% USER WRITTEN FUNCTIONS NEEDED
%------------------------------
% CLSGRF.M	A script .m file called by uicontrol which closes
%		all figure windows
% HYGEOM1.M	Hypergeometric test function
%  -HYGHLP.M		Computes different statistics
%  -QUANTILE.M		Picks the data values outside the quantile
%  -SLVMNU.M		A modified menu function
% JDISP.M	A non-interactive window display function
% PEARSPE2.M	Tests the data using Pearson or Spearman Methods
%  -PEARSP.TAB		Table for critical values
% PROCES1.M	A postprocessing function
%  -INTPRS1.M		An intermediate processing function
%  -INTPRS2.M		An intermediate processing function
%    -PLTEXT1.M			Places text inside figures
%  -INTPRS3.M		An intermediate processing function
%    -PLTEXT1.M			Places text inside figures
%  -PLTEXT1.M		Places text inside figures
% RESCALY.M	A scaling function
% RWCHNG.M	Transforms the data
% RWCULL.M	Picks and plots segment of X series
%  -JDISP.M		A non-interactive window display function
%  -PLTEXT.M		Places text inside figures
%  -SHLMNU.M		A modified menu function
% RWINP.M	Reads the RW file
%  -JDISP.M		A non-interactive window display function
%  -UDISP.M		A user interactive window display function
% RWPICK.M	Picks X and Y series
%  -SVMENU.M		A modified menu function
% SHLMNU.M	A modified menu function
% SIGNT3.M	Sign test function
%  -SGNTBL2.TAB		Table of critical values
% SLHMNU.M	A modified menu function
% SLVMNU.M	A modified menu function
% SVTMNU.M	A modified menu function
% UDISP.M	A user interactive window display function
% USINP.M	A user interactive function to take 'YES' or 'NO'
% XYCOLW	Duplicates X-segment and Y in matrices
%  -UDISP.M		A user interactive window display function
%___________________________________________________________________

% Prompt for the number of files to be read
f1=slhmnu('# of data files to be read ?','1','2','3');

for i=1:f1,
  disp(['Reading the data file name for Z',num2str(i)]);
  [zname,zy,z,zind,fl]=rwinp;		% Read the files
  if fl==-1,   % Failure in reading a file breaks the loop
    break;
  end
  eval(['zname',num2str(i),'=zname;']);
  eval(['zy',num2str(i),'=zy;']);
  eval(['z',num2str(i),'=z;']);
  eval(['zind',num2str(i),'=zind;']);
  eval(['ncol',num2str(i),'=length(zind)-1;']);
end

% If a file cannot be read, return to the command window
if fl==-1,
  jl=jdisp(['Error in Reading the file. Please rerun the Program. ',...
     'Press any key to continue.']);
  pause;
  close(jl);
  return;
end

k1=1;
while k1~=2,  % if k1=2, quit the program
 k1=shlmnu('Select One','X and Y series','Quit');
 if k1==1,
  for i=1:2,
    m1=1;
    % Pick the file Z for X or Y series
    if i==1,
      var='X';
    else
      var='Y';
    end
    disp(['Pick the vector Z for ',var,' series']);
    if f1==2,
      m1=shlmnu([var,' series ?'],'Z1','Z2');
    elseif f1==3,
      m1=shlmnu([var,' series ?'],'Z1','Z2','Z3');
    end

    ms1=num2str(m1);
    zname=eval(['zname',ms1]);
    zy=eval(['zy',ms1]);
    z=eval(['z',ms1]);
    zind=eval(['zind',ms1]);
    ncol=eval(['ncol',ms1]);
    disp(['Choose your ',var,'-series now ']);
    % Pick X and Y series
    [sz,yr,nm]=rwpick(zname,z,zy,zind,var);
    eval([var,'=sz;']);
    eval(['y',var,'=yr;']);
    eval(['nm',var,'=nm(1:6);']);
  end		% End of for loop i=1:2

  xyl=['RW   vs  RW '
       'RW   vs  IND'
       'IND  vs  IND'
       'SKE  vs  SKE'
       'SKE  vs  RW '
       'SKE  vs  IND'];
  rscl=svtmnu('Choose Mode','X          Y',xyl);

  % Set the values of ks for Hypergeometric test
  if rscl<=3,
    ks=1;
  elseif rscl==4,
    ks=2;
  else
    ks=3;
  end

  % Set orgn 'original Series'
  orgn=' ';  

  % Ask if transformation should be done
  ktr=usinp('Would you like to transform the X and Y series ?');	

  if ktr,

    if rscl<=3,
      X=rwchng(X);
      Y=rwchng(Y);
      yX(1)=[];
      yY(1)=[];
      orgn=''' ';
    elseif (rscl==5 | rscl==6),
      Y=rwchng(Y);
      X(1)=[];
      yX(1)=[];
      yY(1)=[];
      orgn=''' ';
    end  
    
    % Rescaling Y time series in mixed mode
    if (rscl==2 | rscl==5 | rscl==6),
      Y=rescaly(X,Y);			
    end

  end

  k2=1;
  while k2<2,  % If k2=3, quit. 
    k2=shlmnu('Select One','X segment','New X and Y series','Quit');
    if k2==1,
     % Pick and plot X segment
     [xseg,tseg] = rwcull(X,Y,yX,yY,nmX,nmY,orgn,rscl);
     lxs=length(xseg);
     [xsegm,YM,yrY,fl] = xycolw(xseg,Y,yY);
     if fl==-1, % Flag indicating wrong xseg length
       return;
     end
     tsm=1;
     while tsm<4,
       % Ask for Spearman/Pearson or Sign Test
       tsm=shlmnu('Specify Test','Pearson/Spearman','Sign',...
                  'Hypergeometric','Quit');
        if tsm==1,
          kps=shlmnu('Choose Test','Pearson','Spearman');
	  [r,npsp] = pearspe2(xsegm,YM,kps);
	  dfc=xyl(rscl,:);
	  kps=[kps npsp];
	  [lnt,lny,h1] = proces1(r,YM,yrY,nmX,nmY,tseg,xseg,orgn,rscl,tsm,kps,dfc);
          close(h1);
       elseif tsm==2,
	  kst=shlmnu('Choose','Sign Departures','Sign First-Diff');
	  [na,ns,nkey]=signt3(xsegm,YM,kst);
	  dfc=xyl(rscl,:);
	  [nkd,nks]=size(nkey);
 	  kst=[kst,max(ns),min(ns),zeros(1,nks-3);nkey];
	  [lnt,lny,h1] = proces1(na,YM,yrY,nmX,nmY,tseg,xseg,orgn,rscl,tsm,kst,dfc);
          close(h1);
       elseif tsm==3,
	 kq=shlmnu('Type of Test ?','Lower','Upper','Two-sided');
	 klq=usinp('Plot Confidence Levels ?');
	 if kq==1,
 	   skq='Lower-Sided ';
	   [N,m,ps,n95,n99,p]=hygeom1(xsegm,YM,kq,ks,klq);
	 elseif kq==2,
	   skq='Upper-Sided ';
	   [N,m,ps,n95,n99,p]=hygeom1(xsegm,YM,kq,ks,klq);
	 elseif kq==3,
	   skq=['Lower-Sided ';'Upper-Sided '];
	   [N,m,ps,n95,n99,p]=hygeom1(xsegm,YM,1,ks,klq);
	   [N1,m1,ps1,n951,n991,p1]=hygeom1(xsegm,YM,2,ks,klq);
     	 end
	 dfc=xyl(rscl,:);
	 dfc=[dfc;skq];
	 if kq<3,
 	   format short e;
	   pb=1-p;
	   format;
	   np=[N' ps' n95' n99'];
	   [lnt,lny,h1] = proces1(m,YM,yrY,nmX,nmY,tseg,xseg,orgn,rscl,tsm,np,dfc,pb);	
	 elseif kq==3,
	   m=[m' m1'];
	   format short e;
	   p=[1-p;1-p1];
	   format;
	   np=[N' ps' n95' n99' N1' ps1' n951' n991'];
	   [lnt,lny,h1] = proces1(m,YM,yrY,nmX,nmY,tseg,xseg,orgn,rscl,tsm,np,dfc,p);	
	 end  
         close(h1);
       elseif tsm==4;
       end			% End of if Loop (tsm)
     end			% End of while Loop (tsm < 4)
    elseif k2==2;
    elseif k2==3,
     k1=2;
    end           % End of k2 if loop
  end             % End of second while loop
  elseif k1==2;
  end              % End of k1 if loop
end                % End of first while loop

% Place a pushbutton in the current figure window asking if the 
% user wanted to close the graphs.
%figure(1);
%psh=uicontrol(gcf,'Style','Pushbutton','String','Close graphs ?',...
%  'Position',[400 350 120 30],'Callback','clsgrf');

close all; % Close all figure windows

clear all; % Erase all data from memory. 

% End of file
