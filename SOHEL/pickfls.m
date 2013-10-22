function pickfls
%
% A function to pick several data files
%
% Writes the file info in a .cat file
%__________________________________________________________

cd c:\mlb\data;

kfl=1;
znam=[];
while kfl==1,
  kfl=shlmnu('Read a file ?','YES','NO/QUIT');
  if kfl==1,
    fln=uigetfile('*.dat','File name ?');

    % If Cancel button is pressed in the input window, flag is 
    % set to -1 and control is returned to calling function
    if fln==0,
      fl=-1;
      return;
    end 

    % open the data file for reading
    fid=fopen(fln,'r');

    % Read the total # of columns in the file from first row
    nstr=fgetl(fid);
    n=str2num(nstr);
    if length(nstr)>2 | isstr(n),
      fl=-1;
      disp('Wrong file type : Need the correct data file name');
      return;
    end
    disp(['The total Number of columns read is : ' nstr]);

    % Read the first line of the file ; this line contains the 
    % name of the samples
    zn=fgetl(fid);

    % Read rest of the data file in %g format and store the data 
    % in a vector
    zv=fscanf(fid,'%g',inf);
    tm=length(zv);
    m=tm/n;         % Number of rows in the data matrix

    % Form an m x n matrix including the years 
    zt=reshape(zv,n,m);
    zt=zt';

    % Check if the # of column read is right
    B=zt(:,1);
    C=ones(length(B)-1,1);
    LC=C==diff(B);
    if LC==1,
      fl=1;
      disp('Right number of columns read');
    else
      fl=-1;
      disp('The number of columns read is wrong : Execution terminated');
      return;
    end

    % Initialize the return variables
    zyr=zeros(m,1);
    zyrs=zeros(n-1,2);

    zyr=zt(:,1);
    L1=zt~=-99;
    L1(:,1)=[];

    for i=1:n-1,
      % Put the beginning and ending years in a matrix
      zyrs(i,:)=[min(zyr(L1(:,i))) max(zyr(L1(:,i)))];
    end

    jhdl=jdisp(['Successful Reading of the data file ' fln]);
    pause(2);
    close(jhdl);	% close the jdisp window

    flnth=length(fln);
    if flnth<=16,
      dum=[];
      for i=1:20-flnth,
        dum=[dum,' '];
      end
    else
      flg=-1;
    end

    if flg==-1,
      rjh=jdisp('The file name is too long. Program terminated');
      pause(2);
      break;
    end

    n1st=num2str(n-1);
    dm1=[];
    for i=1:10-length(n1st),
      dm1=[dm1,' '];
    end

    fnam=[fln,dum,n1st,dm1,num2str(min(zyrs(:,1))),'     ',num2str(max(zyrs(:,2)))];

    fclose(fid);  	% close the input file

    znam=[znam;fnam];

  end			% End of outer if loop
  
end  			% End of while loop

if flg==-1,
  return;
end

[na,nb]=size(znam);

flw=uiputfile('*.cat','Enter a Name');
fid=fopen(flw,'w');

for i=1:na,
  for j=1:nb,
    fprintf(fid,'%s',znam(i,j));
  end
  fprintf(fid,'\n');
end

fclose(fid); 		% Close the write file

% End of file