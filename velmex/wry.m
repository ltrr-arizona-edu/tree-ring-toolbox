% wry.m

% Hard Code
clear;
close all;
x=repmat(NaN,800,1);
clc;

yrstart=1601;

% Create serial port object
instrreset; % Just in case, make sure serial port not assigned to some object
obj=serial('com1');

% Set input buffer so that read completes when one measurement is received
set(obj,'InputBufferSize',16);

% Open serial port and set parameters
fopen(obj);
set(obj,'Parity','none','StopBits',2,'DataBits',8,'BaudRate',9600');
set(obj,'terminator','','TimeOut',10);

% Initialize counters and while control
n=0;
nx=0;
count=0;
kwh1=1;
irev=0; % flag for crank being reversed

% Measurement while loop
while kwh1==1;
    while count==0; % Ready to get a measurement; either this is first measurement, or input has been flushed
        % get data from serial port
        if n==0;
            uiwait(msgbox({'Initialize segment as follows:','  * Positioning the stage',...
                    '  * press RESET -- sets Quick-Check display to zero',...
                    '  * press PRINT -- marks starting point for measuring'},'Message','modal'));
        end;
        [g,count] = fgets(obj); % store measurement (a string) in g 
    end;
    if count==16;  % One measured value is in buffer
        flushinput(obj);   % flush the current measurement from the buffer
        i1=findstr(g,'mm');  % mm string comes after measured value from QUICKCHECK
        
        if~isempty(i1); % "mm" string found; we have data from QuickCheck
            % Cull the measurement from g and convert the measurement from sting to numeric
            isp=i1-1;
            igo=isp-7;
            d = g(igo:isp);
            dnum = str2num(d);
            if n~=0;
                if dnum<0;
                    if irev==1; 
                        irev=0;
                    else;
                        instrreset;
                        error('Quick-Check Display is negative-- invalid except before first ring measured');
                    end;
                end;
            else; % first click, must reset velmex before this
                if dnum~=0;
                    instrreset;
                    error('Must reset Quick-Check to 0 and Press PRINT to start measuring any segment');
                else;
                    disp('OK, ZERO received as baseline from QuickCheck');
                end;
            end;
            dstr=num2str(dnum,'%8.3f');
        else;
            instrreset;
            error(['No  mm string']);
        end;
        
        n=n+1; % Increment index for storing readings from QC; if first reading, n==1, etch
        nx=nx+1; % index into storage vector x
        n1=nx-1;
        
        yrthis = yrstart+n1-1;
        c(n)=dnum; % cumulative counter
        nc=length(c);
        if n>1;
            w=c(nc)-c(nc-1); % measurement is difference of cumulative reading
           
        end;
        
        if n>1 &  w<0; % Check for negative difference -- meaning you reversed crank on velmex
            kmen4=menu('Choose','Stop measuring',...
                'Insert missing value(s)',...
                'Remeasure some part');
            switch kmen4;
            case 1;
                disp('OK, Measurement stopped');
                n=n-1;
                nx=nx-1;
                yr = yrstart+(1:nx)-1;
                kwh1=0;
            case 2; %  insert missing value(s)
                irev=1; % you have reversed the crank
                yrthis=yrthis-1;
                n=n-1;
                nx=nx-1;
                n1=nx-1;
                prompt={'Enter the Year to resume measurements (after the missing values):'};
                def={int2str(yrthis+2)};
                dlgTitle='Input year';
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def);
                yrresume = str2num(answer{1});
                nfill = yrresume-yrthis-1;
                x((n1+1):(n1+1+nfill-1))=NaN;
                nx=nx+nfill-1;
                n1=nx-1;
                %yr = (yrstart:(yrresume-1))';
                % HERE RESUME EDITING
                n=0;
                clear c;
            case 3; % remeasure some part
                irev=1; % you have reversed the crank
                % Jumps to a specified year to begin remeasuring
                % strp off negative trailing measurement
                yrthis=yrthis-1; % adjust current year for stripping
                n=n-1;  % adjust current index
                nx=nx-1; % adjust index into stored series
                n1=nx-1;
                prompt={'Enter start year of segment to begin remeasuring):'};
                def={int2str(yrthis+2)};
                dlgTitle='Input year';
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def);
                yrresume = str2num(answer{1});
                nx = yrresume-yrstart;
                n1=nx-1;
                %yr = (yrstart:(yrresume-1))';
                % HERE RESUME EDITING
                n=0;
                clear c;
                
            otherwise;
            end;
            
        else;
            if n>1;
                x(n1)=w;
                yrthis = yrstart+n1-1;
                nn=length(x);
                str1=num2str([n1 yrthis x(n1)]); 
                disp(str1);
                
                
            end;
        end;
        
        count=0;
    end;
end; % while
x=trailnan(x); % stripp off trailing NaNs
nx = length(x);
yr = (1:nx)' - 1 +yrstart;


instrreset;