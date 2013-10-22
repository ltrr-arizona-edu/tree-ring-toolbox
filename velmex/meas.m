clear t N;
x=repmat(NaN,1000,1);
obj=serial('com1');
fopen(obj);
set(obj,'Parity','none','StopBits',2,'DataBits',8,'BaudRate',9600');
set(obj,'terminator','CR/LF','TimeOut',10);
kwh1=1; % while control
n=0;
while kwh1==1;
    nval=get(obj,'ValuesReceived');
    if nval~=0;
        flushinput(obj);
    end;
    [tline,count] = fgets(obj);
    g=tline;
    i1=findstr(g,'mm');
    if~isempty(i1) & length(i1)==1;
        isp=i1-1;
        igo=isp-7;
        d = g(igo:isp);
        dnum = str2num(d);
        if dnum<=0;
            instrreset;
            error('reading must be greater than zero');
        end;
        dstr=num2str(dnum,'%8.3f');
    else;
        instrreset;
        error(['No  mm string, or more than one, in tline']);
        
    end;
        
            
    n=n+1;
    x(n)=dnum;
    if n>3;
        kwh1=0;
    end;
end; % while
x(isnan(x))=[];
instrreset;