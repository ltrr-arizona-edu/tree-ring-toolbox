function X=subfcn03(X,k)
% update the summary

n= length(X.id); % number of series measured

if k==1; % total rw
    str1='Total Ring Width';
    str2=' ID         Total RW';
elseif k==2; % EWW/LWW
    str1='EWW/LWW';
    str2=['ID            EWW           LWW           EWW+LWW'];
end;

s1=[int2str(n) ' series  have been measured for ' str1];
s1=char(s1,'Spans of measurements are as follows');
s1=char(s1,str2);


for j = 1:n;
    spann=X.span{j};
    if k==1;;
        s2=[X.id{j} '      '  int2str(spann(1,:))];
    elseif k==2;
        s2 = [X.id{j} '     '   int2str(spann(1,:)) '     ' int2str(spann(2,:)) '     '  int2str(spann(3,:))];
    end;
    s1=char(s1,s2);
end;

X.summary=s1;
