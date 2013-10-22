function trivialserial()
%trivialserial tests the display of a message box while reading the serial
%port.
scheck = instrfind;
if ~isempty(scheck)
    fclose(scheck);
    delete(scheck);
    clear scheck;
end
sport = serial('com1');
set(sport, ...
    'Parity', 'none', ...
    'StopBits', 1, ...
    'DataBits', 8, ...
    'BaudRate', 9600', ...
    'terminator', 'CR/LF', ...
    'TimeOut', 600);
fopen(sport);
reclen = 15;
count = 0;
initdlg = msgbox('Please send something to the serial port.','Message');
while count < reclen
    [tline, count] = fgets(sport);
end
if ishandle(initdlg)
    close(initdlg);
end
disp(tline);
clear initdlg;
while count >= reclen && isempty(findstr(tline, '-'))
    [tline, count] = fgets(sport);
    disp(tline);
end
disp('Finished.');
fclose(sport);
delete(sport);
clear sport;
end

