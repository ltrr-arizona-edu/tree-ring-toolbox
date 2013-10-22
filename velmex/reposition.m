function reposition(obj, reclen)
%reposition Prompt the user to reposition the stage, read the first value.
%   Given a serial object, obj, and a minimum record length, reclen (in
%   characters), this prompts the user to reposition the stage and waits
%   for the initial value (which should be zero).
  
initdlg = msgbox({'Initialize segment as follows:', ...
        '  * Position the stage', ...
        '  * press RESET -- sets Quick-Check display to zero', ...
        '  * press PRINT -- marks starting point for measuring'}, ...
        'Message');
    invalid = true;
    while invalid
        try
            count = 0;
            while count < reclen
                [tline, count] = fgets(obj);
            end
            x0 = decodeline(tline);
            invalid = (x0 ~= 0);
            if invalid
                disp('The initial value must be zero, please try again.');
            end
        catch DLEX
            disp(DLEX.message);
            trid = findstr(DLEX.identifier, 'TreeRing');
            if isempty(trid)
                closeserial(sport);
                rethrow(DLEX);
            end
        end
    end
    if ishandle(initdlg)
        close(initdlg);
        pause(1);
        beep();
        pause(1);
        beep();
    end
    disp(['Read the initial value of ', tline]);
    disp('OK, crank and click to record measurements');
end

