function [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu)
% menucell:  menu for selecting elements of a cell array
% [C3,Lout] = menucell (C1,C2,Lin,,nallow,strmenu);
% Last revised 2-11-02
%
% Menu for selecting elements of a cell array.  Originally written to be called by 
% rwlset.m to specify which series should be included in a .rwl file of ring widths
%
%*** INPUT
%
% C1{} cell array of names, which are strings
% C2{} cell array telling whether series are already selected ('-Y') or not ('-N')
% Lin(? x 1)L  logical vector, same size as C1, with 1 for selected, 0 for not, on input
% nallow (1 x 2)i  minimum and maximum allowable numbers of selected series on output
% strmenu (1 x ?)s   Question starting menu (e.g., 'Toggle to Y to accept series')
%
%
%*** OUTPUT
%
% C3{} cell array of series names of selected series
% Lout (? x1)L indicator whether on output series is selected or not (1 or 0) on output
%
%*** REFERENCED -- NONE
%
%*** UW FUNCTIONS CALLED --NONE
%*** TOOLBOXES NEEDED -- NONE
%*** NOTES


nser = length(C1);
D1=C1;
D1{nser+1}='Accept';
D2=C2;
D2{nser+1}='';
D0=cellstr([char(D1) char(D2)]);
D=D0;
Lcurr=Lin;

mess1=['Number of selected series must be in range ' int2str(nallow)];

kwh1 = 1;
while kwh1==1;
    kmen1=menu(strmenu,D);
    if kmen1==nser+1;
        npicked = sum(Lcurr);
        if npicked<nallow(1) | npicked>nallow(2);
            uiwait(msgbox(mess1,'Message','modal'));
        else;
            Lout=Lcurr;
            C3=C1(Lout);
            kwh1=0;
        end;
    else; % You have toggled a series
        if Lcurr(kmen1)==1; % was on, now will be off
            Lcurr(kmen1)=0;
            D2{kmen1}='-N';
        elseif Lcurr(kmen1)==0; 
            Lcurr(kmen1)=1;
            D2{kmen1}='-Y';
        end;
        % Revise menu
        D=cellstr([char(D1) char(D2)]); 
    end;
   
end; 
