% cdnow7 -- change to current directory of choice

pathfull7;
kmen1 = menu('Choose Current Directory',...
    'c:\work1\',...
    'c:\projs\ai3\',...
    'c:\angelika\SALTRIVER\BRF\',...
    'c:\angelika\',...
    'c:\kiyomi\',...
    'c:\occ\',...
    'c:\Amy\');
switch kmen1;
    case 1;
        cd c:\work1\;
    case 2;
        cd c:\projs\ai3\;
    case 3;
        cd c:\angelika\SALTRIVER\BRF\;
    case 4;
        cd c:\angelika\;
    case 5;
        cd c:\kiyomi\;
    case 6;
        cd c:\angelika\occ\;
    case 7;
        cd c:\Amy\;
end;
