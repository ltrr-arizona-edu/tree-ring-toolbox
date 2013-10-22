function yr = rwgoyr(N,lstyr)
% rwgoyr:  start year given rw file with last valid data in row N corresponding to lstyr
% CALL: yr = rwgoyr(N,lstyr);
yr = lstyr -N + 6;