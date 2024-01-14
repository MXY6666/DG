% Copyright (C) Xiangyi Meng.

function [lfn, rfn] = getFaceNumber_rect(fk)
    switch fk
        case 1
            lfn = 2; rfn = 1;
        case 2
            lfn = 4; rfn = 3;
        case 3
            lfn = 7; rfn = 1;
        case 4
            lfn = 8; rfn = 1;
        case 5
            lfn = 2; rfn = 5;
        case 6
            lfn = 2; rfn = 6;
        case 7
            lfn = 11; rfn = 3;
        case 8
            lfn = 12; rfn = 3;
        case 9
            lfn = 4; rfn = 9;
        case 10
            lfn = 4; rfn = 10;
        otherwise
            error('Wrong face kind')
    end
end