% Copyright (C) Xiangyi Meng.

function [lfn, rfn] = getFaceNumber_cuboid(fk)
    switch fk
        case 1
            lfn = 2; rfn = 1;
        case 2
            lfn = 4; rfn = 3;
        case 3
            lfn = 6; rfn = 5;
        case 4
            lfn = 11; rfn = 1;
        case 5
            lfn = 12; rfn = 1;
        case 6
            lfn = 13; rfn = 1;
        case 7
            lfn = 14; rfn = 1;
        case 8
            lfn = 2; rfn = 7;
        case 9
            lfn = 2; rfn = 8;
        case 10
            lfn = 2; rfn = 9;
        case 11
            lfn = 2; rfn = 10;
        case 12
            lfn = 19; rfn = 3;
        case 13
            lfn = 20; rfn = 3;
        case 14
            lfn = 21; rfn = 3;
        case 15
            lfn = 22; rfn = 3;
        case 16
            lfn = 4; rfn = 15;
        case 17
            lfn = 4; rfn = 16;
        case 18
            lfn = 4; rfn = 17;
        case 19
            lfn = 4; rfn = 18;
        case 20
            lfn = 27; rfn = 5;
        case 21
            lfn = 28; rfn = 5;
        case 22
            lfn = 29; rfn = 5;
        case 23
            lfn = 30; rfn = 5;
        case 24
            lfn = 6; rfn = 23;
        case 25
            lfn = 6; rfn = 24;
        case 26
            lfn = 6; rfn = 25;
        case 27
            lfn = 6; rfn = 26;
        otherwise
            error('Wrong face kind')
    end
end

