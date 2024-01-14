function cfl = setCFLNumber(k)

switch k  
    case 1
        cfl = 1;
    case 2
        cfl = 1;
    case 3
        cfl = 1;
    case 4
        cfl = 1;
    otherwise
        error('Please input 1, 2， 3 or 4 as input polynomial degree')
end 

end

