function cfl = setCFLNumber(k)

switch k  
    case 0
        cfl = 0.98;
    case 1      
        cfl = 0.26;
    case 2        
        cfl = 0.18;  
    case 3
        cfl = 0.14;
    case 4
        cfl = 0.1; 
    otherwise
        error('Wrong polynomial degree')
end 

end