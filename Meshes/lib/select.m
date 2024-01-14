% Prompts user to select one of a set of standard domains and sets certain 
% constants based on that selection
function [idmn, a1, b1, a2, b2, a3, b3, a4, b4] = select

disp(' 1) Square')
disp(' 2) Rectangle')
disp(' 3) Parallelogram')
disp(' 4) Trapezoid')
disp(' 5) Quadrilateral')
disp(' 6) Annulus/Horseshoe')
disp(' 7) Swan')
disp(' 8) Dome/Valley')
disp(' 9) Airfoil')
disp('10) Chevron')
disp('11) Horn')
disp('12) Backstep')
disp('13) Plow')
disp('14) C')
disp('15) New Airfoil')
disp('16) Tie')
disp('17) Shell')
disp('18) New Horseshoe')
disp('19) Ring')
disp('20) NACA0012_O-Type')
disp('21) NACA0012_C-Type')
disp('22) Slope')

idmn = input('Input one of the above geometry number : ');
idef = input('1) Default Domain, or 2) Tailor? ');
a1 = []; a2 = []; a3 = []; a4 = [];
b1 = []; b2 = []; b3 = []; b4 = [];

% Unit Square
if (idmn == 1) 
    a1 = 0; b1 = 0;
    if (idef == 1) 
        a3 = 1;
    else
        a3 = input('Input a3 : ');
    end 
    a2 = a3;
    b2 = 0;
    b3 = a3;
    a4 = 0;
    b4 = a3;
end 
% Rectangle
if (idmn == 2) 
    a1 = 0; b1 = 0;
    if (idef == 1) 
        a3 = 2;
        b3 = 1;
    else
        a3 = input('Input a3 : ');
        b3 = input('Input b3 : ');
    end 
    a2 = a3;
    b2 = 0;
    a4 = 0;
    b4 = b3;
end 
% Parallelogram
if (idmn == 3) 
    a1 = 0; b1 = 0;
    if (idef == 1) 
        a2 = 2;
        a3 = 2.5;
        b3 = 1;
    else
        a2 = input('Input a2 : ');
        a3 = input('Input a3 : ');
        b3 = input('Input b3 : ');
    end 
    b2 = 0;
    a4 = a3 - a2;
    b4 = b3;
end 
% Trapezoid
if (idmn == 4) 
    if (idef == 1) 
        a2 = 0.5;
        a3 = 1;
        b3 = 1;
    else
        a2 = input('Input a2 : ');
        a3 = input('Input a3 : ');
        b3 = input('Input b3 : ');
    end 
    a1 = -a2;
    b1 = 0;
    b2 = 0;
    a4 = -a3;
    b4 = b3;
end
% Quadrilateral
if (idmn == 5) 
    a1 = 0; b1 = 0;
    if (idef == 1) 
        a2 = 2;
        b2 = 0;
        a3 = 1.75;
        b3 = 1.5;
        a4 = 0.5;
        b4 = 1;
    else
        a2 = input('Input a2 : ');
        b2 = input('Input b2 : ');
        a3 = input('Input a3 : ');
        b3 = input('Input b3 : ');
        a4 = input('Input a4 : ');
        b4 = input('Input b4 : ');        
    end
end
% Annulus/Horseshoe
if (idmn == 6) 
    if (idef == 1) 
        a1 = 1;
        r = 2;
        ar = 1;
    else
        a1 = input('Input a1 : ');
        r = input('Input r : ');
        ar = input('Input ar : ');
    end
    a2 = ar * a1;
    a3 = r * a1;
    a4 = ar * a3;
end
 % Swan
if (idmn == 7) 
    if (idef == 1) 
        a1 = 0.5;
        a2 = 0.5;
    else
        a1 = input('Input a1 : ');
        a2 = input('Input a2 : ');
    end 
end 
% Dome/Valley
if (idmn == 8) 
    if (idef == 1) 
        a1 = 0.5;
    else
        a1 = input('Input a1 : ');
    end 
end 
% Airfoil
if (idmn == 9) 
    if (idef == 1) 
        a1 = 0.25;
    else
        a1 = input('Input a1 : ');
    end
end
% Chevron
if (idmn == 10) 
    if (idef == 1) 
        a1 = -0.5;
        a2 = -0.5;
    else
        a1 = input('Input a1 : ');
        a2 = input('Input a2 : ');
    end
end 
% Horn
if (idmn == 11) 
    if (idef == 1) 
        a1 = 4;
    else
        a1 = input('Input a1 : ');
    end
end 
% Backstep
if (idmn == 12) 
    if (idef == 1) 
    a1 = 2;
    b1 = 1;
    else
        a1 = input('Input a1 : ');
        b1 = input('Input b1 : ');        
    end
end
% Plow
if (idmn == 13) 
    if (idef == 1) 
        a1 = 1.5;
        b1 = pi / (2 + pi);
    else
        a1 = input('Input a1 : ');
        b1 = input('Input b1 : ');   
    end
end
% C-grid
if (idmn == 14) 
    if (idef == 1) 
        a1 = 0.5;
        b1 = 0.05;
        a2 = 1;
        b2 = 0.5;
    else
        a1 = input('Input a1 : ');
        b1 = input('Input b1 : ');  
        a2 = input('Input a2 : ');
        b2 = input('Input b2 : ');          
    end 
end 
% New Airfoil
if (idmn == 15) 
    if (idef == 1) 
    a1 = 0.5;
    a2 = 2;
    b2 = 1;
    a3 = 1.5;
    b3 = 0.0625;
    else
        a1 = input('Input a1 : ');
        a2 = input('Input a2 : ');
        b2 = input('Input b2 : ');   
        a3 = input('Input a3 : ');
        b3 = input('Input h : ');
    end
end 
% Tie
if (idmn == 16) 
    if (idef == 1) 
        a1 = 1.5;
        a2 = 1;
        a3 = 1;
    else
        a1 = input('Input a1 : ');
        a2 = input('Input a2 : ');   
        a3 = input('Input a3 : ');
    end 
end
% Shell
if (idmn == 17) 
    if (idef == 1) 
        a1 = 0.5;
        a2 = 1;
    else
        a1 = input('Input a1 : ');
        a2 = input('Input a2 : ');  
    end 
end 
% New Horseshoe
if (idmn == 18) 
    if (idef == 1) 
        a1 = 4;
    else
        a1 = input('Input a1 : ');
    end
end
% Ring
if (idmn == 19)
    if (idef == 1)
        a1 = 1;
        a2 = 2;
    else
        a1 = input('Input radius of internal circle : ');
        a2 = input('Input radius of outer circle : ');
    end 
end

% NACA0012 O-type
if (idmn == 20)
    if (idef == 1)
        a1 = 4;
    else
        a1 = input('Input radius of outer circle : ');
    end
end

% NACA0012 C-type
if (idmn == 21)
    if (idef == 1)
        a1 = 2;
        a2 = 4;
        b1 = 0.4;
    else
        a1 = input('Input radius of left circle : ');
        a2 = input('Input length of right rectangle : ');
        b1 = input('Input the proportion of number of elements on the airfoil : ');
    end
end

% Slope 
if (idmn == 22)
    a1 = 0; b1 = 0;
    if (idef == 1)
        a2 = 1;
        a3 = 2;
        b3 = 0.5;
        b4 = 1;
    else
        a2 = input('Input a2 : ');  
        a3 = input('Input a3 : ');
        b3 = input('Input b3 : ');
        b4 = input('Input b4 : ');
    end
end

end