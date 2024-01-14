% Prompts user to select boundary stretch:
%  (i) uniform,
% (ii) quadratic,
%(iii) exponential
function [tb, tt, tl, tr] = ts(m, n, it)

tb = linspace(0, 1, m + 1);
tt = tb;
tl = linspace(0, 1, n + 1);
tr = tl;

if isempty(it)
    idf = input('Default parameterizations on all sides? (1=yes) ');
    if (idf == 1) 
        return
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          bottom boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(it)
    disp('(0) uniform')
    disp('(1) quadratic')
    disp('(2) stretched to one side')
    disp('(3) stretched to both sides')
    itb = input('For bottom boundary, itb = ');
else
    itb = it(1);
end
if (itb == 1)
    again = true;
    while (again)
        xi0 = input('For bottom boundary, input xi0 : ');
        q0 = input('For bottom boundary, input q0 : ');
        r0 = input('For bottom boundary, input r0 : ');
        if (xi0 <= 0) || (xi0 >= 1) 
            disp('Improper value of xi0. Try again. ')
            continue
        end
        if (q0 <= 0) || (q0 >= 1) 
            disp('Improper value of q0. Try again. ')
            continue
        end   
        if (r0 < 0)
            disp('Negative slope not allowed. Try again. ')
            continue
        end
        % monotonicity condition
        nu1 = q0 / xi0;
        nu2 = (1 - q0) / (1 - xi0);
        rmax = 2 * max(nu1, nu2);
        if (r0 >= rmax) 
            disp('Choice of parameters gives non-monotonic weight. Try again. ')
            continue
        end   
        again = false;
    end            
    ind = find(tb > xi0, 1);
    mu1 = tb(1 : ind - 1) / xi0;
    mu2 = (1 - tb(ind : end)) / (1 - xi0);
    tb(1 : ind - 1) = tb(1 : ind - 1) .* (nu1 * (2 - mu1) - r0 * (1 - mu1));
    tb(ind : end) = 1 - (1 - tb(ind : end)) .* (nu2 * (2 - mu2) - r0 * (1 - mu2));
end
if (itb == 2)
    again = true;
    while (again)
        xlam = input('For bottom boundary, input lambda : ');
        if (xlam == 0) || (xlam > 10)
            disp('Invalid input, try again. ')
            continue
        end
        again = false;
    end
    tb = (exp(xlam * tb) - 1) / (exp(xlam) - 1);
end
if (itb == 3)
    again = true;
    while (again)
        xlam = input('For bottom boundary, input lambda : ');
        if (xlam == 0) || (abs(xlam) > 10)
            disp('Invalid input, try again. ')
            continue
        end
        again = false;
    end
    tb = 0.5 * (1 + tanh(xlam * (2 * tb - 1)) / tanh(xlam));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          top boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(it)
    disp('(0) uniform')
    disp('(1) quadratic')
    disp('(2) stretched to one side')
    disp('(3) stretched to both sides')
    itt = input('For top boundary, itt = ');
else
    itt = it(2);
end
if (itt == 1)
    again = true;
    while (again)
        xi0 = input('For top boundary, input xi0 : ');
        q0 = input('For top boundary, input q0 : ');
        r0 = input('For top boundary, input r0 : ');
        if (xi0 <= 0) || (xi0 >= 1) 
            disp('Improper value of xi0. Try again. ')
            continue
        end
        if (q0 <= 0) || (q0 >= 1) 
            disp('Improper value of q0. Try again. ')
            continue
        end   
        if (r0 < 0)
            disp('Negative slope not allowed. Try again. ')
            continue
        end
        % monotonicity condition
        nu1 = q0 / xi0;
        nu2 = (1 - q0) / (1 - xi0);
        rmax = 2 * max(nu1, nu2);
        if (r0 >= rmax) 
            disp('Choice of parameters gives non-monotonic weight. Try again. ')
            continue
        end   
        again = false;
    end            
    ind = find(tt > xi0, 1);
    mu1 = tt(1 : ind - 1) / xi0;
    mu2 = (1 - tt(ind : end)) / (1 - xi0);
    tt(1 : ind - 1) = tt(1 : ind - 1) .* (nu1 * (2 - mu1) - r0 * (1 - mu1));
    tt(ind : end) = 1 - (1 - tt(ind : end)) .* (nu2 * (2 - mu2) - r0 * (1 - mu2));
end
if (itt == 2)
    again = true;
    while (again)
        xlam = input('For top boundary, input lambda : ');
        if (xlam == 0) || (xlam > 10)
            disp('Invalid input, try again. ')
            continue
        end
        again = false;
    end
    tt = (exp(xlam * tt) - 1) / (exp(xlam) - 1);
end
if (itt == 3)
    again = true;
    while (again)
        xlam = input('For top boundary, input lambda : ');
        if (xlam == 0) || (abs(xlam) > 10)
            disp('Invalid input, try again. ')
            continue
        end
        again = false;
    end
    tt = 0.5 * (1 + tanh(xlam * (2 * tt - 1)) / tanh(xlam));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          left boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(it)
    disp('(0) uniform')
    disp('(1) quadratic')
    disp('(2) stretched to one side')
    disp('(3) stretched to both sides')
    itl = input('For left boundary, itl = ');
else
    itl = it(3);
end
if (itl == 1)
    again = true;
    while (again)
        eta0 = input('For left boundary, input eta0 : ');
        q0 = input('For left boundary, input q0 : ');
        r0 = input('For left boundary, input r0 : ');
        if (eta0 <= 0) || (eta0 >= 1) 
            disp('Improper value of eta0. Try again. ')
            continue
        end
        if (q0 <= 0) || (q0 >= 1) 
            disp('Improper value of q0. Try again. ')
            continue
        end   
        if (r0 < 0)
            disp('Negative slope not allowed. Try again. ')
            continue
        end
        % monotonicity condition
        nu1 = q0 / eta0;
        nu2 = (1 - q0) / (1 - eta0);
        rmax = 2 * max(nu1, nu2);
        if (r0 >= rmax) 
            disp('Choice of parameters gives non-monotonic weight. Try again. ')
            continue
        end   
        again = false;
    end            
    ind = find(tl > eta0, 1);
    mu1 = tl(1 : ind - 1) / eta0;
    mu2 = (1 - tl(ind : end)) / (1 - eta0);
    tl(1 : ind - 1) = tl(1 : ind - 1) .* (nu1 * (2 - mu1) - r0 * (1 - mu1));
    tl(ind : end) = 1 - (1 - tl(ind : end)) .* (nu2 * (2 - mu2) - r0 * (1 - mu2));
end
if (itl == 2)
    again = true;
    while (again)
        ylam = input('For left boundary, input lambda : ');
        if (ylam == 0) || (ylam > 10)
            disp('Invalid input, try again. ')
            continue
        end
        again = false;
    end
    tl = (exp(ylam * tl) - 1) / (exp(ylam) - 1);
end
if (itl == 3)
    again = true;
    while (again)
        ylam = input('For left boundary, input lambda : ');
        if (ylam == 0) || (abs(ylam) > 10)
            disp('Invalid input, try again. ')
            continue
        end
        again = false;
    end
    tl = 0.5 * (1 + tanh(ylam * (2 * tl - 1)) / tanh(ylam));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          right boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(it)
    disp('(0) uniform')
    disp('(1) quadratic')
    disp('(2) stretched to one side')
    disp('(3) stretched to both sides')
    itr = input('For right boundary, itr = ');
else
    itr = it(4);
end
if (itr == 1)
    again = true;
    while (again)
        eta0 = input('For right boundary, input eta0 : ');
        q0 = input('For right boundary, input q0 : ');
        r0 = input('For right boundary, input r0 : ');
        if (eta0 <= 0) || (eta0 >= 1) 
            disp('Improper value of eta0. Try again. ')
            continue
        end
        if (q0 <= 0) || (q0 >= 1) 
            disp('Improper value of q0. Try again. ')
            continue
        end   
        if (r0 < 0)
            disp('Negative slope not allowed. Try again. ')
            continue
        end
        % monotonicity condition
        nu1 = q0 / eta0;
        nu2 = (1 - q0) / (1 - eta0);
        rmax = 2 * max(nu1, nu2);
        if (r0 >= rmax) 
            disp('Choice of parameters gives non-monotonic weight. Try again. ')
            continue
        end   
        again = false;
    end            
    ind = find(tr > eta0, 1);
    mu1 = tr(1 : ind - 1) / eta0;
    mu2 = (1 - tr(ind : end)) / (1 - eta0);
    tr(1 : ind - 1) = tr(1 : ind - 1) .* (nu1 * (2 - mu1) - r0 * (1 - mu1));
    tr(ind : end) = 1 - (1 - tr(ind : end)) .* (nu2 * (2 - mu2) - r0 * (1 - mu2));
end
if (itr == 2)
    again = true;
    while (again)
        ylam = input('For right boundary, input lambda : ');
        if (ylam == 0) || (ylam > 10)
            disp('Invalid input, try again. ')
            continue
        end
        again = false;
    end
    tr = (exp(ylam * tr) - 1) / (exp(ylam) - 1);
end
if (itr == 3)
    again = true;
    while (again)
        ylam = input('For right boundary, input lambda : ');
        if (ylam == 0) || (abs(ylam) > 10)
            disp('Invalid input, try again. ')
            continue
        end
        again = false;
    end
    tr = 0.5 * (1 + tanh(ylam * (2 * tr - 1)) / tanh(ylam));
end


end