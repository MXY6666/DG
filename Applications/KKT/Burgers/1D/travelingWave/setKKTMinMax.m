function KKT = setKKTMinMax(KKT, Umin, Umax)

if (KKT.enable)
    if (KKT.positivity)
        KKT.Umin = Umin;
    else
        KKT.Umin = -1.e+20;
    end
    
    if (KKT.maximum)
        KKT.Umax = Umax;
    else
        KKT.Umax = 1.e+20;
    end
    
else
    KKT.Umin = -1.e+20;
    KKT.Umax = 1.e+20;
end

end