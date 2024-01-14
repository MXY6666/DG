function KKT = setKKTMin(KKT, Umin)

if (KKT.positivity)
    KKT.Umin = Umin;
else
    KKT.Umin = -1.e+20;
end    

end