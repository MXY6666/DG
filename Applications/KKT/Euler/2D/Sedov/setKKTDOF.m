function KKT = setKKTDOF(KKT, msh, quad2, bs)

KKT.NDOF = msh.nLElems * bs.nb;
if (KKT.enable)
    KKT.MDOF = KKT.conservation * msh.nLElems;
    KKT.LDOF = KKT.positivity * msh.nLElems * quad2.np;
else
    KKT.MDOF = 0;
    KKT.LDOF = 0;
end

end