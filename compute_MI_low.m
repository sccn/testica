load -mat pmi_stat_save_raw_oct2011
n = 71;
for iDat = 1:length(MIraw)
    MImat = MIraw{iDat};
    for k = 1:n, MImat(k,k) = 0; end;
    MI_low_bnd = 0;
    for k = 1:n
        [mimax, maxcol] = max(max(MImat));
        MI_low_bnd = MI_low_bnd + mimax;
        MImat(:,maxcol) = 0;
    end
    MI_low_bnd_save(iDat) = MI_low_bnd * 1.4429 * 250 / 1000;
end;
