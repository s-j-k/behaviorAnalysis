function lickLatOptoPerAnimalDead(allDataTestsOnly,reinfcolor,optocolor)

SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;

% tests for MGB for each animal
nbsubj=2;
for nbsubj=2:size(allDataTestsOnly,1)
    matrix=allDataTestsOnly{nbsubj,15};
    for i=1:max(matrix(:,SESS))
        mr_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); % reinf hit
        sr_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); 
        mp_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
        sp_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
        mo_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
        so_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
        mto_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==1,LICKL));
        sto_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==1,LICKL));
        mco_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==1,LICKL));
        sco_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==1,LICKL));
        
        mr_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
        sr_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
        mp_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
        sp_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
        mo_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
        so_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
        mto_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==3,LICKL));
        sto_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==3,LICKL));
        mco_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==3,LICKL));
        sco_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==3,LICKL));
        
        mlr_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        slr_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        mlp_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        slp_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        mlo_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        slo_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        mlto_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
        slto_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
        mlco_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
        slco_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
                        
        mlr_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        slr_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        mlp_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        slp_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        mlo_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        slo_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        mlto_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
        slto_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
        mlco_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
        slco_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));        
        
    end
    switch nbsubj
        case 2 % sk198
            mgbDays=[1 1 1 1 1];mgbDays=logical(mgbDays);
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
            
        case 3 %sk203
            mgbDays=[1 1 1 1 1];mgbDays=logical(mgbDays); 
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            MGBmr_hit(MGBmr_hit==0)=NaN;
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
            
        case 4 % sk204
            mgbDays=[1 1 1 1 1];mgbDays=logical(mgbDays);
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            MGBmr_hit(MGBmr_hit==0)=NaN;
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
                        
    end
end

MGBmr_hit(MGBmr_hit==0)=NaN;
MGBsr_hit(MGBsr_hit==0)=NaN;
MGBmp_hit(MGBmp_hit==0)=NaN;
MGBsp_hit(MGBsp_hit==0)=NaN;
MGBmo_hit(MGBmo_hit==0)=NaN;
MGBso_hit(MGBso_hit==0)=NaN;
MGBmto_hit(MGBmto_hit==0)=NaN;
MGBsto_hit(MGBsto_hit==0)=NaN;
MGBmco_hit(MGBmco_hit==0)=NaN;
MGBsco_hit(MGBsco_hit==0)=NaN;

MGBmr_fa(MGBmr_fa==0) = NaN;
MGBsr_fa(MGBsr_fa==0) = NaN;
MGBmp_fa(MGBmp_fa==0) = NaN;
MGBsp_fa(MGBsp_fa==0) = NaN;
MGBmo_fa(MGBmo_fa==0) = NaN;
MGBso_fa(MGBso_fa==0) = NaN;
MGBmto_fa(MGBmto_fa==0) = NaN;
MGBsto_fa(MGBsto_fa==0) = NaN;
MGBmco_fa(MGBmco_fa==0) = NaN;
MGBsco_fa(MGBsco_fa==0) = NaN;

MGBmlr_hit(MGBmlr_hit==0) = NaN;
MGBslr_hit(MGBslr_hit==0) = NaN;
MGBmlp_hit(MGBmlp_hit==0) = NaN;
MGBslp_hit(MGBslp_hit==0) = NaN;
MGBmlo_hit(MGBmlo_hit==0) = NaN;
MGBslo_hit(MGBslo_hit==0) = NaN;
MGBmlto_hit(MGBmlto_hit==0) = NaN;
MGBslto_hit(MGBslto_hit==0) = NaN;
MGBmlco_hit(MGBmlco_hit==0) = NaN;
MGBslco_hit(MGBslco_hit==0) = NaN;

MGBmlr_fa(MGBmlr_fa==0) = NaN;
MGBslr_fa(MGBslr_fa==0) = NaN;
MGBmlp_fa(MGBmlp_fa==0) = NaN;
MGBslp_fa(MGBslp_fa==0) = NaN;
MGBmlo_fa(MGBmlo_fa==0) = NaN;
MGBslo_fa(MGBslo_fa==0) = NaN;
MGBmlto_fa(MGBmlto_fa==0) = NaN;
MGBslto_fa(MGBslto_fa==0) = NaN;
MGBmlco_fa(MGBmlco_fa==0) = NaN;
MGBslco_fa(MGBslco_fa==0) = NaN;

ICmr_hit(ICmr_hit==0)=NaN;
ICsr_hit(ICsr_hit==0)=NaN;
ICmp_hit(ICmp_hit==0)=NaN;
ICsp_hit(ICsp_hit==0)=NaN;
ICmo_hit(ICmo_hit==0)=NaN; 
ICso_hit(ICso_hit==0)=NaN; 
ICmto_hit(ICmto_hit==0)=NaN;
ICsto_hit(ICsto_hit==0)=NaN; 
ICmco_hit(ICmco_hit==0)=NaN; 
ICsco_hit(ICsco_hit==0)=NaN;

ICmr_fa(ICmr_fa==0) = NaN;
ICsr_fa(ICsr_fa==0) = NaN;
ICmp_fa(ICmp_fa==0) = NaN;
ICsp_fa(ICsp_fa==0) = NaN;
ICmo_fa(ICmo_fa==0) = NaN;
ICso_fa(ICso_fa==0) = NaN;
ICmto_fa(ICmto_fa==0) = NaN;
ICsto_fa(ICsto_fa==0) = NaN;
ICmco_fa(ICmco_fa==0) = NaN;
ICsco_fa(ICsco_fa==0) = NaN;

ICmlr_hit(ICmlr_hit==0) = NaN;
ICslr_hit(ICslr_hit==0) = NaN;
ICmlp_hit(ICmlp_hit==0) = NaN;
ICslp_hit(ICslp_hit==0) = NaN;
ICmlo_hit(ICmlo_hit==0) = NaN;
ICslo_hit(ICslo_hit==0) = NaN;
ICmlto_hit(ICmlto_hit==0) = NaN;
ICslto_hit(ICslto_hit==0) = NaN;
ICmlco_hit(ICmlco_hit==0) = NaN;
ICslco_hit(ICslco_hit==0) = NaN;

ICmlr_fa(ICmlr_fa==0) = NaN;
ICslr_fa(ICslr_fa==0) = NaN;
ICmlp_fa(ICmlp_fa==0) = NaN;
ICslp_fa(ICslp_fa==0) = NaN;
ICmlo_fa(ICmlo_fa==0) = NaN;
ICslo_fa(ICslo_fa==0) = NaN;
ICmlto_fa(ICmlto_fa==0) = NaN;
ICslto_fa(ICslto_fa==0) = NaN;
ICmlco_fa(ICmlco_fa==0) = NaN;
ICslco_fa(ICslco_fa==0) = NaN; 


% plots 
jj=2;     
for jj=2:size(allDataTestsOnly,1)
    lickFig=figure(jj+103); 
    subplot(2,3,1); hold on; %MGB FUll
    lickbar=bar([nanmean(MGBmr_hit(jj,:)) nanmean(MGBmo_hit(jj,:)) nanmean(MGBmr_fa(jj,:)) nanmean(MGBmo_fa(jj,:))]); hold on;
    error=[nanmean(MGBsr_hit(jj,:)) nanmean(MGBso_hit(jj,:)) nanmean(MGBsr_fa(jj,:)) nanmean(MGBso_fa(jj,:))];
    lickbar(1).FaceColor='flat'; lickbar(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    for kk=1:numel(lickbar)
        xtips=lickbar(kk).XEndPoints;
        ytips=lickbar(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(MGBmr_hit(jj,:),MGBmo_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(MGBmr_fa(jj,:),MGBmo_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick latency');
    title([allDataTestsOnly{jj,1} ' MGB Dead 1 Inactivation']);
    xticklabels({'hit', 'fa','hit','fa'});

    subplot(2,3,2) % tone MGB
    lickbarT=bar([nanmean(MGBmr_hit(jj,:)) nanmean(MGBmto_hit(jj,:)) nanmean(MGBmr_fa(jj,:)) nanmean(MGBmto_fa(jj,:))]); hold on;
    error=[nanmean(MGBsr_hit(jj,:)) mean(MGBsto_hit(jj,:)) nanmean(MGBsr_fa(jj,:)) nanmean(MGBsto_fa(jj,:))];
    lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    for kk=1:numel(lickbarT)
        xtips=lickbarT(kk).XEndPoints;
        ytips=lickbarT(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(MGBmr_hit(jj,:),MGBmto_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(MGBmr_fa(jj,:),MGBmto_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick latency');
    title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(2,3,3); % choice MGB
    lickbarC=bar([nanmean(MGBmr_hit(jj,:)) nanmean(MGBmco_hit(jj,:)) nanmean(MGBmr_fa(jj,:)) nanmean(MGBmco_fa(jj,:))]); hold on;
    lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(MGBsr_hit(jj,:)) nanmean(MGBsco_hit(jj,:)) nanmean(MGBsr_fa(jj,:)) nanmean(MGBsco_fa(jj,:))];
    for kk=1:numel(lickbarC)
        xtips=lickbarC(kk).XEndPoints;
        ytips=lickbarC(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(MGBmr_hit(jj,:),MGBmco_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(MGBmr_fa(jj,:),MGBmco_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick latency');
    title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});


    subplot(2,3,4); hold on; %IC Full
    lickbar=bar([nanmean(ICmr_hit(jj,:)) nanmean(ICmo_hit(jj,:)) nanmean(ICmr_fa(jj,:)) nanmean(ICmo_fa(jj,:))]); hold on;
    error=[nanmean(ICsr_hit(jj,:)) nanmean(ICso_hit(jj,:)) nanmean(ICsr_fa(jj,:)) nanmean(ICso_fa(jj,:))];
    lickbar(1).FaceColor='flat'; lickbar(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    for kk=1:numel(lickbar)
        xtips=lickbar(kk).XEndPoints;
        ytips=lickbar(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(ICmr_hit(jj,:),ICmo_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(ICmr_fa(jj,:),ICmo_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick latency');
    title([allDataTestsOnly{jj,1} ' IC Full Trial Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(2,3,5) % tone IC
    lickbarT=bar([nanmean(ICmr_hit(jj,:)) nanmean(ICmto_hit(jj,:)) nanmean(MGBmr_fa(jj,:)) nanmean(ICmto_fa(jj,:))]); hold on;
    error=[nanmean(ICsr_hit(jj,:)) mean(ICsto_hit(jj,:)) nanmean(ICsr_fa(jj,:)) nanmean(ICsto_fa(jj,:))];
    lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    for kk=1:numel(lickbarT)
        xtips=lickbarT(kk).XEndPoints;
        ytips=lickbarT(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(ICmr_hit(jj,:),ICmto_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(ICmr_fa(jj,:),ICmto_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick latency');
    title([allDataTestsOnly{jj,1} ' IC Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(2,3,6); % choice IC
    lickbarC=bar([nanmean(ICmr_hit(jj,:)) nanmean(ICmco_hit(jj,:)) nanmean(ICmr_fa(jj,:)) nanmean(ICmco_fa(jj,:))]); hold on;
    lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(ICsr_hit(jj,:)) nanmean(ICsco_hit(jj,:)) nanmean(ICsr_fa(jj,:)) nanmean(ICsco_fa(jj,:))];
    for kk=1:numel(lickbarC)
        xtips=lickbarC(kk).XEndPoints;
        ytips=lickbarC(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(ICmr_hit(jj,:),ICmco_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(ICmr_fa(jj,:),ICmco_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick latency');
    title([allDataTestsOnly{jj,1} ' IC Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    lickFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_IC_LickLat_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_IC_LickLat_Opto.png']);

    % RATE
    lickFigl=figure(jj+13); subplot(2,3,1); % full trial MGB
    lickbarL=bar([nanmean(MGBmlr_hit(jj,:)) nanmean(MGBmlo_hit(jj,:)) nanmean(MGBmlr_fa(jj,:)) nanmean(MGBmlo_fa(jj,:))]); hold on;
    lickbarL(1).FaceColor='flat'; lickbarL(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(MGBslr_hit(jj,:)) nanmean(MGBslo_hit(jj,:)) nanmean(MGBslr_fa(jj,:)) nanmean(MGBslo_fa(jj,:))];
    for kk=1:numel(lickbarL)
        xtips=lickbarL(kk).XEndPoints;
        ytips=lickbarL(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(MGBmlr_hit(jj,:),MGBmlo_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(MGBmlr_fa(jj,:),MGBmlo_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title([allDataTestsOnly{jj,1} ' MGB Full Trial Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(2,3,2) % tone MGB
    lickbarT=bar([nanmean(MGBmlr_hit(jj,:)) nanmean(MGBmlto_hit(jj,:)) nanmean(MGBmlr_fa(jj,:)) nanmean(MGBmlto_fa(jj,:))]); hold on;
    lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(MGBslr_hit(jj,:)) nanmean(MGBslto_hit(jj,:)) nanmean(MGBslr_fa(jj,:)) nanmean(MGBslto_fa(jj,:))];
    for kk=1:numel(lickbarT)
        xtips=lickbarT(kk).XEndPoints;
        ytips=lickbarT(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(MGBmlr_hit(jj,:),MGBmlto_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(MGBmlr_fa(jj,:),MGBmlto_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(2,3,3); % choice MGB
    lickbarC=bar([nanmean(MGBmlr_hit(jj,:)) nanmean(MGBmlco_hit(jj,:)) nanmean(MGBmlr_fa(jj,:)) nanmean(MGBmlco_fa(jj,:))]); hold on;
    lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(MGBslr_hit(jj,:)) nanmean(MGBslco_hit(jj,:)) nanmean(MGBslr_fa(jj,:)) nanmean(MGBslco_fa(jj,:))];
    for kk=1:numel(lickbarC)
        xtips=lickbarC(kk).XEndPoints;
        ytips=lickbarC(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(MGBmlr_hit(jj,:),MGBmlco_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(MGBmlr_fa(jj,:),MGBmlco_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(2,3,4); % full trial IC
    lickbar=bar([nanmean(ICmlr_hit(jj,:)) nanmean(ICmlo_hit(jj,:)) nanmean(ICmlr_fa(jj,:)) nanmean(ICmlo_fa(jj,:))]); hold on;
    error=[nanmean(ICslr_hit(jj,:)) nanmean(ICslo_hit(jj,:)) nanmean(ICslr_fa(jj,:)) nanmean(ICslo_fa(jj,:))];
    lickbar(1).FaceColor='flat'; lickbar(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    for kk=1:numel(lickbar)
        xtips=lickbar(kk).XEndPoints;
        ytips=lickbar(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end

    [h,pHit,ci,stats] = ttest2(ICmlr_hit(jj,:),ICmlo_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(ICmlr_fa(jj,:),ICmlo_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title([allDataTestsOnly{jj,1} ' IC Full Trial Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(2,3,5) % tone IC
    lickbarT=bar([nanmean(ICmlr_hit(jj,:)) nanmean(ICmlto_hit(jj,:)) nanmean(ICmlr_fa(jj,:)) nanmean(ICmlto_fa(jj,:))]); hold on;
    error=[nanmean(ICslr_hit(jj,:)) mean(ICslto_hit(jj,:)) nanmean(ICslr_fa(jj,:)) nanmean(ICslto_fa(jj,:))];
    lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    for kk=1:numel(lickbarT)
        xtips=lickbarT(kk).XEndPoints;
        ytips=lickbarT(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(ICmlr_hit(jj,:),ICmlto_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(ICmlr_fa(jj,:),ICmlto_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title([allDataTestsOnly{jj,1} ' IC Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(2,3,6); % choice IC
    lickbarC=bar([nanmean(ICmlr_hit(jj,:)) nanmean(ICmlco_hit(jj,:)) nanmean(ICmlr_fa(jj,:)) nanmean(ICmlco_fa(jj,:))]); hold on;
    lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(ICslr_hit(jj,:)) nanmean(ICslco_hit(jj,:)) nanmean(ICslr_fa(jj,:)) nanmean(ICslco_fa(jj,:))];
    for kk=1:numel(lickbarC)
        xtips=lickbarC(kk).XEndPoints;
        ytips=lickbarC(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(ICmlr_hit(jj,:),ICmlco_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(ICmlr_fa(jj,:),ICmlco_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title([allDataTestsOnly{jj,1} ' IC Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});
    
    lickFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_IC_LickRate_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_IC_LickRate_Opto.png']);  
end


close all
% LICK LATENCY, AVERAGED ACROSS ANIMALS
% change these variable names, theyre swapped
lickFig=figure(jj+103); 
for jj=2:size(mgbTempTestsOnly,1)
    hitLickRate(jj-1)=nanmean(MGBmr_hit(jj,:));
    ohitLickRate(jj-1)=nanmean(MGBmo_hit(jj,:));
    faLickRate(jj-1)=nanmean(MGBmr_fa(jj,:));
    ofaLickRate(jj-1)=nanmean(MGBmo_fa(jj,:));
    othitLickRate(jj-1)=nanmean(MGBmto_hit(jj,:));
    otfaLickRate(jj-1)=nanmean(MGBmto_fa(jj,:));
    ochitLickRate(jj-1)=nanmean(MGBmco_hit(jj,:));
    ocfaLickRate(jj-1)=nanmean(MGBmco_fa(jj,:));
    
    hitLickRateE(jj-1)=nanmean(MGBsr_hit(jj,:));
    ohitLickRateE(jj-1)=nanmean(MGBso_hit(jj,:));
    faLickRateE(jj-1)=nanmean(MGBsr_fa(jj,:));
    ofaLickRateE(jj-1)=nanmean(MGBso_fa(jj,:));
    othitLickRateE(jj-1)=nanmean(MGBsto_hit(jj,:));
    otfaLickRateE(jj-1)=nanmean(MGBsto_fa(jj,:));
    ochitLickRateE(jj-1)=nanmean(MGBsco_hit(jj,:));
    ocfaLickRateE(jj-1)=nanmean(MGBsco_fa(jj,:));
    
end

subplot(2,3,1); % full trial MGB
lickbar=bar([nanmean(hitLickRate) nanmean(ohitLickRate) nanmean(faLickRate) nanmean(ofaLickRate)]); hold on;
error=[nanmean(hitLickRateE) nanmean(ohitLickRateE) nanmean(faLickRateE) nanmean(ofaLickRateE)];
lickbar(1).FaceColor='flat'; lickbar(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
for kk=1:numel(lickbar)
    xtips=lickbar(kk).XEndPoints;
    ytips=lickbar(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,ohitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,ofaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick latency');
title(['By Animal MGB Full Trial']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,2) % tone MGB
lickbarT=bar([nanmean(hitLickRate) nanmean(othitLickRate) nanmean(faLickRate) nanmean(otfaLickRate)]); hold on;
error=[nanmean(hitLickRateE) nanmean(othitLickRateE) nanmean(faLickRateE) nanmean(otfaLickRateE)];
lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
for kk=1:numel(lickbarT)
    xtips=lickbarT(kk).XEndPoints;
    ytips=lickbarT(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,othitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,otfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick latency');
title(['By Animal MGB Tone']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,3); % choice MGB
lickbarC=bar([nanmean(hitLickRate) nanmean(ochitLickRate) nanmean(faLickRate) nanmean(ocfaLickRate)]); hold on;
lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(hitLickRateE) nanmean(ochitLickRateE) nanmean(faLickRateE) nanmean(ocfaLickRateE)];
for kk=1:numel(lickbarC)
    xtips=lickbarC(kk).XEndPoints;
    ytips=lickbarC(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,ochitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,ocfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick latency');
title(['By Animal MGB Choice']);
xticklabels({'hit', 'hit','fa','fa'});

icIdx=[3,5,6,7];icIdx=icIdx-1;
clear hitLickRate ohitLickRate faLickRate ofaLickRate othitLickRate otfaLickRate ochitLickRate ocfaLickRate
for jj=2:length(icIdx)
    hitLickRate(jj)=nanmean(ICmr_hit(icIdx(jj),:));
    ohitLickRate(jj)=nanmean(ICmo_hit(icIdx(jj),:));
    faLickRate(jj)=nanmean(ICmr_fa(icIdx(jj),:));
    ofaLickRate(jj)=nanmean(ICmo_fa(icIdx(jj),:));
    othitLickRate(jj)=nanmean(ICmto_hit(icIdx(jj),:));
    otfaLickRate(jj)=nanmean(ICmto_fa(icIdx(jj),:));
    ochitLickRate(jj)=nanmean(ICmco_hit(icIdx(jj),:));
    ocfaLickRate(jj)=nanmean(ICmco_fa(icIdx(jj),:));
    
    hitLickRateE(jj)=nanmean(ICsr_hit(icIdx(jj),:));
    ohitLickRateE(jj)=nanmean(ICso_hit(icIdx(jj),:));
    faLickRateE(jj)=nanmean(ICsr_fa(icIdx(jj),:));
    ofaLickRateE(jj)=nanmean(ICso_fa(icIdx(jj),:));
    othitLickRateE(jj)=nanmean(ICsto_hit(icIdx(jj),:));
    otfaLickRateE(jj)=nanmean(ICsto_fa(icIdx(jj),:));
    ochitLickRateE(jj)=nanmean(ICsco_hit(icIdx(jj),:));
    ocfaLickRateE(jj)=nanmean(ICsco_fa(icIdx(jj),:));
    
end

subplot(2,3,4); % full trial IC
lickbar=bar([nanmean(hitLickRate) nanmean(ohitLickRate) nanmean(faLickRate) nanmean(ofaLickRate)]); hold on;
error=[nanmean(hitLickRateE) nanmean(ohitLickRateE) nanmean(faLickRateE) nanmean(ofaLickRateE)];
lickbar(1).FaceColor='flat'; lickbar(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
for kk=1:numel(lickbar)
    xtips=lickbar(kk).XEndPoints;
    ytips=lickbar(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,ohitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,ofaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick latency');
title(['By Animal IC Full Trial']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,5) % tone IC
lickbarT=bar([nanmean(hitLickRate) nanmean(othitLickRate) nanmean(faLickRate) nanmean(otfaLickRate)]); hold on;
error=[nanmean(hitLickRateE) mean(othitLickRateE) nanmean(faLickRateE) nanmean(otfaLickRateE)];
lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
for kk=1:numel(lickbarT)
    xtips=lickbarT(kk).XEndPoints;
    ytips=lickbarT(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,othitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,otfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick latency');
title(['By Animal IC Tone']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,6); % choice IC
lickbarC=bar([nanmean(hitLickRate) nanmean(ochitLickRate) nanmean(faLickRate) nanmean(ocfaLickRate)]); hold on;
lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(hitLickRateE) nanmean(ochitLickRateE) nanmean(faLickRateE) nanmean(ocfaLickRateE)];
for kk=1:numel(lickbarC)
    xtips=lickbarC(kk).XEndPoints;
    ytips=lickbarC(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,ochitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,ocfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick latency');
title(['By Animal IC Choice']);
xticklabels({'hit', 'hit','fa','fa'});

lickFig.Position(3:4)=[725 475];
saveas(gcf,['ByAnimal_T_MGB_IC_LickLat_Opto']);
saveas(gcf,['ByAnimal_T_MGB_IC_LickLat_Opto.png']);
close all
% RATE
% need to edit this part! 10/18/24
lickFigl=figure(jj+13); 
for jj=2:size(mgbTempTestsOnly,1)
    hitLickRate(jj-1)=nanmean(MGBmlr_hit(jj,:));
    ohitLickRate(jj-1)=nanmean(MGBmlo_hit(jj,:));
    faLickRate(jj-1)=nanmean(MGBmlr_fa(jj,:));
    ofaLickRate(jj-1)=nanmean(MGBmlo_fa(jj,:));
    othitLickRate(jj-1)=nanmean(MGBmlto_hit(jj,:));
    otfaLickRate(jj-1)=nanmean(MGBmlto_fa(jj,:));
    ochitLickRate(jj-1)=nanmean(MGBmlco_hit(jj,:));
    ocfaLickRate(jj-1)=nanmean(MGBmlco_fa(jj,:));
    
    hitLickRateE(jj-1)=nanmean(MGBslo_hit(jj,:));
    ohitLickRateE(jj-1)=nanmean(MGBslo_hit(jj,:));
    faLickRateE(jj-1)=nanmean(MGBslr_fa(jj,:));
    ofaLickRateE(jj-1)=nanmean(MGBslo_fa(jj,:));
    othitLickRateE(jj-1)=nanmean(MGBslto_hit(jj,:));
    otfaLickRateE(jj-1)=nanmean(MGBslto_fa(jj,:));
    ochitLickRateE(jj-1)=nanmean(MGBslco_hit(jj,:));
    ocfaLickRateE(jj-1)=nanmean(MGBslco_fa(jj,:));    
    
end

subplot(2,3,1); % full trial MGB
lickbarL=bar([nanmean(hitLickRate) nanmean(ohitLickRate) nanmean(faLickRate) nanmean(ofaLickRate)]); hold on;
lickbarL(1).FaceColor='flat'; lickbarL(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(hitLickRateE) nanmean(ohitLickRateE) nanmean(faLickRateE) nanmean(ofaLickRateE)];
for kk=1:numel(lickbarL)
    xtips=lickbarL(kk).XEndPoints;
    ytips=lickbarL(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,ohitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,ofaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick rate');
title(['By Animal MGB Full Trial']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,2) % tone MGB
lickbarT=bar([nanmean(hitLickRate) nanmean(othitLickRate) nanmean(faLickRate) nanmean(otfaLickRate)]); hold on;
lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(hitLickRateE) nanmean(othitLickRateE) nanmean(faLickRateE) nanmean(otfaLickRateE)];
for kk=1:numel(lickbarT)
    xtips=lickbarT(kk).XEndPoints;
    ytips=lickbarT(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,othitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,otfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick rate');
title(['By Animal MGB Tone']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,3); % choice MGB
lickbarC=bar([nanmean(hitLickRate) nanmean(ochitLickRate) nanmean(faLickRate) nanmean(ocfaLickRate)]); hold on;
lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(hitLickRateE) nanmean(ochitLickRateE) nanmean(faLickRateE) nanmean(ocfaLickRateE)];
for kk=1:numel(lickbarC)
    xtips=lickbarC(kk).XEndPoints;
    ytips=lickbarC(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,ochitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,ocfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick rate');
title(['By Animal MGB Choice']);
xticklabels({'hit', 'hit','fa','fa'});

% do the IC

icIdx=[3,5,6,7];
for jj=1:length(icIdx)
    hitLickRate(jj)=nanmean(ICmlo_hit(icIdx(jj),:));
    ohitLickRate(jj)=nanmean(ICmlo_hit(icIdx(jj),:));
    faLickRate(jj)=nanmean(ICmlr_fa(icIdx(jj),:));
    ofaLickRate(jj)=nanmean(ICmlo_fa(icIdx(jj),:));
    othitLickRate(jj)=nanmean(ICmlto_hit(icIdx(jj),:));
    otfaLickRate(jj)=nanmean(ICmlto_fa(icIdx(jj),:));
    ochitLickRate(jj)=nanmean(ICmlco_hit(icIdx(jj),:));
    ocfaLickRate(jj)=nanmean(ICmlco_fa(icIdx(jj),:));
    %error
    hitLickRateE(jj)=nanmean(ICslo_hit(icIdx(jj),:));
    ohitLickRateE(jj)=nanmean(ICslo_hit(icIdx(jj),:));
    faLickRateE(jj)=nanmean(ICslr_fa(icIdx(jj),:));
    ofaLickRateE(jj)=nanmean(ICslo_fa(icIdx(jj),:));
    othitLickRateE(jj)=nanmean(ICslto_hit(icIdx(jj),:));
    otfaLickRateE(jj)=nanmean(ICslto_fa(icIdx(jj),:));
    ochitLickRateE(jj)=nanmean(ICslco_hit(icIdx(jj),:));
    ocfaLickRateE(jj)=nanmean(ICslco_fa(icIdx(jj),:));    
    
end

subplot(2,3,4); % full trial IC
lickbarL=bar([nanmean(hitLickRate) nanmean(ohitLickRate) nanmean(faLickRate) nanmean(ofaLickRate)]); hold on;
lickbarL(1).FaceColor='flat'; lickbarL(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(hitLickRateE) nanmean(ohitLickRateE) nanmean(faLickRateE) nanmean(ofaLickRateE)];
for kk=1:numel(lickbarL)
    xtips=lickbarL(kk).XEndPoints;
    ytips=lickbarL(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,ohitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,ofaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick rate');
title(['By Animal IC Full Trial']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,5) % tone IC
lickbarT=bar([nanmean(hitLickRate) nanmean(othitLickRate) nanmean(faLickRate) nanmean(otfaLickRate)]); hold on;
lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(hitLickRateE) nanmean(othitLickRateE) nanmean(faLickRateE) nanmean(otfaLickRateE)];
for kk=1:numel(lickbarT)
    xtips=lickbarT(kk).XEndPoints;
    ytips=lickbarT(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,othitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,otfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick rate');
title(['By Animal IC Tone']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,6); % choice IC
lickbarC=bar([nanmean(hitLickRate) nanmean(ochitLickRate) nanmean(faLickRate) nanmean(ocfaLickRate)]); hold on;
lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(hitLickRateE) nanmean(ochitLickRateE) nanmean(faLickRateE) nanmean(ocfaLickRateE)];
for kk=1:numel(lickbarC)
    xtips=lickbarC(kk).XEndPoints;
    ytips=lickbarC(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,ochitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,ocfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick rate');
title(['By Animal IC Choice']);
xticklabels({'hit', 'hit','fa','fa'});

lickFigl.Position(3:4)=[725 475];
saveas(gcf,['ByAnimal_T_MGB_IC_LickRate_Opto']);
saveas(gcf,['ByAnimal_T_MGB_IC_LickRate_Opto.png']);  
close all

% by session
lickFig1=figure(jj+108); 

% latency
clear hitLickRate ohitLickRate faLickRate ofaLickRate othitLickRate otfaLickRate ...
    ochitLickRate ocfaLickRate hitLickRateE ohitLickRateE faLickRateE ...
    ofaLickRateE othitLickRateE otfaLickRateE ochitLickRateE ocfaLickRateE 
hitLickRate=NaN; ohitLickRate=NaN; faLickRate=NaN; ofaLickRate=NaN; othitLickRate=NaN;
otfaLickRate=NaN; ochitLickRate=NaN; ocfaLickRate=NaN; hitLickRateE=NaN; ohitLickRateE=NaN;
faLickRateE=NaN; ofaLickRateE=NaN; othitLickRateE=NaN; otfaLickRateE=NaN; ochitLickRateE=NaN;
ocfaLickRateE=NaN; 
for jj=2:size(mgbTempTestsOnly,1)
    hitLickRate=cat(2,hitLickRate,MGBmr_hit(jj,:));
    ohitLickRate=cat(2,ohitLickRate,MGBmo_hit(jj,:));
    faLickRate=cat(2,faLickRate,MGBmr_fa(jj,:));
    ofaLickRate=cat(2,ofaLickRate,MGBmo_fa(jj,:));
    othitLickRate=cat(2,othitLickRate,MGBmto_hit(jj,:));
    otfaLickRate=cat(2,otfaLickRate,MGBmto_fa(jj,:));
    ochitLickRate=cat(2,ochitLickRate,MGBmco_hit(jj,:));
    ocfaLickRate=cat(2,ocfaLickRate,MGBmco_fa(jj,:));
    
    hitLickRateE=cat(2,hitLickRateE,MGBsr_hit(jj,:));
    ohitLickRateE=cat(2,ohitLickRateE,MGBso_hit(jj,:));
    faLickRateE=cat(2,faLickRateE,MGBsr_fa(jj,:));
    ofaLickRateE=cat(2,ofaLickRateE,MGBso_fa(jj,:));
    othitLickRateE=cat(2,othitLickRateE,MGBsto_hit(jj,:));
    otfaLickRateE=cat(2,otfaLickRateE,MGBsto_fa(jj,:));
    ochitLickRateE=cat(2,ochitLickRateE,MGBsco_hit(jj,:));
    ocfaLickRateE=cat(2,ocfaLickRateE,MGBsco_fa(jj,:));
    
end

subplot(2,3,1); % full trial MGB
lickbar=bar([nanmean(hitLickRate) nanmean(ohitLickRate) nanmean(faLickRate) nanmean(ofaLickRate)]); hold on;
error=[nanmean(hitLickRateE) nanmean(ohitLickRateE) nanmean(faLickRateE) nanmean(ofaLickRateE)];
lickbar(1).FaceColor='flat'; lickbar(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
for kk=1:numel(lickbar)
    xtips=lickbar(kk).XEndPoints;
    ytips=lickbar(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,ohitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,ofaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick latency');
title(['By Animal MGB Full Trial']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,2) % tone MGB
lickbarT=bar([nanmean(hitLickRate) nanmean(othitLickRate) nanmean(faLickRate) nanmean(otfaLickRate)]); hold on;
error=[nanmean(hitLickRateE) nanmean(othitLickRateE) nanmean(faLickRateE) nanmean(otfaLickRateE)];
lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
for kk=1:numel(lickbarT)
    xtips=lickbarT(kk).XEndPoints;
    ytips=lickbarT(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,othitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,otfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick latency');
title(['By Animal MGB Tone']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,3); % choice MGB
lickbarC=bar([nanmean(hitLickRate) nanmean(ochitLickRate) nanmean(faLickRate) nanmean(ocfaLickRate)]); hold on;
lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(hitLickRateE) nanmean(ochitLickRateE) nanmean(faLickRateE) nanmean(ocfaLickRateE)];
for kk=1:numel(lickbarC)
    xtips=lickbarC(kk).XEndPoints;
    ytips=lickbarC(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,ochitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,ocfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick latency');
title(['By Animal MGB Choice']);
xticklabels({'hit', 'hit','fa','fa'});

% for rate
clear hitLickRate ohitLickRate faLickRate ofaLickRate othitLickRate otfaLickRate ...
    ochitLickRate ocfaLickRate hitLickRateE ohitLickRateE faLickRateE ...
    ofaLickRateE othitLickRateE otfaLickRateE ochitLickRateE ocfaLickRateE 
hitLickRate=NaN; ohitLickRate=NaN; faLickRate=NaN; ofaLickRate=NaN; othitLickRate=NaN;
otfaLickRate=NaN; ochitLickRate=NaN; ocfaLickRate=NaN; hitLickRateE=NaN; ohitLickRateE=NaN;
faLickRateE=NaN; ofaLickRateE=NaN; othitLickRateE=NaN; otfaLickRateE=NaN; ochitLickRateE=NaN;
ocfaLickRateE=NaN;
for jj=2:size(mgbTempTestsOnly,1)
    hitLickRate=cat(2,hitLickRate,MGBmlr_hit(jj,:));
    ohitLickRate=cat(2,ohitLickRate,MGBmlo_hit(jj,:));
    faLickRate=cat(2,faLickRate,MGBmlr_fa(jj,:));
    ofaLickRate=cat(2,ofaLickRate,MGBmlo_fa(jj,:));
    othitLickRate=cat(2,othitLickRate,MGBmlto_hit(jj,:));
    otfaLickRate=cat(2,otfaLickRate,MGBmlto_fa(jj,:));
    ochitLickRate=cat(2,ochitLickRate,MGBmlco_hit(jj,:));
    ocfaLickRate=cat(2,ocfaLickRate,MGBmlco_fa(jj,:));
    
    hitLickRateE=cat(2,hitLickRateE,MGBslr_hit(jj,:));
    ohitLickRateE=cat(2,ohitLickRateE,MGBslo_hit(jj,:));
    faLickRateE=cat(2,faLickRateE,MGBslr_fa(jj,:));
    ofaLickRateE=cat(2,ofaLickRateE,MGBslo_fa(jj,:));
    othitLickRateE=cat(2,othitLickRateE,MGBslto_hit(jj,:));
    otfaLickRateE=cat(2,otfaLickRateE,MGBslto_fa(jj,:));
    ochitLickRateE=cat(2,ochitLickRateE,MGBslco_hit(jj,:));
    ocfaLickRateE=cat(2,ocfaLickRateE,MGBslco_fa(jj,:));
end

subplot(2,3,4); % full trial MGB
lickbarL=bar([nanmean(hitLickRate) nanmean(ohitLickRate) nanmean(faLickRate) nanmean(ofaLickRate)]); hold on;
lickbarL(1).FaceColor='flat'; lickbarL(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(hitLickRateE) nanmean(ohitLickRateE) nanmean(faLickRateE) nanmean(ofaLickRateE)];
for kk=1:numel(lickbarL)
    xtips=lickbarL(kk).XEndPoints;
    ytips=lickbarL(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,ohitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,ofaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick rate');
title(['By Animal MGB Full Trial']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,5) % tone MGB
lickbarT=bar([nanmean(hitLickRate) nanmean(othitLickRate) nanmean(faLickRate) nanmean(otfaLickRate)]); hold on;
lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(hitLickRateE) nanmean(othitLickRateE) nanmean(faLickRateE) nanmean(otfaLickRateE)];
for kk=1:numel(lickbarT)
    xtips=lickbarT(kk).XEndPoints;
    ytips=lickbarT(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,othitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,otfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick rate');
title(['By Animal MGB Tone']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,6); % choice MGB
lickbarC=bar([nanmean(hitLickRate) nanmean(ochitLickRate) nanmean(faLickRate) nanmean(ocfaLickRate)]); hold on;
lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(hitLickRateE) nanmean(ochitLickRateE) nanmean(faLickRateE) nanmean(ocfaLickRateE)];
for kk=1:numel(lickbarC)
    xtips=lickbarC(kk).XEndPoints;
    ytips=lickbarC(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,ochitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,ocfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick rate');
title(['By Animal MGB Choice']);
xticklabels({'hit', 'hit','fa','fa'});


saveas(gcf,['BySess_T_MGB_LickRateLat_Opto']);
saveas(gcf,['BySess_T_MGB_LickRateLat_Opto.png']);  


end