function lickLatOpto(mgbTempTestsOnly,tempTestsOnly,allDataTestsOnly)

SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;

% tests for MGB by session
nbsubj=2;
for nbsubj=2:size(mgbTempTestsOnly,1)
    matrix=allDataTestsOnly{nbsubj,26};
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
        case 2
            mgbDays = [1 1 0 0 1 1 0 0 1 1 1 0 0];mgbDays=logical(mgbDays);
            icDays = [0 0 1 1 0 0 1 1 0 0 0 1 1];icDays=logical(icDays);
            expRange=1:8; 
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
        case 3
            expRange=22:32;
            mgbDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 1 1 0 0 1 1 0 0 ...
                0 0 0];mgbDays=logical(mgbDays);
            icDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                1 1 1 0 0 1 1 0 0 1 1 ...
                0 0 0];icDays=logical(icDays);
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
        case 4
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 ...
                1 0 1 1 0 1 1 0 1 0 0 1];
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 ... 
                0 1 0 0 1 0 0 1 0 1 1 0];
            icDays=logical(icDays);
            expRange=20:length(mgbDays); 
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
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

%for IC by session
for nbsubj=2:size(tempTestsOnly,1)
    matrix=tempTestsOnly{nbsubj,26};
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
        case 2
            mgbDays = [1 1 0 0 1 1 0 0 1 1 1 0 0];mgbDays=logical(mgbDays);
            icDays = [0 0 1 1 0 0 1 1 0 0 0 1 1];icDays=logical(icDays);
            expRange=1:8; 
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
        case 3
            expRange=22:32;
            mgbDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 1 1 0 0 1 1 0 0 ...
                0 0 0];mgbDays=logical(mgbDays);
            icDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                1 1 1 0 0 1 1 0 0 1 1 ...
                0 0 0];icDays=logical(icDays);
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
        case 4
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 ...
                1 0 1 1 0 1 1 0 1 0 0 1];
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 ... 
                0 1 0 0 1 0 0 1 0 1 1 0];
            icDays=logical(icDays);
            expRange=20:length(mgbDays); 
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
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

% plots for MGB
jj=2;     
for jj=2:4
    lickFig=figure(jj+103); 
    subplot(1,3,1); % full trial MGB
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
    title([allDataTestsOnly{jj,1} ' MGB Full Trial Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(1,3,2) % tone MGB
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

    subplot(1,3,3); % choice MGB
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
    lickFig.Position(3:4)=[725 275];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_LickLat_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_LickLat_Opto.png']);
    
    % LATENCY
    lickFigl=figure(jj+13); subplot(1,3,1); % full trial MGB
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

    subplot(1,3,2) % tone MGB
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

    subplot(1,3,3); % choice MGB
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
    
% now do the IC
    
    
    
    
    lickFigl.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_LickRate_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_LickRate_Opto.png']);    
end

end