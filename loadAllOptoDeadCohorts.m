function allCohorts = loadAllOptoDeadCohorts(cohortRange)
% takes input of an array which is the cohort number, and loads data from
% those cohorts specified
    for ii=1:length(cohortRange)
        if cohortRange(ii)==6
                cd('O:\sjk\Behavior\MGBIC_6')
                load('deadSummaryData.mat');
                cohort6=optomeanMat;
                
        else disp('Cohort not found')
            
        end
    end
    cd('O:\sjk\Behavior')
    
    % now make a concatenated variable of the cohorts called, and call each
    % of them to create the allCohorts variable 
    allCohorts={};
    
    % the variables have different columns but they shouldn't. need to
    % rerun the extraction for cohort1&2 to have the rates column, and
    % cohort 3
    
    %here is a shitty workaround
    cohort1(1,27)={'Rates Variable Placeholder'};
    cohort2(1,27)={'Rates Variable Placeholder'};
    cohort3(:,28)=[];
    cohort2(1,:)=[];cohort3(1,:)=[];cohort4(1,:)=[];cohort5(1,:)=[];
    cohort6(1,:)=[];
    allCohorts=vertcat(cohort1, cohort2, cohort3, cohort4, cohort5, cohort6);
    allCohorts(26:27,:)=[]; 
    allCohorts(28:31,:)=[]; 
    allCohorts(30:33,:)=[];
    allCohorts(36:39,:)=[];
    save('allOptoDeadCohortData.mat','allCohorts','cohort1','cohort2','cohort3','cohort4','cohort5','cohort6');
    cd('O:\sjk\Figures\MGB IC Opto')
    clear cohort1 cohort2 cohort3 cohort4 cohort5 optomeanMat
end