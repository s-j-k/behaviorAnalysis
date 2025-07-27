function allCohorts = loadAllOptoCohorts(cohortRange)
% takes input of an array which is the cohort number, and loads data from
% those cohorts specified
    for ii=1:length(cohortRange)
        if cohortRange(ii)==1
                cd('O:\sjk\Behavior\OptoMGBIC_1')
                load('summaryData.mat');
                cohort1=optomeanMat;
        elseif cohortRange(ii)==2
                cd('O:\sjk\Behavior\OptoMGBIC_2')
                load('summaryData.mat');
                cohort2=optomeanMat;
        elseif cohortRange(ii)==3
                cd('O:\sjk\Behavior\OptoMGBIC_3')
                load('summaryData.mat');
                cohort3=optomeanMat;
        elseif cohortRange(ii)==4
                cd('O:\sjk\Behavior\MGBIC_4')
                load('summaryData.mat');
                cohort4=optomeanMat;
        elseif cohortRange(ii)==5
                cd('O:\sjk\Behavior\MGBIC_5')
                load('summaryData.mat');
                cohort5=optomeanMat;
        elseif cohortRange(ii)==6
                cd('O:\sjk\Behavior\MGBIC_6')
                load('summaryData.mat');
                cohort6=optomeanMat;
        elseif cohortRange(ii)==7
                cd('O:\sjk\Behavior\cohort_7')
                load('summaryData.mat');
                cohort7=optomeanMat;
        elseif cohortRange(ii)==8
                cd('O:\sjk\Behavior\cohort_8')
                load('summaryData.mat');
                cohort8=optomeanMat;        
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
    cohort7(1,:)=[];
    allCohorts=vertcat(cohort1, cohort2, cohort3, cohort4, cohort5, cohort6, cohort7, cohort8);
    allCohorts(26:27,:)=[]; 
    allCohorts(28:31,:)=[]; 
    allCohorts(30:33,:)=[];
    allCohorts(36:41,:)=[]; % keep 189-191, 195, 199-205
    allCohorts(38:41,:)=[]; %keep 198, 203-205
    allCohorts(40:47,:)=[];
    allCohorts(44:45,:)=[];
    allCohorts(60:68,:)=[];
%     allCohorts(66:67,:)=[];
    save('allOptoCohortData.mat','allCohorts','cohort1','cohort2','cohort3','cohort4','cohort5','cohort6','cohort7','cohort8');
    cd('O:\sjk\Figures\MGB IC Opto')
    clear cohort1 cohort2 cohort3 cohort4 cohort5 cohort6 cohort7 cohort8 optomeanMat
end