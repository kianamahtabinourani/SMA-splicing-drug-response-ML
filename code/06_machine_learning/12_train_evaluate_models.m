%% 12 Train Evaluate Models
% Refactored from the original research script.
% Numeric methods, thresholds, and analysis statements retained.
% Only file paths, filename, and explanatory comments were adjusted.

%% %% Advanced ML Regression with Corrected Bayesian Optimization
clc; clear; close all;
%% Portable project paths (no scientific calculations changed)
projectRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
% Place external input data under data/raw; see README.md.


%% ===================== SETTINGS & PATHS =====================
outDirML = fullfile(projectRoot, "Machine learning_Optimized");
if ~exist(outDirML,"dir"), mkdir(outDirML); end

% Use the responsive/non-responsive tables
RG_SMN2_resp = fullfile(projectRoot, "Sequence functional analysis", "RG_SMN2_responsive_sequence_functional.xlsx");
RG_SMN2_nonr = fullfile(projectRoot, "Sequence functional analysis", "RG_SMN2_nonresponsive_sequence_functional.xlsx");
NVS_SMN2_resp = fullfile(projectRoot, "Sequence functional analysis", "NVS_SMN2_responsive_sequence_functional.xlsx");
NVS_SMN2_nonr = fullfile(projectRoot, "Sequence functional analysis", "NVS_SMN2_nonresponsive_sequence_functional.xlsx");

targetsToRun = ["RG","NVS"];
rng(1); 

featWanted = ["ESE_density","ESS_density","ESE_minus_ESS","ESE_to_ESS_ratio", ...
              "sf2_max","sc35_max","srp40_max","srp55_max", ...
              "sf2_hits","sc35_hits","srp40_hits","srp55_hits", ...
              "donor_mismatch_count","donor10_GC","exon7_GC", ...
              "upIntron200_GC","downIntron200_GC","upIntron_last30_pyFrac", ...
              "PSI_DMSO_mean"];

allResults = table();

%% ===================== MAIN LOOP =====================
for t = 1:numel(targetsToRun)
    drug = targetsToRun(t);
    fprintf("\n--- Processing %s ---\n", drug);

    if drug == "RG"
        T = [readtable(RG_SMN2_resp); readtable(RG_SMN2_nonr)];
        yName = "dPSI_RG_mean";
    else
        T = [readtable(NVS_SMN2_resp); readtable(NVS_SMN2_nonr)];
        yName = "dPSI_NVS_mean";
    end

    % 1. Feature Check (Verify names match exactly)
    feat = featWanted(ismember(featWanted, string(T.Properties.VariableNames)));
    fprintf("  Using %d out of %d requested features.\n", numel(feat), numel(featWanted));
    
    X = table2array(T(:, feat));
    y = T.(yName);

    % 2. Median Imputation (Fill NaNs to keep dataset size)
    for c = 1:size(X,2)
        col = X(:,c);
        col(isnan(col)) = median(col, 'omitnan');
        X(:,c) = col;
    end
    goodY = ~isnan(y);
    X = X(goodY,:); y = y(goodY);
    fprintf("  Final dataset size: %d samples.\n", numel(y));

    % 3. Corrected Optimization Options
    % Changed 'bayesian' to 'bayesopt'
    opt = struct('Optimizer','bayesopt', 'MaxObjectiveEvaluations', 20, 'ShowPlots', false);

    %% --- 1. Optimized Random Forest ---
    fprintf("  Optimizing Random Forest...\n");
    rfMdl = fitrensemble(X, y, 'Method', 'Bag', ...
        'OptimizeHyperparameters', {'NumLearningCycles', 'MinLeafSize'}, ...
        'HyperparameterOptimizationOptions', opt);
    yhat_RF = kfoldPredict(crossval(rfMdl, 'KFold', 5));

    %% --- 2. Optimized GBM (LSBoost) ---
    fprintf("  Optimizing Gradient Boosting...\n");
    gbmMdl = fitrensemble(X, y, 'Method', 'LSBoost', ...
        'OptimizeHyperparameters', {'NumLearningCycles', 'LearnRate', 'MinLeafSize'}, ...
        'HyperparameterOptimizationOptions', opt);
    yhat_GBM = kfoldPredict(crossval(gbmMdl, 'KFold', 5));

    %% --- 3. Optimized SVM ---
    fprintf("  Optimizing SVM...\n");
    svmMdl = fitrsvm(X, y, 'Standardize', true, ...
        'OptimizeHyperparameters', {'BoxConstraint', 'KernelScale'}, ...
        'HyperparameterOptimizationOptions', opt);
    yhat_SVM = kfoldPredict(crossval(svmMdl, 'KFold', 5));

    %% --- 4. Deep Neural Network (Stable for small data) ---
    fprintf("  Training DNN...\n");
    % Small layers [10 5] + strong Lambda (0.1) prevents mathematical explosion
    dnnMdl = fitrnet(X, y, 'LayerSizes', [10 5], ...
        'Standardize', true, 'Lambda', 0.1, 'IterationLimit', 1000);
    yhat_DNN = kfoldPredict(crossval(dnnMdl, 'KFold', 5));

    %% ===================== Metrics & Results =====================
    resRF  = calculateMetrics(drug, "RandomForest", y, yhat_RF);
    resGBM = calculateMetrics(drug, "GradientBoosting", y, yhat_GBM);
    resSVM = calculateMetrics(drug, "SVM", y, yhat_SVM);
    resDNN = calculateMetrics(drug, "DNN", y, yhat_DNN);
    
    allResults = [allResults; resRF; resGBM; resSVM; resDNN];
    disp([resRF; resGBM; resSVM; resDNN]);
    
    % Save comparison plot
    fig = figure('Position', [100 100 1200 350], 'Visible', 'off');
    tiledlayout(1,4,'Padding','compact');
    mdls = {"RF", "GBM", "SVM", "DNN"}; preds = {yhat_RF, yhat_GBM, yhat_SVM, yhat_DNN};
    for i = 1:4
        nexttile;
        scatter(y, preds{i}, 15, 'filled', 'MarkerFaceAlpha', 0.5); hold on;
        plot([min(y) max(y)], [min(y) max(y)], 'r--');
        title(sprintf('%s (R^2=%.2f)', mdls{i}, calculateR2(y, preds{i})));
        xlabel('Measured'); ylabel('Predicted'); grid on;
    end
    saveas(fig, fullfile(outDirML, drug + "_Comparison_Optimized.png"));
end

writetable(allResults, fullfile(outDirML, "Performance_Summary_Optimized.xlsx"));
fprintf("\nDONE. Results saved in: %s\n", outDirML);

%% ===================== FUNCTIONS =====================
function out = calculateMetrics(drug, modelName, y, yhat)
    mse = mean((y - yhat).^2);
    R2 = 1 - (sum((y - yhat).^2) / sum((y - mean(y)).^2));
    rho = corr(y, yhat, 'Type', 'Spearman');
    out = table(string(drug), string(modelName), R2, sqrt(mse), mse, rho, ...
        'VariableNames', {'Drug','Model','R2','RMSE','MSE','SpearmanR'});
end

function r2 = calculateR2(y, yh)
    r2 = 1 - (sum((y - yh).^2) / sum((y - mean(y)).^2));
end


-----------------------------------------------------------------------------

--- Processing RG ---
  Using 19 out of 19 requested features.
  Final dataset size: 285 samples.
  Optimizing Random Forest...
|=====================================================================================================|
| Iter | Eval   | Objective:  | Objective   | BestSoFar   | BestSoFar   | NumLearningC-|  MinLeafSize |
|      | result | log(1+loss) | runtime     | (observed)  | (estim.)    | ycles        |              |
|=====================================================================================================|
|    1 | Best   |      3.7932 |      1.1039 |      3.7932 |      3.7932 |           94 |           50 |
|    2 | Best   |      3.6702 |      1.8908 |      3.6702 |      3.6787 |          203 |            5 |
|    3 | Accept |      3.9013 |      0.1548 |      3.6702 |      3.6705 |           15 |           70 |
|    4 | Best   |      3.6187 |     0.19558 |      3.6187 |       3.651 |           16 |            3 |
|    5 | Accept |      3.6697 |     0.17587 |      3.6187 |      3.6187 |           10 |            2 |
|    6 | Accept |      3.6896 |     0.12411 |      3.6187 |      3.6576 |           10 |            5 |
|    7 | Accept |      3.6926 |     0.12656 |      3.6187 |      3.6187 |           10 |            1 |
|    8 | Accept |      3.6598 |     0.34573 |      3.6187 |      3.6522 |           34 |            2 |
|    9 | Accept |      3.7155 |     0.21775 |      3.6187 |      3.6712 |           21 |            3 |
|   10 | Accept |      3.6842 |     0.10896 |      3.6187 |      3.6699 |           10 |           16 |
|   11 | Accept |      3.6449 |     0.12313 |      3.6187 |      3.6685 |           10 |            9 |
|   12 | Accept |      3.7069 |     0.12811 |      3.6187 |      3.6708 |           10 |            7 |
|   13 | Accept |      3.6903 |     0.11922 |      3.6187 |      3.6706 |           10 |           11 |
|   14 | Accept |      3.7098 |     0.10774 |      3.6187 |      3.6709 |           10 |           26 |
|   15 | Accept |      3.6578 |     0.11767 |      3.6187 |      3.6679 |           10 |            2 |
|   16 | Accept |      4.0263 |      0.0989 |      3.6187 |      3.6676 |           10 |          142 |
|   17 | Accept |      3.6915 |     0.12154 |      3.6187 |      3.6719 |           10 |            2 |
|   18 | Accept |      3.7219 |     0.12659 |      3.6187 |      3.6738 |           10 |            4 |
|   19 | Accept |       3.689 |      0.1534 |      3.6187 |      3.6735 |           10 |            9 |
|   20 | Accept |      3.9949 |     0.13496 |      3.6187 |      3.6724 |           10 |          105 |

__________________________________________________________
Optimization completed.
MaxObjectiveEvaluations of 20 reached.
Total function evaluations: 20
Total elapsed time: 9.7147 seconds
Total objective function evaluation time: 5.6753

Best observed feasible point:
    NumLearningCycles    MinLeafSize
    _________________    ___________

           16                 3     

Observed objective function value = 3.6187
Estimated objective function value = 3.6781
Function evaluation time = 0.19558

Best estimated feasible point (according to models):
    NumLearningCycles    MinLeafSize
    _________________    ___________

           34                 2     

Estimated objective function value = 3.6724
Estimated function evaluation time = 0.36527

  Optimizing Gradient Boosting...
|====================================================================================================================|
| Iter | Eval   | Objective:  | Objective   | BestSoFar   | BestSoFar   | NumLearningC-|    LearnRate |  MinLeafSize |
|      | result | log(1+loss) | runtime     | (observed)  | (estim.)    | ycles        |              |              |
|====================================================================================================================|
|    1 | Best   |      4.0325 |      0.2525 |      4.0325 |      4.0325 |           27 |    0.0012633 |           41 |
|    2 | Best   |      3.7251 |     0.12114 |      3.7251 |      3.7407 |           12 |      0.86858 |           36 |
|    3 | Accept |      3.7277 |      1.5203 |      3.7251 |      3.7251 |          230 |     0.074036 |           41 |
|    4 | Accept |      3.9555 |     0.61302 |      3.7251 |      3.7251 |           79 |      0.10827 |            1 |
|    5 | Accept |      3.7701 |     0.10455 |      3.7251 |      3.7251 |           10 |      0.29136 |           56 |
|    6 | Best   |      3.6664 |     0.10027 |      3.6664 |      3.6664 |           10 |      0.24588 |           27 |
|    7 | Accept |      3.7266 |     0.10511 |      3.6664 |      3.6664 |           10 |      0.31318 |           15 |
|    8 | Accept |      3.6944 |     0.11806 |      3.6664 |      3.6665 |           10 |      0.10614 |           24 |
|    9 | Accept |       3.667 |     0.11547 |      3.6664 |      3.6664 |           10 |      0.37548 |           24 |
|   10 | Accept |       3.682 |     0.13615 |      3.6664 |      3.6724 |           14 |      0.29969 |           26 |
|   11 | Best   |      3.6644 |     0.11688 |      3.6644 |      3.6647 |           10 |      0.27163 |           26 |
|   12 | Accept |       4.049 |     0.10515 |      3.6644 |      3.6701 |           10 |    0.0010098 |            1 |
|   13 | Best   |      3.6635 |     0.10529 |      3.6635 |      3.6643 |           10 |      0.25396 |           27 |
|   14 | Accept |      4.1309 |     0.11309 |      3.6635 |      3.6642 |           10 |       0.8278 |            3 |
|   15 | Accept |      3.6862 |     0.10062 |      3.6635 |      3.6705 |           10 |      0.35447 |           27 |
|   16 | Best   |      3.6605 |     0.12231 |      3.6605 |      3.6689 |           10 |      0.24743 |           24 |
|   17 | Accept |      3.6605 |     0.10611 |      3.6605 |      3.6602 |           10 |      0.24408 |           24 |
|   18 | Accept |      3.8518 |     0.10855 |      3.6605 |      3.6602 |           10 |      0.79948 |           20 |
|   19 | Accept |      3.6886 |     0.10258 |      3.6605 |      3.6602 |           10 |      0.11275 |           38 |
|   20 | Accept |      3.8228 |     0.10674 |      3.6605 |      3.6602 |           10 |     0.044454 |           33 |

__________________________________________________________
Optimization completed.
MaxObjectiveEvaluations of 20 reached.
Total function evaluations: 20
Total elapsed time: 6.0438 seconds
Total objective function evaluation time: 4.2739

Best observed feasible point:
    NumLearningCycles    LearnRate    MinLeafSize
    _________________    _________    ___________

           10             0.24743         24     

Observed objective function value = 3.6605
Estimated objective function value = 3.6602
Function evaluation time = 0.12231

Best estimated feasible point (according to models):
    NumLearningCycles    LearnRate    MinLeafSize
    _________________    _________    ___________

           10             0.24743         24     

Estimated objective function value = 3.6602
Estimated function evaluation time = 0.10824

  Optimizing SVM...
|=====================================================================================================|
| Iter | Eval   | Objective:  | Objective   | BestSoFar   | BestSoFar   | BoxConstraint|  KernelScale |
|      | result | log(1+loss) | runtime     | (observed)  | (estim.)    |              |              |
|=====================================================================================================|
|    1 | Best   |      9.7271 |      9.3408 |      9.7271 |      9.7271 |       255.42 |    0.0025801 |
|    2 | Best   |       7.861 |    0.054595 |       7.861 |       7.963 |     0.020714 |      0.11179 |
|    3 | Best   |       4.071 |    0.028215 |       4.071 |      4.0713 |       218.35 |          678 |
|    4 | Best   |      4.0523 |     0.02949 |      4.0523 |      4.0551 |    0.0011097 |       296.97 |
|    5 | Accept |      4.3612 |      0.0315 |      4.0523 |      4.0556 |       415.87 |       364.67 |
|    6 | Best   |      4.0522 |    0.023165 |      4.0522 |      4.0246 |    0.0010027 |       992.34 |
|    7 | Accept |      4.0522 |    0.024043 |      4.0522 |      4.0468 |    0.0012551 |       999.98 |
|    8 | Accept |      4.0522 |    0.024504 |      4.0522 |      4.0503 |       1.2262 |       997.53 |
|    9 | Accept |      4.0529 |      2.3392 |      4.0522 |      4.0507 |    0.0010007 |       9.1342 |
|   10 | Best   |      4.0521 |     0.02536 |      4.0521 |       4.051 |    0.0010026 |        571.7 |
|   11 | Accept |      4.1012 |    0.023798 |      4.0521 |      4.0472 |       934.77 |       996.57 |
|   12 | Best   |       4.052 |    0.023858 |       4.052 |       4.043 |        0.001 |       86.013 |
|   13 | Accept |      4.0526 |    0.026034 |       4.052 |      4.0435 |    0.0010118 |       39.635 |
|   14 | Accept |      4.0526 |    0.023321 |       4.052 |      4.0422 |    0.0010003 |       28.526 |
|   15 | Accept |      4.0526 |    0.024439 |       4.052 |      4.0506 |    0.0010117 |       22.831 |
|   16 | Accept |      4.0532 |    0.023614 |       4.052 |      4.0528 |       35.883 |       999.76 |
|   17 | Accept |      4.0521 |     0.02339 |       4.052 |       4.049 |      0.41046 |       771.08 |
|   18 | Best   |      4.0518 |    0.025587 |      4.0518 |      4.0526 |    0.0010052 |       157.26 |
|   19 | Accept |      4.0524 |     0.02491 |      4.0518 |      4.0462 |        6.096 |       834.59 |
|   20 | Accept |      4.0527 |    0.025055 |      4.0518 |       4.046 |    0.0010047 |        17.59 |

__________________________________________________________
Optimization completed.
MaxObjectiveEvaluations of 20 reached.
Total function evaluations: 20
Total elapsed time: 13.5165 seconds
Total objective function evaluation time: 12.1649

Best observed feasible point:
    BoxConstraint    KernelScale
    _____________    ___________

      0.0010052        157.26   

Observed objective function value = 4.0518
Estimated objective function value = 4.0523
Function evaluation time = 0.025587

Best estimated feasible point (according to models):
    BoxConstraint    KernelScale
    _____________    ___________

        6.096          834.59   

Estimated objective function value = 4.046
Estimated function evaluation time = 0.024176

  Training DNN...
    Drug          Model              R2        RMSE      MSE      SpearmanR
    ____    __________________    ________    ______    ______    _________

    "RG"    "RandomForest"          0.2729    6.3166    39.899      0.42751
    "RG"    "GradientBoosting"     0.24905    6.4194    41.208      0.32888
    "RG"    "SVM"                 -0.03047    7.5197    56.547    0.0022182
    "RG"    "DNN"                  0.25223    6.4057    41.033      0.42701


--- Processing NVS ---
  Using 19 out of 19 requested features.
  Final dataset size: 285 samples.
  Optimizing Random Forest...
|=====================================================================================================|
| Iter | Eval   | Objective:  | Objective   | BestSoFar   | BestSoFar   | NumLearningC-|  MinLeafSize |
|      | result | log(1+loss) | runtime     | (observed)  | (estim.)    | ycles        |              |
|=====================================================================================================|
|    1 | Best   |      4.4992 |      3.0107 |      4.4992 |      4.4992 |          346 |           19 |
|    2 | Best   |      4.3366 |     0.17474 |      4.3366 |      4.3451 |           13 |            1 |
|    3 | Accept |      4.5833 |     0.41795 |      4.3366 |      4.3366 |           44 |           52 |
|    4 | Accept |      4.3581 |     0.80806 |      4.3366 |      4.3366 |           77 |            6 |
|    5 | Accept |      4.3791 |     0.12606 |      4.3366 |      4.3595 |           10 |            2 |
|    6 | Accept |      4.4112 |     0.11994 |      4.3366 |      4.3702 |           10 |            1 |
|    7 | Best   |      4.2893 |     0.12699 |      4.2893 |      4.3389 |           10 |            3 |
|    8 | Accept |      4.4184 |     0.11555 |      4.2893 |      4.3662 |           10 |            3 |
|    9 | Accept |      4.6947 |      0.1076 |      4.2893 |      4.3623 |           10 |          142 |
|   10 | Accept |      4.4124 |       0.118 |      4.2893 |      4.3674 |           10 |            4 |
|   11 | Accept |      4.3118 |     0.11878 |      4.2893 |        4.36 |           10 |            1 |
|   12 | Best   |      4.2517 |     0.12062 |      4.2517 |      4.3332 |           10 |            1 |
|   13 | Accept |      4.4887 |     0.12531 |      4.2517 |      4.3315 |           10 |           10 |
|   14 | Accept |      4.5235 |     0.12581 |      4.2517 |      4.3311 |           10 |           31 |
|   15 | Accept |      4.5812 |     0.12157 |      4.2517 |       4.331 |           10 |           91 |
|   16 | Accept |       4.509 |     0.12564 |      4.2517 |      4.3301 |           10 |           14 |
|   17 | Accept |      4.3746 |       0.151 |      4.2517 |      4.3379 |           10 |            1 |
|   18 | Accept |      4.6276 |     0.11174 |      4.2517 |      4.3377 |           10 |           69 |
|   19 | Accept |      4.4255 |     0.12343 |      4.2517 |      4.3375 |           10 |            7 |
|   20 | Accept |      4.5519 |      0.1177 |      4.2517 |      4.3372 |           10 |           23 |

__________________________________________________________
Optimization completed.
MaxObjectiveEvaluations of 20 reached.
Total function evaluations: 20
Total elapsed time: 9.021 seconds
Total objective function evaluation time: 6.3672

Best observed feasible point:
    NumLearningCycles    MinLeafSize
    _________________    ___________

           10                 1     

Observed objective function value = 4.2517
Estimated objective function value = 4.3372
Function evaluation time = 0.12062

Best estimated feasible point (according to models):
    NumLearningCycles    MinLeafSize
    _________________    ___________

           13                 1     

Estimated objective function value = 4.3372
Estimated function evaluation time = 0.16054

  Optimizing Gradient Boosting...
|====================================================================================================================|
| Iter | Eval   | Objective:  | Objective   | BestSoFar   | BestSoFar   | NumLearningC-|    LearnRate |  MinLeafSize |
|      | result | log(1+loss) | runtime     | (observed)  | (estim.)    | ycles        |              |              |
|====================================================================================================================|
|    1 | Best   |      4.3053 |      2.8509 |      4.3053 |      4.3053 |          338 |     0.020048 |           21 |
|    2 | Accept |      4.7957 |      0.2223 |      4.3053 |      4.3322 |           22 |    0.0016567 |           19 |
|    3 | Accept |      4.6844 |      0.3666 |      4.3053 |      4.3054 |           42 |      0.92685 |          128 |
|    4 | Accept |      4.4533 |     0.24633 |      4.3053 |      4.3309 |           24 |      0.43275 |            6 |
|    5 | Accept |      4.6844 |     0.11672 |      4.3053 |      4.3054 |           10 |      0.69474 |          142 |
|    6 | Accept |      4.6059 |     0.12799 |      4.3053 |      4.3054 |           10 |     0.039101 |           17 |
|    7 | Accept |      4.7258 |     0.13882 |      4.3053 |      4.3054 |           10 |      0.01563 |           23 |
|    8 | Best   |      4.2005 |     0.13827 |      4.2005 |      4.2062 |           10 |      0.41099 |            1 |
|    9 | Accept |       4.544 |     0.12649 |      4.2005 |      4.2005 |           10 |      0.99394 |            1 |
|   10 | Best   |      4.1596 |     0.11012 |      4.1596 |      4.1596 |           10 |       0.3196 |            1 |
|   11 | Best   |      4.1422 |     0.11019 |      4.1422 |      4.1423 |           10 |      0.18447 |            1 |
|   12 | Accept |      4.1833 |     0.16372 |      4.1422 |      4.1424 |           18 |      0.20183 |            1 |
|   13 | Accept |      4.2658 |     0.10753 |      4.1422 |      4.1424 |           10 |     0.088503 |            1 |
|   14 | Accept |      4.2242 |     0.12497 |      4.1422 |      4.1425 |           10 |      0.19598 |            3 |
|   15 | Best   |      4.1022 |     0.12545 |      4.1022 |      4.1217 |           10 |       0.2245 |            1 |
|   16 | Accept |      4.1745 |     0.12414 |      4.1022 |      4.1375 |           10 |      0.23044 |            1 |
|   17 | Accept |      4.8078 |     0.12456 |      4.1022 |      4.1376 |           10 |    0.0010433 |            1 |
|   18 | Accept |      4.7623 |     0.12942 |      4.1022 |      4.1376 |           10 |    0.0062113 |            1 |
|   19 | Accept |      4.8147 |     0.12714 |      4.1022 |      4.1378 |           10 |    0.0010012 |          136 |
|   20 | Accept |      4.6921 |     0.12138 |      4.1022 |      4.1373 |           10 |      0.13048 |          142 |

__________________________________________________________
Optimization completed.
MaxObjectiveEvaluations of 20 reached.
Total function evaluations: 20
Total elapsed time: 7.9466 seconds
Total objective function evaluation time: 5.703

Best observed feasible point:
    NumLearningCycles    LearnRate    MinLeafSize
    _________________    _________    ___________

           10             0.2245           1     

Observed objective function value = 4.1022
Estimated objective function value = 4.1373
Function evaluation time = 0.12545

Best estimated feasible point (according to models):
    NumLearningCycles    LearnRate    MinLeafSize
    _________________    _________    ___________

           10             0.2245           1     

Estimated objective function value = 4.1373
Estimated function evaluation time = 0.12305

  Optimizing SVM...
|=====================================================================================================|
| Iter | Eval   | Objective:  | Objective   | BestSoFar   | BestSoFar   | BoxConstraint|  KernelScale |
|      | result | log(1+loss) | runtime     | (observed)  | (estim.)    |              |              |
|=====================================================================================================|
|    1 | Best   |      4.8007 |    0.029823 |      4.8007 |      4.8007 |       13.867 |       522.16 |
|    2 | Accept |      8.6637 |         9.6 |      4.8007 |      4.9952 |      0.29815 |    0.0036538 |
|    3 | Accept |       4.801 |    0.028903 |      4.8007 |      4.8006 |    0.0036207 |       508.01 |
|    4 | Accept |      6.7643 |      9.4645 |      4.8007 |      4.8008 |       480.74 |    0.0048022 |
|    5 | Accept |      4.8011 |    0.024092 |      4.8007 |      4.8008 |     0.029622 |       833.95 |
|    6 | Best   |      4.7969 |     0.02505 |      4.7969 |       4.797 |        992.7 |       991.27 |
|    7 | Accept |       4.801 |    0.023494 |      4.7969 |       4.797 |    0.0010109 |        988.3 |
|    8 | Best   |      4.7556 |    0.024278 |      4.7556 |      4.7554 |       993.74 |       154.69 |
|    9 | Best   |      4.7354 |    0.035563 |      4.7354 |      4.7361 |       997.55 |       20.672 |
|   10 | Best   |       4.728 |    0.026029 |       4.728 |      4.7301 |       998.56 |       48.664 |
|   11 | Accept |      4.7876 |    0.024581 |       4.728 |      4.7294 |        12.53 |       58.317 |
|   12 | Accept |      4.7294 |    0.028737 |       4.728 |      4.7288 |       999.63 |       42.429 |
|   13 | Accept |      4.7282 |    0.028736 |       4.728 |      4.7286 |       995.22 |        40.21 |
|   14 | Accept |      4.7292 |    0.027417 |       4.728 |      4.7286 |          997 |       38.121 |
|   15 | Accept |      4.8011 |    0.024744 |       4.728 |      4.7288 |    0.0010023 |       25.212 |
|   16 | Accept |      4.8011 |    0.023438 |       4.728 |       4.729 |    0.0010038 |       87.438 |
|   17 | Accept |      6.7693 |     0.19137 |       4.728 |      4.7274 |        999.7 |       4.4543 |
|   18 | Accept |      4.8008 |    0.023912 |       4.728 |      4.7271 |    0.0010072 |       195.26 |
|   19 | Accept |      4.7295 |    0.034704 |       4.728 |      4.7294 |        995.7 |       27.252 |
|   20 | Accept |      4.7311 |    0.030409 |       4.728 |      4.7294 |       992.21 |       84.959 |

__________________________________________________________
Optimization completed.
MaxObjectiveEvaluations of 20 reached.
Total function evaluations: 20
Total elapsed time: 21.0954 seconds
Total objective function evaluation time: 19.7198

Best observed feasible point:
    BoxConstraint    KernelScale
    _____________    ___________

       998.56          48.664   

Observed objective function value = 4.728
Estimated objective function value = 4.7295
Function evaluation time = 0.026029

Best estimated feasible point (according to models):
    BoxConstraint    KernelScale
    _____________    ___________

       995.22           40.21   

Estimated objective function value = 4.7294
Estimated function evaluation time = 0.028216

  Training DNN...
    Drug           Model              R2         RMSE      MSE      SpearmanR
    _____    __________________    _________    ______    ______    _________

    "NVS"    "RandomForest"          0.26174    8.8778    78.814     0.53917 
    "NVS"    "GradientBoosting"      0.24989    8.9487     80.08     0.52405 
    "NVS"    "SVM"                 -0.051004    10.593     112.2     0.42379 
    "NVS"    "DNN"                   0.19177    9.2889    86.285     0.41181 


DONE. Results saved in: C:\Users\LENOVO\Desktop\SMA project\Machine learning_Optimized