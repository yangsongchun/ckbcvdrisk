********************************************************************************
* Programme: The CKB-CVD models (pragmatic and flexibly updated 10-year 
*            cardiovascular risk prediction models)
* Programmer: Songchun Yang (Peking University Health Science Center)
* Created date: 2022-01-21
* Version: 1.0
* Citation: Yang S, Han Y, Yu C, et al. Development of a model to predict 10-year 
*           risk of ischemic and hemorrhagic stroke and ischemic heart disease 
*           using the China Kadoorie Biobank. Neurology 2022; 98(23).
********************************************************************************

capture program drop ckbcvdrisk
program define ckbcvdrisk, rclass
    syntax [, wkihd(real 1) wbihd(real 0) wkis(real 1) wbis(real 0) wkhs(real 1) wbhs(real 0) ///
            mkihd(real 1) mbihd(real 0) mkis(real 1) mbis(real 0) mkhs(real 1) mbhs(real 0)]
    tempvar lp_ihd lp_is lp_hs age sbp dbp wai

    ** Step 1. Centralization
    gen `age' = age - 55
    gen `sbp' = sbp - 120
    gen `dbp' = dbp - 80
    gen `wai' = wai - 80

    ** Step 2. Log hazard (Linear prediction)
    * Ischemic heart disease
    qui {
    gen `lp_ihd' = 0.0666511 * `age' + ///
                   0.0022094 * `sbp' + ///
                   0.0052584 * `dbp' + ///
                   0.4381349 * hpt + ///
                   0.1601188 * sms + ///
                   0.4633056 * dia + ///
                   0.0146348 * `wai' + ///
                   -0.0000259 * `age' * `sbp' + ///
                   -0.0001553 * `age' * `dbp' + ///
                   -0.0101449 * `age' * hpt + ///
                   -0.0067592 * `age' * sms + ///
                   -0.0115995 * `age' * dia + ///
                   -0.0006227 * `age' * `wai' if sex == 0

    replace `lp_ihd' = 0.0716554 * `age' + ///
                       0.0052816 * `sbp' + ///
                       0.0091778 * `dbp' + ///
                       0.4660422 * hpt + ///
                       0.2024671 * sms + ///
                       0.5971099 * dia + ///
                       0.0125664 * `wai' + ///
                       0.0000168 * `age' * `sbp' + ///
                       -0.0006826 * `age' * `dbp' + ///
                       -0.0053771 * `age' * hpt + ///
                       -0.0037854 * `age' * sms + ///
                       -0.0183296 * `age' * dia + ///
                       -0.0006205 * `age' * `wai' if sex == 1

    * Ischemic stroke
    gen `lp_is' = 0.0740794 * `age' + ///
                  0.0087970 * `sbp' + ///
                  0.0100736 * `dbp' + ///
                  0.4315641 * hpt + ///
                  0.2590549 * sms + ///
                  0.5991041 * dia + ///
                  0.0100013 * `wai' + ///
                  -0.0004318 * `age' * `sbp' + ///
                  -0.0001013 * `age' * `dbp' + ///
                  -0.0089527 * `age' * hpt + ///
                  -0.0155619 * `age' * sms + ///
                  -0.0208146 * `age' * dia + ///
                  -0.0004689 * `age' * `wai' if sex == 0

    replace `lp_is' = 0.0835158 * `age' + ///
                      0.0124603 * `sbp' + ///
                      0.0140617 * `dbp' + ///
                      0.4656016 * hpt + ///
                      0.2512960 * sms + ///
                      0.6059128 * dia + ///
                      0.0077610 * `wai' + ///
                      -0.0004323 * `age' * `sbp' + ///
                      -0.0003044 * `age' * `dbp' + ///
                      -0.0116475 * `age' * hpt + ///
                      -0.0050969 * `age' * sms + ///
                      -0.0127534 * `age' * dia + ///
                      -0.0001151 * `age' * `wai' if sex == 1

    * Hemorrhagic stroke
    gen `lp_hs' = 0.0820293 * `age' + ///
                  0.0147824 * `sbp' + ///
                  0.0280161 * `dbp' + ///
                  0.4797150 * hpt + ///
                  0.2458049 * sms + ///
                  0.4047413 * dia + ///
                  -0.0084531 * `wai' + ///
                  -0.0003590 * `age' * `sbp' + ///
                  -0.0007176 * `age' * `dbp' + ///
                  -0.0089115 * `age' * hpt + ///
                  -0.0062883 * `age' * sms + ///
                  -0.0094574 * `age' * dia + ///
                  -0.0001766 * `age' * `wai' if sex == 0

    replace `lp_hs' = 0.0820196 * `age' + ///
                      0.0168301 * `sbp' + ///
                      0.0287048 * `dbp' + ///
                      0.3949390 * hpt + ///
                      0.0962472 * sms + ///
                      0.3715960 * dia + ///
                      -0.0081349 * `wai' + ///
                      -0.0002353 * `age' * `sbp' + ///
                      -0.0009919 * `age' * `dbp' + ///
                      -0.0100336 * `age' * hpt + ///
                      -0.0009236 * `age' * sms + ///
                      -0.0078866 * `age' * dia + ///
                      -0.0004154 * `age' * `wai' if sex == 1


    ** Step 3. 10-year risk before recalibration
    cap drop ckb_ihd_risk ckb_is_risk ckb_hs_risk ckb_cvd_risk
    gen ckb_ihd_risk = 1 - 0.8996251^exp(`lp_ihd') if sex == 0
        replace ckb_ihd_risk = 1 - 0.9287638^exp(`lp_ihd') if sex == 1
    gen ckb_is_risk = 1 - 0.9103525^exp(`lp_is') if sex == 0
        replace ckb_is_risk = 1 - 0.9275716^exp(`lp_is') if sex == 1
    gen ckb_hs_risk = 1 - 0.9858385^exp(`lp_hs') if sex == 0
        replace ckb_hs_risk = 1 - 0.9825860^exp(`lp_hs') if sex == 1
    }

    gen ckb_cvd_risk = 1 - (1-ckb_ihd_risk)*(1-ckb_is_risk)*(1-ckb_hs_risk)        
    di as txt "The original 10-year predicted risk of ischemic heart disease is generated as 'ckb_ihd_risk'."
    di as txt "The original 10-year predicted risk of ischemic stroke is generated as 'ckb_is_risk'."
    di as txt "The original 10-year predicted risk of hemorrhagic stroke is generated as 'ckb_hs_risk'."
    di as txt "The original 10-year predicted risk of total cardiovascular disease is generated as 'ckb_cvd_risk'."

    ** Step 4. 10-year risk after recalibration
    cap drop ckb_ihd_recalrisk ckb_is_recalrisk ckb_hs_recalrisk ckb_cvd_recalrisk
    qui {
    gen ckb_ihd_recalrisk = 1 - exp(-exp(`wbihd' + `wkihd' * ln(-ln(1 - ckb_ihd_risk)))) if sex == 0
        replace ckb_ihd_recalrisk = 1 - exp(-exp(`mbihd' + `mkihd' * ln(-ln(1 - ckb_ihd_risk)))) if sex == 1
    gen ckb_hs_recalrisk = 1 - exp(-exp(`wbhs' + `wkhs' * ln(-ln(1 - ckb_hs_risk)))) if sex == 0
        replace ckb_hs_recalrisk = 1 - exp(-exp(`mbhs' + `mkhs' * ln(-ln(1 - ckb_hs_risk)))) if sex == 1
    gen ckb_is_recalrisk = 1 - exp(-exp(`wbis' + `wkis' * ln(-ln(1 - ckb_is_risk)))) if sex == 0
        replace ckb_is_recalrisk = 1 - exp(-exp(`mbis' + `mkis' * ln(-ln(1 - ckb_is_risk)))) if sex == 1
    }

    gen ckb_cvd_recalrisk = 1 - (1-ckb_ihd_recalrisk)*(1-ckb_is_recalrisk)*(1-ckb_hs_recalrisk)
    di as txt "The recalibrated 10-year predicted risk of ischemic heart disease is generated as 'ckb_ihd_recalrisk'."
    di as txt "The recalibrated 10-year predicted risk of ischemic stroke is generated as 'ckb_is_recalrisk'."
    di as txt "The recalibrated 10-year predicted risk of hemorrhagic stroke is generated as 'ckb_hs_recalrisk'."
    di as txt "The recalibrated 10-year predicted risk of total cardiovascular disease is generated as 'ckb_cvd_recalrisk'."
end
