SELECT 
    R.run_id, 
    S.molecule_hash, 
    AVG(S.final_score) as final_score, 
    AVG(S.qed_score) as qed_score, 
    AVG(S.mol_wt) as avg_molwt,
    AVG(S.logp) as avg_logp,
    AVG(S.sa_score) as avg_sa_score,
    MAX(S.toxicity_label) as toxicity_label, 
    'xgboost' as score_method
FROM `gcda-apac-sc.bioops_platform.screening_results` S
JOIN `gcda-apac-sc.bioops_platform.molecule_registry` R ON S.molecule_hash = R.molecule_hash
WHERE R.run_id = 'prod_100k_v1' AND S.final_score IS NOT NULL
GROUP BY 1, 2
