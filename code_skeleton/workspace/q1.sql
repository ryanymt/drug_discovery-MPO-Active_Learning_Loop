SELECT 
    R.run_id, 
    S.molecule_hash, 
    S.final_score, 
    S.qed_score, 
    S.mol_wt,
    S.logp,
    S.sa_score,
    S.toxicity_label, 
    S.score_method
FROM `gcda-apac-sc.bioops_platform.screening_results` S
JOIN `gcda-apac-sc.bioops_platform.molecule_registry` R ON S.molecule_hash = R.molecule_hash
WHERE R.run_id = 'cycle2_10k' AND S.final_score IS NOT NULL
