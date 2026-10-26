---
title: "Cost-Sensitive Ordinal Classification Methods to Predict SARS-CoV-2 Pneumonia Severity"
authors:
  - "Fernando García-García"
  - admin
  - "Pedro Pablo España Yandiola"
  - "Isabel Urrutia Landa"
  - "Joaquín Martínez-Minaya"
  - "Miren Hayet-Otero"
  - "Mónica Nieves Ermecheo"
  - "José María Quintana"
  - "Rosario Menéndez"
  - "Antoni Torres"
  - "Rafael Zalacain Jorge"
date: "2024-05-01"
publishDate: "2024-05-01"
doi: "10.1109/jbhi.2024.3363765"
publication_types: ["2"]
publication: "*IEEE Journal of Biomedical and Health Informatics*"
publication_short: ""
abstract: "OBJECTIVE: To study the suitability of costsensitive ordinal artificial intelligence-machine learning (AIML) strategies in the prognosis of SARS-CoV-2 pneumonia severity. MATERIALS & METHODS: Observational, retrospective, longitudinal, cohort study in 4 hospitals in Spain. Information regarding demographic and clinical status was supplemented by socioeconomic data and air pollution exposures. We proposed AI-ML algorithms for ordinal classification via ordinal decomposition and for cost-sensitive learning via resampling techniques. For performance-based model selection, we defined a custom score including per-class sensitivities and asymmetric misprognosis costs. 260 distinct AI-ML models were evaluated via 10 repetitions of 5×5 nested cross-validation with hyperparameter tuning. Model selection was followed by the calibration of predicted probabilities. Final overall performance was compared against five well-established clinical severity scores and against a 'standard' (non-cost sensitive, non-ordinal) AI-ML baseline. In our best model, we also evaluated its explainability with respect to each of the input variables. RESULTS: The study enrolled n = 1548 patients: 712 experienced low, 238 medium, and 598 high clinical severity. d = 131 variables were collected, becoming d ' = 148 features after categorical encoding. Model selection resulted in our best-performing AI-ML pipeline having: a) no imputation of missing data, b) no feature selection (i.e. using the full set of d ' features), c) 'Ordered Partitions' ordinal decomposition, d) cost-based reimbalance, and e) a Histogram-based Gradient Boosting classifier. This best model (calibrated) obtained a median accuracy of 68.1% [67.3%, 68.8%] (95% confidence interval), a balanced accuracy of 57.0% [55.6%, 57.9%], and an overall area under the curve (AUC) 0.802 [0.795, 0.808]. In our dataset, it outperformed all five clinical severity scores and the 'standard' AI-ML baseline. DISCUSSION & CONCLUSION: We conducted an exhaustive exploration of AI-ML methods designed for both ordinal and cost-sensitive classification, motivated by a real-world application domain (clinical severity prognosis) in which these topics arise naturally. Our model with the best classification performance exploited successfully the ordering information of ground truth classes, coping with imbalance and asymmetric costs. However, these ordinal and cost-sensitive aspects are seldom explored in the literature."
featured: false
url_pdf: "https://bird.bcamath.org/bitstream/20.500.11824/1790/1/manuscript-repo.pdf"
url_code: ""
url_dataset: ""
url_poster: ""
url_project: ""
url_slides: ""
url_source: ""
url_video: ""
orcid_put_code: 159080471
---
