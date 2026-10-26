---
title: "Extracting relevant predictive variables for COVID-19 severity prognosis: An exhaustive comparison of feature selection techniques"
authors:
  - "Miren Hayet-Otero"
  - "Fernando García-García"
  - admin
  - "Joaquín Martínez-Minaya"
  - "Pedro Pablo España Yandiola"
  - "Isabel Urrutia Landa"
  - "Mónica Nieves Ermecheo"
  - "José María Quintana"
  - "Rosario Menéndez"
  - "Antoni Torres"
  - "Rafael Zalacain Jorge"
  - "Inmaculada Arostegui"
date: "2023-04-13"
publishDate: "2023-04-13"
doi: "10.1371/journal.pone.0284150"
publication_types: ["2"]
publication: "*PLOS ONE*"
publication_short: ""
abstract: "With the COVID-19 pandemic having caused unprecedented numbers of infections and deaths, large research efforts have been undertaken to increase our understanding of the disease and the factors which determine diverse clinical evolutions. Here we focused on a fully data-driven exploration regarding which factors (clinical or otherwise) were most informative for SARS-CoV-2 pneumonia severity prediction via machine learning (ML). In particular, feature selection techniques (FS), designed to reduce the dimensionality of data, allowed us to characterize which of our variables were the most useful for ML prognosis. We conducted a multi-centre clinical study, enrolling n = 1548 patients hospitalized due to SARS-CoV-2 pneumonia: where 792, 238, and 598 patients experienced low, medium and high-severity evolutions, respectively. Up to 106 patient-specific clinical variables were collected at admission, although 14 of them had to be discarded for containing ⩾60% missing values. Alongside 7 socioeconomic attributes and 32 exposures to air pollution (chronic and acute), these became d = 148 features after variable encoding. We addressed this ordinal classification problem both as a ML classification and regression task. Two imputation techniques for missing data were explored, along with a total of 166 unique FS algorithm configurations: 46 filters, 100 wrappers and 20 embeddeds. Of these, 21 setups achieved satisfactory bootstrap stability (⩾0.70) with reasonable computation times: 16 filters, 2 wrappers, and 3 embeddeds. The subsets of features selected by each technique showed modest Jaccard similarities across them. However, they consistently pointed out the importance of certain explanatory variables. Namely: patient’s C-reactive protein (CRP), pneumonia severity index (PSI), respiratory rate (RR) and oxygen levels –saturation Sp O2, quotients Sp O2/RR and arterial Sat O2/Fi O2–, the neutrophil-to-lymphocyte ratio (NLR) –to certain extent, also neutrophil and lymphocyte counts separately–, lactate dehydrogenase (LDH), and procalcitonin (PCT) levels in blood. A remarkable agreement has been found a posteriori between our strategy and independent clinical research works investigating risk factors for COVID-19 severity. Hence, these findings stress the suitability of this type of fully data-driven approaches for knowledge extraction, as a complementary to clinical perspectives."
featured: false
url_pdf: "http://arxiv.org/pdf/1202.1056"
url_code: ""
url_dataset: ""
url_poster: ""
url_project: ""
url_slides: ""
url_source: ""
url_video: ""
orcid_put_code: 133058090
---
