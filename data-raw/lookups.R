
# data-raw/lookups.R

## code to prepare `ESCDAGs_causal_criteria_df`
name <- c("Temporality",
          "Face-validity",
          "Recourse to theory"
          )
question <- c("Does the variable to the left of the arrow precede the variable on the right?",
              "Is the posited relationship plausible?",
              "Is the posited relationship supported by theory?"
              )
description <- c("\'Of the Bradford Hill criteria, temporality is the only one not requiring extensive qualification or not yet disproven. (Glass et al., 2013; DOI: https://doi.org/10.1146/annurev-publhealth-031811-124606) It states that effect cannot precede cause. For example, in Figure 1(A) (Ferguson et al., 2020; DOI: https://doi.org/10.1093/ije/dyz150), adolescent substance use cannot precede historical parental alcohol use, so the relationship would not be temporal. Unless the directed edge is not temporal, we proceed to causal criterion 2.\' (Ferguson et al., 2020, p. 326)",
                 "\'Face-validity is related to the Bradford Hill criterion of (biologic) plausibility. Nested within the wider causal criteria scheme, the face-validity criterion is a rapid means of using reviewer background knowledge to identify implausible relationships, given the temporality established in criterion 1. For example, in Figure 1(A) it is plausible that directed edges originate from sex, but implausible that historical parental alcohol use could influence adolescent sex assignment despite temporal ordering.\' (Ferguson et al., 2020, p. 326)",
                 "\'The recourse to theory criterion considers background and expert knowledge more overtly. It subsumes the temporality and face-validity criteria and continues to cement a platform for the counterfactual thought experiment. Where the face-validity criterion is concerned with the researcher’s own knowledge, the step assesses whether there is formal theoretical support for the relationship. The decision log for this criterion requires the reviewer to state briefly what theory applies (if any) with space for a reference. As noted above, lack of theory does not equate to lack of effect. As such the purpose of this criterion is not so much falsification as preparation for the next step.\' (Ferguson et al., 2020, p. 326)"
                 )
source <- c("Ferguson et al., 2020, 'Evidence synthesis for constructing directed acyclic graphs (ESC-DAGs): a novel and systematic method for building directed acyclic graphs', pp. 324-326, DOI: https://doi.org/10.1093/ije/dyz150)",
            "Ferguson et al., 2020, 'Evidence synthesis for constructing directed acyclic graphs (ESC-DAGs): a novel and systematic method for building directed acyclic graphs', pp. 324-326, DOI: https://doi.org/10.1093/ije/dyz150)",
            "Ferguson et al., 2020, 'Evidence synthesis for constructing directed acyclic graphs (ESC-DAGs): a novel and systematic method for building directed acyclic graphs', pp. 324-326, DOI: https://doi.org/10.1093/ije/dyz150)"
            )
required <- c("yes", "yes", "no")

ESCDAGs <- data.frame(name,
                      question,
                      description,
                      source,
                      required,
                      stringsAsFactors = FALSE)

#usethis::use_data(ESCDAGs, overwrite = TRUE)
