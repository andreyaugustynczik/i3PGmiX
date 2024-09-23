# i3PGmiX
i3PGmiX forest growth model is an enhanced version of 3PGmix, built upon the r3PG package (Trotsiuk et al. 2020). i3PGmiX, expands the representation of GPP, NPP, soil processes and natural disturbances in 3PGmix, improving its climate sensitivity and the representation of carbon and water balance in the model.

i3PGmiX integrates the P-model (Stocker et al. 2020) and the P-hydro models (Joshi et al. 2022) of GPP. The latter models are optimality-based Farquhar–von Caemmerer–Berry (FvBC) models of photosynthesis that compute optimal leaf internal to ambient CO2 ratios to balance the costs of maintaining carboxylation and transpiration. Growth respiration is given by fixed fraction of the assimilation and maintenance respiration is a function of fine roots and sapwood biomass, nitrogen content and temperature (Collalti et al. 2016). Furthermore, dark respiration from the GPP models is used to assess the foliage maintenance respiration in the model. 

The management options were also expanded in i3PGmiX, allowing to define thinning interventions not only based on the number of remaining trees, but also based on targets of removal shares of the basal area or the standing volume, as well as absolute target values for the removal of basal area or volume.

Besides the improvements in the representation related to the vegetation and forest management, i3PGmiX includes different soil models, allowing to account for the interactions and feedbacks between soil and vegetation. Specifically, i3PGmiX includes the Yasso20 and ICBM/2N soil carbon models  (Xenakis et al. 2008). 


