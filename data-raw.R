#
library(meffil)
manifest_EPIC <- meffil.featureset(featureset = "epic")
manifest_450k <- meffil.featureset(featureset = "450k")

use_data(
	manifest_EPIC,
	manifest_450k,
	internal = TRUE,
	overwrite = TRUE
)
