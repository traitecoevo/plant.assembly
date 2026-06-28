# Shared model_support for integration tests that run the SCM. Keep
# max_patch_lifetime small (30) so each solve takes a few seconds, matching the
# convention used in test-community.R / test-support-fitness.R.
assembly_model_support <- function(max_patch_lifetime = 30) {
  list(p = plant_default_assembly_pars(max_patch_lifetime = max_patch_lifetime),
       plant_control = plant_default_assembly_control())
}
