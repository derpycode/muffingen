#
#
# --- CLIMATE ---------------------------------------------------------
#
# scaling for atmospheric CO2 radiative forcing, calculated as the multiplier of 278 ppm CO2
ea_radfor_scl_co2=1.0
#
# --- BIOLOGICAL NEW PRODUCTION ---------------------------------------
#
# Set no biological production
bg_par_bio_prodopt='NONE'
#
# --- DATA SAVING -----------------------------------------------------
#
# Set BASIC + tracer + full physics save level (12)
# (assume default time-series and time-slice save intervals)
bg_par_data_save_level=12
#
# --- MISC ------------------------------------------------------------
#
# activate single-tracer ventillation age tracer
# NOTE: you *must* have the 'red' color tracer incuded in the base-config and hence gm_ocn_select_4=.true.
bg_ctrl_force_ocn_age1=.true.
#
# --- END -------------------------------------------------------------
#
#