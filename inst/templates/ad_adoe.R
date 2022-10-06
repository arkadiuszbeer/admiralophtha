# Name: ADOE
#
# Label: Ophthalmology Exam Analysis Dataset
#
# Input: adsl, oe
library(admiral)
library(admiral.test) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

# ---- Load source datasets ----

# Use e.g. `haven::read_sas()` to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data


# data("admiral_oe") # when this exists!
# data("admiral_adsl") # can't use this yet as ADSL needs STUDYEYE

load("inst/templates/admiral_adsl.rda")
oe <- load("inst/templates/admiral_oe.rda")

adsl <- adsl %>%
  as.data.frame() %>%
  mutate(STUDYEYE = sample(c("LEFT", "RIGHT"), n(), replace = TRUE)) %>%
  convert_blanks_to_na()

oe <- convert_blanks_to_na(admiral_oe) %>% ungroup()

# ---- Lookup tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~OETESTCD, ~OELAT, ~STUDYEYE, ~PARAMCD, ~PARAM, ~PARAMN,
  "DRSSR", "RIGHT", "RIGHT", "SDRSS", "Diabetic Retinopathy Sev Recode Value-STUDY EYE", 1,
  "DRSSR", "LEFT", "LEFT", "SDRSS", "Diabetic Retinopathy Sev Recode Value-STUDY EYE", 1,
  "DRSSR", "RIGHT", "LEFT", "FDRSS", "Diabetic Retinopathy Sev Recode Value-FELLOW EYE", 2,
  "DRSSR", "LEFT", "RIGHT", "FDRSS", "Diabetic Retinopathy Sev Recode Value-FELLOW EYE", 2,
  "VACSCORE", "RIGHT", "RIGHT", "SBCVA", "Visual Acuity Score-STUDY EYE", 3,
  "VACSCORE", "LEFT", "LEFT", "SBCVA", "Visual Acuity Score-STUDY EYE", 3,
  "VACSCORE", "RIGHT", "LEFT", "FBCVA", "Visual Acuity Score-FELLOW EYE", 4,
  "VACSCORE", "LEFT", "RIGHT", "FBCVA", "Visual Acuity Score-FELLOW EYE", 4
)
attr(param_lookup$OETESTCD, "label") <- "Ophthalmology Test Short Name"

# ---- Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, TRT01A, TRT01P, STUDYEYE)

adoe <- oe %>%
  # Join ADSL with OE (need TRTSDT and STUDYEYE for ADY and PARAMCD derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  # Calculate ADT, ADY
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = OEDTC,
    flag_imputation = "none"
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))

adoe <- adoe %>%
  # Add PARAM, PARAMCD
  derive_vars_merged(
    dataset_add = param_lookup,
    new_vars = vars(PARAM, PARAMCD),
    by_vars = vars(OETESTCD, OELAT, STUDYEYE)
  ) %>%
  # Calculate AVAL and AVALC
  mutate(
    AVAL = OESTRESN,
    AVALC = OESTRESC
  )

# Derive visit info - no ATPT and ATPTN yet as SDTM variables not in test data
adoe <- adoe %>%
  mutate(
    # ATPTN = OETPTNUM,
    # ATPT = OETPT,
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN") ~ "Screening",
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),

    AVISITN = case_when(
      AVISIT == "Baseline" ~ 0,
      !is.na(VISITNUM) ~ round(VISITNUM,0)
    )

  )

# Derive DTYPE and BASETYPE
adoe <- adoe %>%
  mutate(
    DTYPE = NA_character_,
    BASETYPE = "LAST"
  )

# Derive Treatment flags
adoe <- adoe %>%
  # Calculate ONTRTFL
  derive_var_ontrtfl(
    start_date = ADT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    filter_pre_timepoint = AVISIT == "Baseline"
  ) %>%
  # Calculate ABLFL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
      order = vars(ADT, VISITNUM, OESEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(BASETYPE))
  )


# Derive baseline information
adoe <- adoe %>%
  # Calculate BASE
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE
  ) %>%
  # Calculate BASEC
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVALC,
    new_var = BASEC
  ) %>%
  # Calculate CHG
  derive_var_chg() %>%
  # Calculate PCHG
  derive_var_pchg()


# ANL01FL: Flag last result within an AVISIT and ATPT for post-baseline records
adoe <- adoe %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      new_var = ANL01FL,
      by_vars = vars(USUBJID, PARAMCD, AVISIT, DTYPE),
      order = vars(ADT, AVAL),
      mode = "last"
    ),
    filter = !is.na(AVISITN) & ONTRTFL == "Y"
  )

# Get treatment information
advs <- advs %>%
  # Assign TRTA, TRTP
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  ) %>%
  # Create End of Treatment Record
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(STUDYID, USUBJID, PARAMCD, ATPTN),
      order = vars(ADT),
      new_var = EOTFL,
      mode = "last"
    ),
    filter = (4 < VISITNUM &
                VISITNUM <= 13 & ANL01FL == "Y" & is.na(DTYPE))
  ) %>%
  filter(EOTFL == "Y") %>%
  mutate(
    AVISIT = "End of Treatment",
    AVISITN = 99
  ) %>%
  union_all(advs) %>%
  select(-EOTFL)

# Get ASEQ and AVALCATx and add PARAM/PARAMN
advs <- advs %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = vars(STUDYID, USUBJID),
    order = vars(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN, DTYPE),
    check_type = "error"
  ) %>%
  # Derive AVALCA1N and AVALCAT1
  mutate(AVALCA1N = format_avalcat1n(param = PARAMCD, aval = AVAL)) %>%
  derive_vars_merged(dataset_add = avalcat_lookup, by_vars = vars(PARAMCD, AVALCA1N)) %>%
  # Derive PARAM and PARAMN
  derive_vars_merged(dataset_add = select(param_lookup, -VSTESTCD), by_vars = vars(PARAMCD))


# Add all ADSL variables
advs <- advs %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = vars(STUDYID, USUBJID)
  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

# ---- Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
save(advs, file = file.path(dir, "advs.rda"), compress = "bzip2")
