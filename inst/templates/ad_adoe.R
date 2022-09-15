# Name: ADOE
#
# Label: Ophthalmology Exam Analysis Dataset Dataset
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

adsl <- adsl_sas
oe <- oe

oe <- convert_blanks_to_na(oe)

# ---- Lookup tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~OETESTCD, ~OELAT, ~STUDYEYE, ~PARAMCD, ~PARAM, ~PARAMN,
  "CSUBTH", "RIGHT", "RIGHT", "SCSUBTH", "Center Subfield Thickness (mmHg)-STUDY EYE", 1,
  "CSUBTH", "LEFT", "LEFT", "SCSUBTH", "Center Subfield Thickness (mmHg)-STUDY EYE", 1,
  "CSUBTH", "RIGHT", "LEFT", "FCSUBTH", "Center Subfield Thickness (mmHg)-FELLOW EYE", 2,
  "CSUBTH", "LEFT", "RIGHT", "FCSUBTH", "Center Subfield Thickness (mmHg)-FELLOW EYE", 2,
  "ASCAN", "RIGHT", "RIGHT", "SASCAN", "A-Scan-STUDY EYE", 3,
  "ASCAN", "LEFT", "LEFT", "SASCAN", "A-Scan-STUDY EYE", 3,
  "ASCAN", "RIGHT", "LEFT", "FASCAN", "A-Scan-FELLOW EYE", 4,
  "ASCAN", "LEFT", "RIGHT", "FASCAN", "A-Scan-FELLOW EYE", 4,
  "BSCAN", "RIGHT", "RIGHT", "SBSCAN", "B-Scan-STUDY EYE", 5,
  "BSCAN", "LEFT", "LEFT", "SBSCAN", "B-Scan-STUDY EYE", 5,
  "BSCAN", "RIGHT", "LEFT", "FBSCAN", "B-Scan-FELLOW EYE", 6,
  "BSCAN", "LEFT", "RIGHT", "FBSCAN", "B-Scan-FELLOW EYE", 6
)
attr(param_lookup$OETESTCD, "label") <- "Ophthalmology Test Short Name"

# # Assign ANRLO/HI, A1LO/HI
# range_lookup <- tibble::tribble(
#   ~PARAMCD, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
#   "SYSBP", 90, 130, 70, 140,
#   "DIABP", 60, 80, 40, 90,
#   "PULSE", 60, 100, 40, 110,
#   "TEMP", 36.5, 37.5, 35, 38
# )
# # ASSIGN AVALCAT1
# avalcat_lookup <- tibble::tribble(
#   ~PARAMCD, ~AVALCA1N, ~AVALCAT1,
#   "HEIGHT", 1, ">100 cm",
#   "HEIGHT", 2, "<= 100 cm"
# )
#
# # ---- User defined functions ----
#
# # Here are some examples of how you can create your own functions that
# #  operates on vectors, which can be used in `mutate()`.
# format_avalcat1n <- function(param, aval) {
#   case_when(
#     param == "HEIGHT" & aval > 140 ~ 1,
#     param == "HEIGHT" & aval <= 140 ~ 2
#   )
# }

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


# Derive visit info - requires updating once we get oe test data
adoe <- adoe %>%
  # Derive Timing
  mutate(
    ATPTN = OETPTNUM,
    ATPT = OETPT,
    AVISIT = case_when(
      str_detect(str_to_upper(VISIT), "SCREEN|RETRIEVAL|AMBUL|TIMEPOINT|DELAYED") ~ NA_character_,
      VISIT == "Day 1" ~ "Baseline",
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = as.numeric(case_when(
      AVISIT == "Baseline" ~ "0",
      AVISIT == "Early Termination" ~ "299",
      AVISIT == "Unscheduled" ~ "999",
      str_detect(AVISIT, "Day") ~ str_trim(str_replace(AVISIT, "Day", "")),
      TRUE ~ NA_character_
    ))
  )

adoe <- adoe %>%
  # Calculate ONTRTFL
  derive_var_ontrtfl(
    start_date = ADT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    filter_pre_timepoint = AVISIT == "Baseline"
  )

############### END OF WORK SO FAR ####################
# need to derive: ANL01FL, ANL02FL, LAST01FL, WORS01FL, and some of the stuff below


# Derive baseline flags
advs <- advs %>%
  # Calculate BASETYPE
  derive_var_basetype(
    basetypes = rlang::exprs(
      "LAST: AFTER LYING DOWN FOR 5 MINUTES" = ATPTN == 815,
      "LAST: AFTER STANDING FOR 1 MINUTE" = ATPTN == 816,
      "LAST: AFTER STANDING FOR 3 MINUTES" = ATPTN == 817,
      "LAST" = is.na(ATPTN)
    )
  ) %>%
  # Calculate ABLFL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
      order = vars(ADT, VISITNUM, VSSEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = (!is.na(AVAL) &
                ADT <= TRTSDT & !is.na(BASETYPE) & is.na(DTYPE))
  )

# Derive baseline information
advs <- advs %>%
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
  # Calculate BNRIND
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ANRIND,
    new_var = BNRIND
  ) %>%
  # Calculate CHG
  derive_var_chg() %>%
  # Calculate PCHG
  derive_var_pchg()


# ANL01FL: Flag last result within an AVISIT and ATPT for post-baseline records
advs <- advs %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      new_var = ANL01FL,
      by_vars = vars(USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
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
