######################################################################################
### DATA PREPARATION
# Rename variables: height → Ht, sbp → msbp, dbp → mdbp
# Truncate age and name as age_trun, categorize sex: 1=male, 2=female
# Download 'ht.lms.csv' and 'bp.refcn.csv' to your working directory
######################################################################################

# Load required library
library(data.table)

# 1. Read reference data and convert to data.table
ht.ref <- fread('ht.lms.csv')
bp.ref <- fread('bp.refcn.csv')
setDT(data)  # Convert original data to data.table

# Merge with height reference data
d1 <- merge(data, ht.ref, by = c("age_trun", "sex"), all.x = TRUE)
setDT(d1)  # Ensure d1 is data.table

# 2. Calculate extreme height boundaries (in-place modification)
d1[, `:=`(
  sd2ps = ht_m * ((1 + 2 * ht_l * ht_s) ^ (1 / ht_l)),
  sd3ps = ht_m * ((1 + 3 * ht_l * ht_s) ^ (1 / ht_l)),
  sd3ng = ht_m * ((1 - 3 * ht_l * ht_s) ^ (1 / ht_l)),
  sd2ng = ht_m * ((1 - 2 * ht_l * ht_s) ^ (1 / ht_l))
)]

# 3. Calculate height z-scores
d1[, zht_raw := ((Ht / ht_m) ^ ht_l - 1) / (ht_s * ht_l)]
d1[, zht := fcase(
  zht_raw > 3,  3 + (Ht - sd3ps) / (sd3ps - sd2ps),
  zht_raw < -3, -3 - (sd3ng - Ht) / (sd2ng - sd3ng),
  default = zht_raw
)]

# 4. Convert z-score to height percentiles
d1[, htpct := fcase(
  zht < -1.4632,  5L,
  zht < -0.97802, 10L,
  zht < -0.33724, 25L,
  zht <  0.33724, 50L,
  zht <  0.97802, 75L,
  zht <  1.4632,  90L,
  zht >= 1.4632,  95L,
  default = NA_integer_
)]

# 5. Merge with blood pressure reference data
d2 <- merge(d1, bp.ref, by = c("age_trun", "sex", "htpct"), all.x = TRUE)
setDT(d2)

# Added: cutoff for elevated BP in children/adolescents aged 3–15 y
# According to the guideline: P90–P95 or >=120/80, whichever is lower
d2[, `:=`(
  sbp90_cut = pmin(sbp90, 120),
  dbp90_cut = pmin(dbp90, 80)
)]

# 6. Blood pressure grading
# 6.1 Childhood BP grading (3–15 y)
# Order from high to low because fcase returns the first matched condition
d2[, hbp_ch := fcase(
  msbp >= sbp99 + 5 | mdbp >= dbp99 + 5, 4L,
  (sbp95 <= msbp & msbp < sbp99 + 5) | (dbp95 <= mdbp & mdbp < dbp99 + 5), 3L,
  (sbp90_cut <= msbp & msbp < sbp95) | (dbp90_cut <= mdbp & mdbp < dbp95), 2L,
  msbp < sbp90_cut & mdbp < dbp90_cut, 1L,
  default = NA_integer_
)]

# 6.2 Adolescent BP grading (>=16 y)
# Also use high-to-low ordering to avoid misclassification
d2[, hbp_ad := fcase(
  msbp >= 180 | mdbp >= 110, 5L,
  (160 <= msbp & msbp < 180) | (100 <= mdbp & mdbp < 110), 4L,
  (140 <= msbp & msbp < 160) | (90 <= mdbp & mdbp < 100), 3L,
  (120 <= msbp & msbp < 140) | (80 <= mdbp & mdbp < 90), 2L,
  msbp < 120 & mdbp < 80, 1L,
  default = NA_integer_
)]

# 7. Determine final BP grade (age-specific)
d2[, bp_grade := fifelse(age_trun < 16, hbp_ch, hbp_ad)]

# 8. Define hypertension status based on P95 / 140/90
# <16 y: P95
# >=16 y: 140/90
d2[, `:=`(
  sbp_high95 = fifelse(age_trun < 16, msbp >= sbp95, msbp >= 140),
  dbp_high95 = fifelse(age_trun < 16, mdbp >= dbp95, mdbp >= 90)
)]

# 9. Hypertension phenotype
# 0 = non-HTN
# 1 = ISH
# 2 = IDH
# 3 = SDH
d2[, bp_type := fcase(
  !sbp_high95 & !dbp_high95, 0L,
  sbp_high95 & !dbp_high95, 1L,
  !sbp_high95 &  dbp_high95, 2L,
  sbp_high95 &  dbp_high95, 3L,
  default = NA_integer_
)]

# Optional: labels
d2[, bp_grade_label := fcase(
  bp_grade == 1L, "Normal BP",
  bp_grade == 2L, "Elevated BP",
  bp_grade == 3L, "Grade 1 HTN",
  bp_grade == 4L, "Grade 2 HTN",
  bp_grade == 5L, "Grade 3 HTN",
  default = NA_character_
)]

d2[, bp_type_label := fcase(
  bp_type == 0L, "non-HTN",
  bp_type == 1L, "ISH",
  bp_type == 2L, "IDH",
  bp_type == 3L, "SDH",
  default = NA_character_
)]


# The final dataset is d2