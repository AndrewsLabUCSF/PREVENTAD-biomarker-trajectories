
# BIOMARKER OUTCOME LABELS ------------------------------------------------

# Short labels — used in trajectory/facet plots
OUTCOME_LABELS <- c(
  ab_ratio               = "Aβ42/40",
  gfap_log               = "GFAP",
  nfl_log                = "NfL",
  ptau181_log            = "ptau181",
  ptau217_log            = "ptau217",
  ptau217_ab42_ratio_log = "ptau217/Aβ42 ratio"
)

# Long labels with (log) annotation — used in forest/coefficient plots
OUTCOME_LABELS_LOG <- c(
  ab_ratio               = "Aβ42/40",
  gfap_log               = "GFAP (log)",
  nfl_log                = "NfL (log)",
  ptau181_log            = "ptau181 (log)",
  ptau217_log            = "ptau217 (log)",
  ptau217_ab42_ratio_log = "ptau217/Aβ42 (log)"
)


# BIOMARKER ORDERING (for gt tables) --------------------------------------

BIOMARKER_ORDER_LME <- c(
  "Aβ42/40", "GFAP", "NfL", "ptau181", "ptau217", "ptau217/Aβ42 ratio"
)

BIOMARKER_ORDER_LME_LOG <- c(
  "Aβ42/40", "GFAP (log)", "NfL (log)",
  "ptau181 (log)", "ptau217 (log)", "ptau217/Aβ42 (log)"
)


# SPARSE OUTCOMES ---------------------------------------------------------

SPARSE_OUTCOMES <- c("ptau217_ab42_ratio_log")


# MODEL COLOURS -----------------------------------------------------------

COLOR_COGDRISK <- "#374E55FF"
COLOR_LIBRA2   <- "#DF8F44FF"

MODEL_COLORS <- c(CogDrisk = COLOR_COGDRISK, LIBRA2 = COLOR_LIBRA2)
