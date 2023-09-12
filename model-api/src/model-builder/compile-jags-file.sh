#!/usr/bin/env bash

# User- or pipeline-supplied args.
LIKELIHOOD="${1}"
DETERMINISTIC_MODEL="${2}"
GROUP_LEVEL_EFFECTS="${3}"
DESIGN_MATRIX="${4}"
MOD_DIR="${5}"
KEEP_PRED_SWITCHES="${6:-drop_pred_switches}"
GLE_ZEROS="${7:-none}"
VAR_LEVEL="${8:-stratum}"
VAR_TYPE="${9:-fixed}"
DERIVED_QUANTITIES="${10:-none}"
CENSORING="${11:-no}"
TRUNCATION="${12:-no}"
PARAMETERIZATION="${13:-centered}"
DEFLECTIONS="${14:-no}"
B_DROP="${15:-keepit}"
G_DROP="${16:-keepit}"
FINITE_POP="${17:-no}"
OMIT_PRED_STRAT_BLOCK="${18:-no}"
GET_MARG_EFFECTS="${19:-no}"
RANDOM_BETA="${20:-no}" # new!
APPLY_OFFSET="${21:-no}"

# Variables for oft-recycled paths.
CODE_BLOCKS="model-api/src/model-builder/code-blocks"
PRIORS_BLOCK="${CODE_BLOCKS}/prior"
GEN_QUANT_BLOCK="${CODE_BLOCKS}/derived-quantities"
MOD_CHECK_BLOCK="${CODE_BLOCKS}/model-checking"
SLI="${GEN_QUANT_BLOCK}/stratum-level-inference"

if [ "${APPLY_OFFSET}" == "yes" ]; then
	EXPOSURE="* exposure[n]"
else
	EXPOSURE=""
fi

# Dropouts suffix switch.
DROP_SUFFIX=""
if [[ "${B_DROP}" = *"-int-only"* && "${G_DROP}" == "keepit" ]]; then
	DROP_SUFFIX="-b-drop"
elif [[ "${B_DROP}" = *"-int-only"* && "${G_DROP}" = *"-int-only"* ]]; then
	DROP_SUFFIX="-bg-drop"
elif [[ "${B_DROP}" == "keepit" && "${G_DROP}" = *"-int-only"* ]]; then
	DROP_SUFFIX="-g-drop"
fi
B_DROP_ADDENDUM=""
if [[ "${B_DROP}" == "-int-only-fixed" ]]; then
	B_DROP_ADDENDUM="-fixed-b0"
fi
G_DROP_ADDENDUM=""
if [[ "${G_DROP}" == "-int-only-fixed" ]]; then
	G_DROP_ADDENDUM="-fixed-g0"
fi

# Other likelihood-related derived quantities.
L_RELATED_DQS_BLOCK="${CODE_BLOCKS}/likelihood/related-dqs"
LIKELIHOOD_RELATED_DQS_SWITCH_FILE="${L_RELATED_DQS_BLOCK}/.empty"
GAP_SIZE_THRESH_SWITCH_FILE="${L_RELATED_DQS_BLOCK}/.empty"
if [ "${CENSORING}" == "yes" ]; then  # hacky, do fix this
	LIKELIHOOD_RELATED_DQS_SWITCH_FILE="${L_RELATED_DQS_BLOCK}/prop-gt-thresh"
	GAP_SIZE_THRESH_SWITCH_FILE="${L_RELATED_DQS_BLOCK}/prop-landscape"
fi

if [ "${RANDOM_BETA}" == "yes" ]; then
	GROUP_LEVEL_EFFECTS="${GROUP_LEVEL_EFFECTS}-Beta"
fi

HAT_SITE_VARIANCE="${GEN_QUANT_BLOCK}/hat-site-means/.empty"
HAT_STRATUM_VARIANCE="${GEN_QUANT_BLOCK}/stratum-level-inference/.empty"
if [ "${LIKELIHOOD}" == "negative-binomial" ]; then
	HAT_SITE_VARIANCE="hat.site.sigma.y[j, t, k] <- sigma.ss[j, k]"
	HAT_STRATUM_VARIANCE="hat.strat.sigma.y[t, k] <- sigma.ss[j.draw[k], k]"
fi

# Covariate-related conditionals.
PRED_SITE_IDX="1:length(x.pred)"
STRATUM_SWITCH_FILE_SUFFIX=""
VAR_MOD_KIND=""
if [ "${LIKELIHOOD}" == "hurdle-ordinal-latent-beta" ]; then
	STRATUM_SWITCH_FILE_SUFFIX="-hurdle"
	VAR_MOD_KIND="-01"
fi
if [ "${LIKELIHOOD}" == "lognormal" ]; then
	STRATUM_SWITCH_FILE_SUFFIX="-${LIKELIHOOD}"
fi
FINITE_POP_SUFFIX=""
BONUS_SUFFIX=""  # TODO: deal with this later -- a tough problem for sure!
if [ "${FINITE_POP}" == "yes" ]; then
	PRED_SITE_IDX="in.sample.idx"
	FINITE_POP_SUFFIX="-finite-pop"
	if [ "${OMIT_PRED_STRAT_BLOCK}" == "yes" ]; then
			BONUS_SUFFIX="-tbd"
	fi
	HAT_STRATUM_OOS_SWITCH_FILE="${SLI}/hat-gles/.empty"
else
	HAT_STRATUM_OOS_SWITCH_FILE="${SLI}/hat-gles/${GROUP_LEVEL_EFFECTS}"
fi

if [ "${GROUP_LEVEL_EFFECTS}" == "b0-b1" ]; then
	HAT_STRAT_LIN_PRED="B.tilde[1, k] + B.tilde[2, k] * x.hat[t]"
	PRED_STRAT_LIN_PRED="B.tilde[1, k] + B.tilde[2, k] * x.pred[i.pred[t, k]]"
	DRIVER_STRAT_LIN_PRED="B.tilde[1, k] + X.driver[s, ] %*% Beta[which.drivers]"
	_B0_4_SITE_FROM_ZONE_='B.tilde[1,'
	_B1_4_SITE_FROM_ZONE_='B.tilde[2,'
elif [ "${GROUP_LEVEL_EFFECTS}" == "b0" ]; then
	HAT_STRAT_LIN_PRED="B0.tilde[k] + B1[k] * x.hat[t]"
	PRED_STRAT_LIN_PRED="B0.tilde[k] + B1[k] * x.pred[i.pred[t, k]]"
	DRIVER_STRAT_LIN_PRED="B0.tilde[k] + X.driver[s, ] %*% Beta[which.drivers]"
	_B0_4_SITE_FROM_ZONE_='B0.tilde['
	_B1_4_SITE_FROM_ZONE_='B1['
elif [ "${GROUP_LEVEL_EFFECTS}" == "b0-Beta" ]; then
	PLACEHOLDER="todo"
	# PLACEHOLDER for work to come. Will need a block for b0-b1-Beta as well.
else
	GROUP_LEVEL_EFFECTS='no-gle'
	HAT_STRAT_LIN_PRED="B0[k] + B1[k] * x.hat[t]"
	PRED_STRAT_LIN_PRED="B0[k] + B1[k] * x.pred[i.pred[t, k]]"
	DRIVER_STRAT_LIN_PRED="B0[k] + X.driver[s, ] %*% Beta[which.drivers]"
	_B0_4_SITE_FROM_ZONE_='B0['
	_B1_4_SITE_FROM_ZONE_='B1['
fi
DRIVER_MARGINAL_EFFECTS=".empty"
THIS_BETA="Beta"
RANDOM_BETA_DEF=".empty"
if [ "${DESIGN_MATRIX}" == "null" ]; then
	X='without-covariates'
	X_DESC="not present"
	ACP_TEMPLATE_FILE="${PRIORS_BLOCK}/.empty"
else
	X='with-covariates'
	X_DESC="present"
	ACP_TEMPLATE_FILE="${PRIORS_BLOCK}/_with-covariates-template_"
	if [ "${LIKELIHOOD}" == "hurdle-ordinal-latent-beta" ]; then
		ACP_TEMPLATE_FILE="${PRIORS_BLOCK}/with-covariates-hurdle-ordinal-latent-beta${DROP_SUFFIX}"
	fi
	if [ "${GET_MARG_EFFECTS}" == "yes" ]; then
		DRIVER_MARGINAL_EFFECTS="driver-stratum${STRATUM_SWITCH_FILE_SUFFIX}${DROP_SUFFIX}${FINITE_POP_SUFFIX}"
		if [ "${LIKELIHOOD}" == "ordinal-latent-normal" ]; then
			DRIVER_MARGINAL_EFFECTS="driver-stratum-oln${FINITE_POP_SUFFIX}"
		fi
		if [ "${LIKELIHOOD}" == "gen-pois" ]; then
			DRIVER_MARGINAL_EFFECTS="driver-stratum-gen-pois${FINITE_POP_SUFFIX}"
		fi
		if [ "${LIKELIHOOD}" == "zero-inflated-binomial" ] || \
			 [ "${LIKELIHOOD}" == "zero-inflated-beta-binomial" ]; then
			 DRIVER_MARGINAL_EFFECTS="driver-stratum-zi-binom${FINITE_POP_SUFFIX}"
  	fi
		if [ "${LIKELIHOOD}" == "zero-inflated-poisson" ] || \
			 [ "${LIKELIHOOD}" == "zero-inflated-negative-binomial" ]; then
			 DRIVER_MARGINAL_EFFECTS="driver-stratum-zi-counts${FINITE_POP_SUFFIX}"
		fi
	fi

	if [ "${RANDOM_BETA}" == "yes" ]; then
		THIS_BETA="Beta.n[n, ]"
		RANDOM_BETA_DEF="Beta-n"
	fi

fi

PSM_TEMPLATE_FIILENAME=".empty"
PRED_PARK_SWITCH_FILE=".empty"
HAT_STRATUM_SWITCH_FILE="hat-stratum${STRATUM_SWITCH_FILE_SUFFIX}${DROP_SUFFIX}${FINITE_POP_SUFFIX}"
PRED_STRATUM_SWITCH_FILE=".empty"
if [ "${KEEP_PRED_SWITCHES}" != 'drop_pred_switches' ]; then
	PSM_TEMPLATE_FIILENAME="TBD"
	PRED_PARK_SWITCH_FILE="pred-park"
	PRED_STRATUM_SWITCH_FILE="pred-stratum${STRATUM_SWITCH_FILE_SUFFIX}${DROP_SUFFIX}${FINITE_POP_SUFFIX}${BONUS_SUFFIX}"
fi

HAT_SITE_MEANS_TEMPFILE="${MOD_DIR}/hat-site-means-block"
PRED_SITE_MEANS_TEMPFILE="${MOD_DIR}/pred-site-means-block"
GLE_PRIORS_TEMPFILE="${MOD_DIR}/gle-priors-block"
ADDED_COV_PRIORS_TEMPFILE="${MOD_DIR}/added-cov-priors-block"

HAT_STRATUM_ADDITIONAL_QUANTS_FILE="${SLI}/.empty"
PRED_STRATUM_ADDITIONAL_QUANTS_FILE="${SLI}/.empty"

# Create the directory to which the model definition will be written.
mkdir -p "${MOD_DIR}"
# echo ${MOD_DIR}

# Likelihood-related conditionals.
LINK="("
Y_VEC="y"
CENSOR_SWITCH_LIKELIHOOD_FILE="${CODE_BLOCKS}/likelihood/.empty"
CENSOR_SWITCH_CHECKING_FILE="${CODE_BLOCKS}/model-checking/.empty"
TRUNCATION_STRING=""
HAT_QUANT_SUFFIX="mean"
DISTRIBUTION_PARAMS_SWITCH_FILE="${LIKELIHOOD}"
MISC_INFERENCE_FILE=".empty"
CHANGE_IN_ODDS=".empty"
Y_SWITCH='y'
HIL_TEMPLATE_FILE="${CODE_BLOCKS}/likelihood/.empty"
HAT_PARK_FILE="hat-park"

MAX_R_SWITCH_FILE=".empty"
COEF_PRIOR_INTERCEPT="4E-06"  #"1/500^2"  # was ".00001"
COEF_PRIOR_SLOPE="1E-04"  # "1/100^2"
VAR_PRIOR_UL="100"
TAU_SITE_FILE=".empty"

# Hierarchical site-level variances priors (defaults for inv logit models).
MU_G_SHAPE="0.01" #"49"
MU_G_RATE="0.01" #"28"
SIGMA_G_SHAPE="0.01" #"1"
SIGMA_G_RATE="0.01" #"4"
MU_SIGMA_PRIOR_DIST="dgamma"
SIGMA_SIGMA_PRIOR_DIST="dgamma"

# Fixed stratum-level variances constants.
SIGMA_K_UPPER="10"
SIGMA_B0_PRIOR="dgamma(0.01, 0.01)"  # was: "dunif(1E-6, 1000)"
SIGMA_B1_PRIOR="dgamma(0.01, 0.01)"  # was: "dunif(1E-6, 100)"
if [ "${DETERMINISTIC_MODEL}" == "exponential" ]; then
	 	SIGMA_B0_PRIOR="dunif(1E-6, 100)"
		SIGMA_B1_PRIOR="dunif(1E-6, 100)"
		COEF_PRIOR_INTERCEPT="0.0001"
		COEF_PRIOR_SLOPE="0.001"
fi

B_SUFFIX=""
if [[ "${B_DROP}" = *"-int-only"* ]]; then
	B_SUFFIX="-int-only"  # was: "${B_DROP}"
fi
G_SUFFIX=""
HSM_TEMPLATE_FIILENAME="_template_"
HAT_PNORM="${GEN_QUANT_BLOCK}/zone-means/.empty"
PRED_PNORM="${GEN_QUANT_BLOCK}/zone-means/.empty"
PHI_DEF="${CODE_BLOCKS}/deterministic-model/.empty"
if [ "${LIKELIHOOD}" == "poisson" ]; then
	DISTRIBUTION_PARAMS_SWITCH_FILE=".empty"
elif [[ "${LIKELIHOOD}" = *"hurdle"* ]]; then
	Y_SWITCH="y.beta"  # may be extended to, e.g., count data models
	SIGMA_K_UPPER="1"
	HIL_TEMPLATE_FILE="${CODE_BLOCKS}/likelihood/hurdle/_bernoulli-template_"
	if [[ "${G_DROP}" = *"-int-only"* ]]; then
		G_SUFFIX="-int-only" # was: G_SUFFIX="${G_DROP}"
	fi
	PHI_DEF="${CODE_BLOCKS}/deterministic-model/${X}/hurdle-inverse-logit${G_SUFFIX}"
#elif [ "${LIKELIHOOD}" == "zero-inflated-poisson" ] || \
#	 [ "${LIKELIHOOD}" == "zero-inflated-negative-binomial" ]; then
#	 	DETERMINISTIC_MODEL="${DETERMINISTIC_MODEL}-mbh"
elif [ "${LIKELIHOOD}" == "zero-inflated-beta-binomial" ]; then
	 	SIGMA_B0_PRIOR="dunif(1E-6, 100)"
		SIGMA_B1_PRIOR="dunif(1E-6, 50)"
elif [ "${LIKELIHOOD}" == "beta" ]; then
	DETERMINISTIC_MODEL="restricted-inverse-logit"
elif [ "${LIKELIHOOD}" == "binomial" ] || \
	 [ "${LIKELIHOOD}" == "beta-binomial" ] || \
	 [ "${LIKELIHOOD}" == "zero-inflated-binomial" ] || \
	 [ "${LIKELIHOOD}" == "zero-inflated-beta-binomial" ]; then
		MISC_INFERENCE_FILE="binomial"
		CHANGE_IN_ODDS="${GROUP_LEVEL_EFFECTS}"
elif [ "${LIKELIHOOD}" == "ordinal-latent-normal" ]; then
	HAT_PNORM="${GEN_QUANT_BLOCK}/zone-means/hat-pnorm"
	PRED_PNORM="${GEN_QUANT_BLOCK}/zone-means/pred-pnorm"
	HAT_QUANT_SUFFIX="mu"
	DETERMINISTIC_MODEL="linear" # was: linear-oln
	MU_G_SHAPE="0" #".01"
	MU_G_RATE="3" #".01"
	SIGMA_G_SHAPE="0" #".01"
	SIGMA_G_RATE="3" #".01"
	MU_SIGMA_PRIOR_DIST="dunif" #"dgamma"
	SIGMA_SIGMA_PRIOR_DIST="dunif" #"dgamma"
	VAR_PRIOR_UL="10"
	TAU_SITE_FILE="tau-site"
fi

if [ "${DETERMINISTIC_MODEL}" == "restricted-linear" ]; then
		LINK="max(0.00001,"
	  MU_G_SHAPE=".01"
	 	MU_G_RATE=".01"
	 	SIGMA_G_SHAPE=".01"  #"0"
	 	SIGMA_G_RATE=".01"  #"100"
		SIGMA_K_UPPER="500"
elif [ "${DETERMINISTIC_MODEL}" == "exponential" ]; then
	MU_G_SHAPE=".01"
	MU_G_RATE=".01"
	SIGMA_G_SHAPE=".01"  #"0"
	SIGMA_G_RATE=".01"  #"100"
	SIGMA_K_UPPER="500"
	MU_SIGMA_PRIOR_DIST="dunif"  #"dgamma"
	SIGMA_SIGMA_PRIOR_DIST="dunif"  #"dgamma"
fi
if [ "${LIKELIHOOD}" == "gamma" ]; then
	SIGMA_K_UPPER="500"
fi

# If a uniform prior is used for site level variances, redefine the numerical
# arguments to the prior distribution.
if [ "${MU_SIGMA_PRIOR_DIST}" == "dunif" ]; then
	MU_G_SHAPE="0"
	MU_G_RATE="500" # was 100
	if [ "${LIKELIHOOD}" == "ordinal-latent-normal" ]; then
		MU_G_RATE="5"
	fi
fi
if [ "${SIGMA_SIGMA_PRIOR_DIST}" == "dunif" ]; then
	SIGMA_G_SHAPE="0"
	SIGMA_G_RATE="500" # was 100
	if [ "${LIKELIHOOD}" == "ordinal-latent-normal" ]; then
		SIGMA_G_RATE="5"
	fi
fi

DET_FUN_PREFIX=""
if [ "${DETERMINISTIC_MODEL}" == "linear" ]; then
	if [[ "${LIKELIHOOD}" = *"poisson"* ]] || \
		 [[ "${LIKELIHOOD}" = *"negative-binomial"* ]]; then
			 DET_FUN_PREFIX="restricted-"
			 LINK="max(0.00001,"
	fi
fi

if [ "${CENSORING}" == "yes" ]; then
	Y_VEC="y"  #was "censor.limit.vec"
	CENSOR_SWITCH_LIKELIHOOD_FILE="${CODE_BLOCKS}/likelihood/censored-data"
	CENSOR_SWITCH_CHECKING_FILE="${CODE_BLOCKS}/model-checking/censored-data"
fi
if [ "${TRUNCATION}" != "no" ]; then
	TRUNCATION_STRING="T(${TRUNCATION})"
fi

# Conditionals related to the deterministic model.
if [[ "${DETERMINISTIC_MODEL}" = *"exponential"* ]]; then
	LINK="exp("
elif [[ "${DETERMINISTIC_MODEL}" = *"inverse-logit"* ]]; then
	LINK="ilogit("
	COEF_PRIOR_INTERCEPT="pow(1.5, -2)"
	COEF_PRIOR_SLOPE="16"  #"1/0.25^2"
	VAR_PRIOR_UL="1"  # was "10"
	if [ "${PARAMETERIZATION}" == "non-centered" ]; then
		VAR_PRIOR_UL="1.5"
	fi
elif [ "${DETERMINISTIC_MODEL}" == "monomolecular" ]; then
	MAX_R_SWITCH_FILE="max-r"
	HSM_TEMPLATE_FIILENAME="monomolecular"
	PSM_TEMPLATE_FIILENAME="monomolecular"
	HAT_QUANT_SUFFIX="lp"
	HAT_STRATUM_ADDITIONAL_QUANTS_FILE="${SLI}/additional-quantities/hat-stratum/monomolecular"
	PRED_STRATUM_ADDITIONAL_QUANTS_FILE="${SLI}/additional-quantities/pred-stratum/monomolecular"
fi

DEFLECTIONS_FILE="${GEN_QUANT_BLOCK}/deflections/.empty"
MU_D_BLOCK="${GEN_QUANT_BLOCK}/deflections/mu-D-Beta"
HAT_STRAT_DEFL_MEAN="${SLI}/.empty"
PRED_STRAT_DEFL_MEAN="${SLI}/.empty"
TREND_ZONE_FILE=".empty"
if [ "${DEFLECTIONS}" == "yes" ]; then
	DEFLECTIONS_FILENAME="deflections-Beta"
	if [ "${RANDOM_BETA}" == "yes" ]; then
		DEFLECTIONS_FILENAME="deflections-Beta-tilde"
		MU_D_BLOCK="${GEN_QUANT_BLOCK}/deflections/mu-D-Beta-tilde"
	fi
	DEFLECTIONS_FILE="${GEN_QUANT_BLOCK}/deflections/${DEFLECTIONS_FILENAME}"
	if [ "${FINITE_POP}" == "yes" ]; then
		MU_D_BLOCK="${GEN_QUANT_BLOCK}/deflections/.empty"
		HAT_STRAT_DEFL_MEAN="${SLI}/hat-stratum-defl-mean"
		PRED_STRAT_DEFL_MEAN="${SLI}/pred-stratum-defl-mean"
		TREND_ZONE_FILE="zone"
	fi
fi

# Misc. control flow for hurdle models (defaults followed by replacements).
GLE_ZEROS_SUFFIX=""
GLE_ZEROS_JAGSFILE_SUFFIX=""
EPSILON_FILE="${DETERMINISTIC_MODEL}"
HAT_ZI_ZONE_QUANTS="${SLI}/.empty"
PRED_ZI_ZONE_QUANTS="${SLI}/.empty"
if [ "${LIKELIHOOD}" == "zero-inflated-binomial" ] || \
   [ "${LIKELIHOOD}" == "zero-inflated-beta-binomial" ]; then
		 HAT_QUANT_SUFFIX="mu.p"
		 HAT_PARK_FILE="zi-hat-park"
		 HSM_TEMPLATE_FIILENAME="_zi-template_"
		 EPSILON_FILE="zi-inverse-logit"
		 HAT_ZI_ZONE_QUANTS="${SLI}/hat-zone-z"
		 PRED_ZI_ZONE_QUANTS="${SLI}/pred-zone-z"
fi
BAYESIAN_P_VALUES_SWITCH_FILE="bayesian-p-values"
MODEL_CHECKING_COMMENT_FILE="_comment_"
if [[ "${LIKELIHOOD}" = *"hurdle"* ]]; then
	GLE_ZEROS_SUFFIX="-${GLE_ZEROS}"
	GLE_ZEROS_JAGSFILE_SUFFIX="_${GLE_ZEROS}"
	EPSILON_FILE="over-hurdle-latent-beta"  # may be extended to, e.g., count data models
	BAYESIAN_P_VALUES_SWITCH_FILE=".empty"
	MODEL_CHECKING_COMMENT_FILE=".empty"
	HSM_TEMPLATE_FIILENAME="hurdle-ordinal-latent-beta${DROP_SUFFIX}"
elif [ "${LIKELIHOOD}" == "ordinal-latent-normal" ]; then
	HSM_TEMPLATE_FIILENAME="ordinal-latent-normal"
fi
if [ "${PSM_TEMPLATE_FIILENAME}" == "TBD" ]; then
	PSM_TEMPLATE_FIILENAME="${HSM_TEMPLATE_FIILENAME}"
fi

# was: if [ "${VAR_TYPE}" == "hier" ]; then
SIGMA_TILDE_FIILENAME="${SLI}/.empty"
if [[ "${VAR_TYPE}" == "hier" && "${LIKELIHOOD}" != *"poisson"* ]]; then
	SIGMA_TILDE_FIILENAME="${SLI}/sigma-tilde-hier"
elif [[ "${VAR_TYPE}" == "fixed" && \
			"${VAR_LEVEL}" == "stratum" && \
			"${LIKELIHOOD}" == "ordinal-latent-normal" ]]; then
	SIGMA_TILDE_FIILENAME="${SLI}/sigma-tilde-fixed"
fi
# echo $SIGMA_TILDE_FIILENAME

if [ "${LIKELIHOOD}" == "gen-pois" ]; then
	HAT_QUANT_SUFFIX="mean.log.y"
	BAYESIAN_P_VALUES_SWITCH_FILE=".empty"
	HAT_QUANT_SUFFIX="eta"
	ADDITIONAL_QUANTS_PATH="${SLI}/additional-quantities"
	HAT_STRATUM_ADDITIONAL_QUANTS_FILE="${ADDITIONAL_QUANTS_PATH}/hat-stratum/${LIKELIHOOD}"
	PRED_STRATUM_ADDITIONAL_QUANTS_FILE="${ADDITIONAL_QUANTS_PATH}/pred-stratum/${LIKELIHOOD}"
fi

# Although the following models have quasi site-level hierarchical variances,
# there is not a separate group-level variance parameter
if [ "${LIKELIHOOD}" == "poisson" ] || \
   [ "${LIKELIHOOD}" == "gen-pois" ] || \
	 [ "${LIKELIHOOD}" == "negative-binomial-simple" ] || \
   [ "${LIKELIHOOD}" == "binomial" ] || \
   [ "${LIKELIHOOD}" == "zero-inflated-binomial" ]; then
	 SIGMA_TILDE_FIILENAME="${SLI}/.empty"
fi

TREND_STRAT_FILE=strat${FINITE_POP_SUFFIX}
if [ "${DERIVED_QUANTITIES}" == "no-trend" ]; then
	TREND_STRAT_FILE=.empty
	HAT_PARK_FILE=.empty # mostly a hack for MISC results
fi



MOD_FIILENAME=model.jags  #"${1}_${2}_${3}${GLE_ZEROS_JAGSFILE_SUFFIX}_${X}.jags"
MOD_FILE="${MOD_DIR}/${MOD_FIILENAME}"
# If a JAGS file already exists, delete it.
rm -f "${MOD_FILE}"

HSM_TEMPLATE_FILE="${GEN_QUANT_BLOCK}/hat-site-means/${HSM_TEMPLATE_FIILENAME}"
PSM_TEMPLATE_FILE="${GEN_QUANT_BLOCK}/pred-site-means/${PSM_TEMPLATE_FIILENAME}"
GP_TEMPLATE_SUFFIX=""
if [ "${PARAMETERIZATION}" == "non-centered" ]; then
	GP_TEMPLATE_SUFFIX="-nc"
fi
GP_TEMPLATE_FILE="${PRIORS_BLOCK}/group-level-effects/_${GROUP_LEVEL_EFFECTS}-template${GP_TEMPLATE_SUFFIX}${DROP_SUFFIX}_${B_DROP_ADDENDUM}"


STRATUM_NEW_OBS_SWITCH_FILE="${SLI}/new-obs/${LIKELIHOOD}"
if [ "${LIKELIHOOD}" == "ordinal-latent-normal" ]; then
	STRATUM_NEW_OBS_SWITCH_FILE="${STRATUM_NEW_OBS_SWITCH_FILE}${FINITE_POP_SUFFIX}"  # TODO: move into appropriate block above
fi
# PRED_SITE_OOS_SWITCH_FILE="${GEN_QUANT_BLOCK}/pred-site-means/pred-gles/${GROUP_LEVEL_EFFECTS}"  # TODO: not sure this exists!
PRED_SITE_OOS_SWITCH_FILE="${GEN_QUANT_BLOCK}/pred-site-means/.empty"
if [ "${LIKELIHOOD}" == "hurdle-ordinal-latent-beta" ]; then
	HAT_STRATUM_OOS_SWITCH_FILE="${SLI}/hat-gles-mix/${GROUP_LEVEL_EFFECTS}${GLE_ZEROS_SUFFIX}${DROP_SUFFIX}"
	PRED_SITE_OOS_SWITCH_FILE="${GEN_QUANT_BLOCK}/pred-site-means/pred-gles-mix/${GROUP_LEVEL_EFFECTS}${GLE_ZEROS_SUFFIX}${DROP_SUFFIX}"
fi
SITE_Z_SWITCH_FILE="${GEN_QUANT_BLOCK}/hat-site-means/.empty"
STRATUM_Z_SWITCH_FILE="${SLI}/.empty"
if [ "${LIKELIHOOD}" == "zero-inflated-binomial" ] || \
	[ "${LIKELIHOOD}" == "zero-inflated-beta-binomial" ]; then
		SITE_Z_SWITCH_FILE="${GEN_QUANT_BLOCK}/hat-site-means/site-z"
		STRATUM_Z_SWITCH_FILE="${SLI}/strat-z"
		HAT_STRATUM_OOS_SWITCH_FILE="${SLI}/hat-gles-zi/${GROUP_LEVEL_EFFECTS}"
		PRED_SITE_OOS_SWITCH_FILE="${GEN_QUANT_BLOCK}/pred-site-means/pred-gles-zi/${GROUP_LEVEL_EFFECTS}"
fi

SITE_EPS_SWITCH_FILE=".empty"
BB_OOS_SWITCH_FILE="${SLI}/new-obs/.empty"  # TODO: when pred.oos is working, remove!
HAT_STRATUM_EXTRA_VARIANCE_SWITCH_FILE="${SLI}/.empty"
PRED_STRATUM_EXTRA_VARIANCE_SWITCH_FILE="${SLI}/.empty"
if [ "${LIKELIHOOD}" == "beta-binomial" ] || \
	 [ "${LIKELIHOOD}" == "zero-inflated-beta-binomial" ]; then
		 HAT_STRATUM_EXTRA_VARIANCE_SWITCH_FILE="${SLI}/extra-variance/hat-dispersed-binomial${FINITE_POP_SUFFIX}"
		 PRED_STRATUM_EXTRA_VARIANCE_SWITCH_FILE="${SLI}/extra-variance/pred-dispersed-binomial${FINITE_POP_SUFFIX}"
		 HSM_TEMPLATE_FILE="${GEN_QUANT_BLOCK}/hat-site-means/${LIKELIHOOD}"
		 if [ "${DESIGN_MATRIX}" != "null" ]; then
		 		PSM_TEMPLATE_FILE="${GEN_QUANT_BLOCK}/pred-site-means/${LIKELIHOOD}"
		 fi
		 SITE_EPS_SWITCH_FILE="site-eps"
		 if [ "${FINITE_POP}" != "yes" ]; then
		 		BB_OOS_SWITCH_FILE="${SLI}/new-obs/bb-oos-${GROUP_LEVEL_EFFECTS}"
		 fi
fi

# if [ "${LIKELIHOOD}" == "gen-pois" ] || \
# 	 [ "${LIKELIHOOD}" == "lognormal" ] || \
#  	 [ "${LIKELIHOOD}" == "negative-binomial" ]; then
if [ "${LIKELIHOOD}" == "gen-pois" ] || \
	 [ "${LIKELIHOOD}" == "lognormal" ]; then
		 HSM_TEMPLATE_FILE="${GEN_QUANT_BLOCK}/hat-site-means/${LIKELIHOOD}"
		 if [ "${DESIGN_MATRIX}" != "null" ]; then
			 PSM_TEMPLATE_FILE="${GEN_QUANT_BLOCK}/pred-site-means/${LIKELIHOOD}"
		 fi
fi

# Override for variables defined above (affects gen-pois mods, for example).
if [ "${KEEP_PRED_SWITCHES}" == 'drop_pred_switches' ]; then
	PRED_STRATUM_ADDITIONAL_QUANTS_FILE="${SLI}/.empty"
	PSM_TEMPLATE_FILE="${GEN_QUANT_BLOCK}/pred-site-means/.empty"
fi

# Last tagged release.
RELEASE=$(cd model-api; git describe)

# Replace header in general model file with metadata for a given invocation.
awk '/# == top/ {f=1} !f; \
	 /# == bottom/ {print m1; print m2; print m3; print m4; print m5; print m6; f=0}' \
	 m1="# Likelihood: ${LIKELIHOOD}" \
	 m2="# Design matrix for additional covariates, X: ${X_DESC}" \
	 m3="# Deterministic model: ${DETERMINISTIC_MODEL}" \
	 m4="# Group-level effects: ${GROUP_LEVEL_EFFECTS}${GLE_ZEROS_JAGSFILE_SUFFIX}" \
	 m5="# Variance structure: ${VAR_TYPE}-${VAR_LEVEL}" \
	 m6="# Release: ${RELEASE:0:6}" \
	 model-api/src/model-builder/status-and-trend-general.jags > "${MOD_FILE}"
cp "${HSM_TEMPLATE_FILE}" "${HAT_SITE_MEANS_TEMPFILE}"
cp "${PSM_TEMPLATE_FILE}" "${PRED_SITE_MEANS_TEMPFILE}"
cp "${GP_TEMPLATE_FILE}" "${GLE_PRIORS_TEMPFILE}"
cp "${ACP_TEMPLATE_FILE}" "${ADDED_COV_PRIORS_TEMPFILE}"

# Function to consolidate the sed commands.
TEMPFILE=$(mktemp ${MOD_DIR}/tmp.XXXXXXXXX)
model_sub () {
    sed "$1" "${MOD_FILE}" > "${TEMPFILE}" \
		&& mv "${TEMPFILE}" "${MOD_FILE}"
}

template_sub () {
	TEMPLATE_TEMPFILE=$(mktemp ${MOD_DIR}/tmp.XXXXXXXXX)
    sed "$1" "$2" > "${TEMPLATE_TEMPFILE}" && mv "${TEMPLATE_TEMPFILE}" "$2"
}

# ==== derived quantities override ============================================
DERIVED_QUANTITIES_FILE=model-api/src/model-builder/status-and-trend-derived-quantities.jags
if [ "${DERIVED_QUANTITIES}" == "none" ]; then
	DERIVED_QUANTITIES_FILE=model-api/src/model-builder/.empty
fi
model_sub "/COMPLETE_SET_OF_DERIVED_QUANTITIES/{
	r ${DERIVED_QUANTITIES_FILE}
	/COMPLETE_SET_OF_DERIVED_QUANTITIES/d
}"

# ==== priors block substitution ==============================================
# Priors for the parameter(s) of the probability distribution.
model_sub "/DISTRIBUTION_PARAMS_SWITCH/{
	r ${PRIORS_BLOCK}/${DISTRIBUTION_PARAMS_SWITCH_FILE}${GLE_ZEROS_SUFFIX}${DROP_SUFFIX}${G_DROP_ADDENDUM}
	/DISTRIBUTION_PARAMS_SWITCH/d
}"
# Prior for `maxR` (if an asympototic function is umodel_sub).
model_sub "/MAX_R_SWITCH/{
	r ${PRIORS_BLOCK}/${MAX_R_SWITCH_FILE}
	/MAX_R_SWITCH/d
}"
# Priors for the parameters associated with each additional covariate.
# template_sub "s/VAGUE_COEF_PRIOR_INTERCEPT/${COEF_PRIOR_INTERCEPT}/g" "${ADDED_COV_PRIORS_TEMPFILE}"
template_sub "s/VAGUE_COEF_PRIOR_SLOPE/${COEF_PRIOR_SLOPE}/g" "${ADDED_COV_PRIORS_TEMPFILE}"
model_sub "/COVARIATE_PARAMS_SWITCH/{
	r ${ADDED_COV_PRIORS_TEMPFILE}
	/COVARIATE_PARAMS_SWITCH/d
}"
# Priors for any group-level effects.
template_sub "s/VAGUE_COEF_PRIOR_INTERCEPT/${COEF_PRIOR_INTERCEPT}/g" "${GLE_PRIORS_TEMPFILE}"
template_sub "s/VAGUE_COEF_PRIOR_SLOPE/${COEF_PRIOR_SLOPE}/g" "${GLE_PRIORS_TEMPFILE}"
template_sub "s/VAGUE_VAR_PRIOR_UL/${VAR_PRIOR_UL}/g" "${GLE_PRIORS_TEMPFILE}"
model_sub "/GROUP_LEVEL_EFFECTS_SWITCH/{
	r ${GLE_PRIORS_TEMPFILE}
	/GROUP_LEVEL_EFFECTS_SWITCH/d
}"

# ==== likelihood block substitution ==========================================
# Deterministic model.
model_sub "/RANDOM_BETA_SWITCH/{
	r ${CODE_BLOCKS}/deterministic-model/${RANDOM_BETA_DEF}
	/RANDOM_BETA_SWITCH/d
}"
model_sub "/DETERMINISTIC_MOD_SWITCH/{
	r ${CODE_BLOCKS}/deterministic-model/${X}/${DET_FUN_PREFIX}${DETERMINISTIC_MODEL}${B_SUFFIX}
	/DETERMINISTIC_MOD_SWITCH/d
}"
model_sub "s/OFFSET/${EXPOSURE}/"
# Likelihood.
model_sub "/LIKELIHOOD_SWITCH/{
	r ${CODE_BLOCKS}/likelihood/${LIKELIHOOD}
	/LIKELIHOOD_SWITCH/d
}"
model_sub "s/Y_SWITCH/${Y_SWITCH}/"
model_sub "/L_FOR_PA_HURDLE_SWITCH/{
	r ${HIL_TEMPLATE_FILE}
	/L_FOR_PA_HURDLE_SWITCH/d
}"
model_sub "/PHI_SWITCH/{
	r ${PHI_DEF}
	/PHI_SWITCH/d
}"
model_sub "/CENSOR_SWITCH_LIKELIHOOD/{
	r ${CENSOR_SWITCH_LIKELIHOOD_FILE}
	/CENSOR_SWITCH_LIKELIHOOD/d
}"
template_sub "s/TRUNCATION_SWITCH/${TRUNCATION_STRING}/" "${MOD_FILE}"
# Residuals.
model_sub "/RESIDUALS_SWITCH/{
	r ${GEN_QUANT_BLOCK}/epsilon/${EPSILON_FILE}
	/RESIDUALS_SWITCH/d
}"
model_sub "/EXTRA_VARIANCE_SWITCH/{
	r ${CODE_BLOCKS}/deterministic-model/${X}/extra-variance/${DETERMINISTIC_MODEL}
	/EXTRA_VARIANCE_SWITCH/d
}"
model_sub "/DEFLECTIONS_SWITCH/{
	r ${DEFLECTIONS_FILE}
	/DEFLECTIONS_SWITCH/d
}"
model_sub "/MU_D_BLOCK/{
	r ${MU_D_BLOCK}
	/MU_D_BLOCK/d
}"

# ==== derived quantities block substitution ==================================
# Inference.
template_sub "s/LINK/${LINK}/" "${HAT_SITE_MEANS_TEMPFILE}"
model_sub "/HAT_SITE_MEAN_SWITCH/{
	r ${HAT_SITE_MEANS_TEMPFILE}
	/HAT_SITE_MEAN_SWITCH/d
}"
model_sub "/HAT_Y_SIM_SWITCH/{
	r ${GEN_QUANT_BLOCK}/y-sim/${LIKELIHOOD}
	/HAT_Y_SIM_SWITCH/d
}"
template_sub "s/NEW_OBS_PREFIX/hat/g" "${MOD_FILE}"
template_sub "s/INDICES/j, t, k/g" "${MOD_FILE}"
# template_sub "s/SECIDNI/j, k, t/g" "${MOD_FILE}"
template_sub "s/K_IND/k/g" "${MOD_FILE}"
template_sub "s/SL_IND/j, k/g" "${MOD_FILE}"
template_sub "s/TABHOLDER /				/g" "${MOD_FILE}"
model_sub "/CHANGE_IN_ODDS_SWITCH/{
	r ${GEN_QUANT_BLOCK}/change-in-odds/${CHANGE_IN_ODDS}
	/CHANGE_IN_ODDS_SWITCH/d
}"
model_sub "/MISC_PARK_LEVEL_INFERENCE_SWITCH/{
	r ${GEN_QUANT_BLOCK}/park-level-inference/misc/${MISC_INFERENCE_FILE}
	/MISC_PARK_LEVEL_INFERENCE_SWITCH/d
}"
model_sub "/HAT_PARK_SWITCH/{
	r ${GEN_QUANT_BLOCK}/park-level-inference/${HAT_PARK_FILE}${FINITE_POP_SUFFIX}
	/HAT_PARK_SWITCH/d
}"
model_sub "/HAT_STRATUM_SWITCH/{
	r ${SLI}/${HAT_STRATUM_SWITCH_FILE}
	/HAT_STRATUM_SWITCH/d
}"
model_sub "/PRED_STRATUM_SWITCH/{
	r ${SLI}/${PRED_STRATUM_SWITCH_FILE}
	/PRED_STRATUM_SWITCH/d
}"
model_sub "/PRED_STRATUM_NEW_OBS_SWITCH/{
	r ${STRATUM_NEW_OBS_SWITCH_FILE}
	/PRED_STRATUM_NEW_OBS_SWITCH/d
}"
template_sub "s/NEW_OBS_PREFIX_TOUPPER/PRED/g" "${MOD_FILE}"
template_sub "s/NEW_OBS_PREFIX/pred/g" "${MOD_FILE}"
model_sub "/BB_OOS_SWITCH/{
	r ${SLI}/new-obs/.empty
	/BB_OOS_SWITCH/d
}"
model_sub "/OLN_OOS_SWITCH/{
	r ${SLI}/new-obs/.empty
	/OLN_OOS_SWITCH/d
}"
model_sub "/HAT_STRATUM_OOS_SWITCH/{
	r ${HAT_STRATUM_OOS_SWITCH_FILE}
	/HAT_STRATUM_OOS_SWITCH/d
}"
model_sub "/BETA_TILDE_DRAW/{
	r ${SLI}/${GROUP_LEVEL_EFFECTS}-draw${GP_TEMPLATE_SUFFIX}${B_DROP_ADDENDUM}
	/BETA_TILDE_DRAW/d
}"
model_sub "/GAMMA_TILDE_DRAW/{
	r ${SLI}/${GLE_ZEROS}-draw${G_DROP_ADDENDUM}
	/GAMMA_TILDE_DRAW/d
}"
model_sub "/HAT_STRATUM_NEW_OBS_SWITCH/{
	r ${STRATUM_NEW_OBS_SWITCH_FILE}
	/HAT_STRATUM_NEW_OBS_SWITCH/d
}"
template_sub "s/NEW_OBS_PREFIX_TOUPPER/HAT/g" "${MOD_FILE}"
template_sub "s/NEW_OBS_PREFIX/hat/g" "${MOD_FILE}"
model_sub "/BB_OOS_SWITCH/{
	r ${BB_OOS_SWITCH_FILE}
	/BB_OOS_SWITCH/d
}"
model_sub "/OLN_OOS_SWITCH/{
	r ${SLI}/new-obs/oln-oos
	/OLN_OOS_SWITCH/d
}"
model_sub "/PRED_PARK_SWITCH/{
	r ${GEN_QUANT_BLOCK}/park-level-inference/${PRED_PARK_SWITCH_FILE}${BONUS_SUFFIX}
	/PRED_PARK_SWITCH/d
}"
model_sub "/TREND_STRAT_SWITCH/{
	r ${GEN_QUANT_BLOCK}/trend/${TREND_STRAT_FILE}
	/TREND_STRAT_SWITCH/d
}"
model_sub "/TREND_ZONE_SWITCH/{
	r ${GEN_QUANT_BLOCK}/trend/${TREND_ZONE_FILE}
	/TREND_ZONE_SWITCH/d
}"
model_sub "/GLE_DEP_PARK_LEVEL_INFERENCE_SWITCH/{
	r ${GEN_QUANT_BLOCK}/park-level-inference/${GROUP_LEVEL_EFFECTS}
	/GLE_DEP_PARK_LEVEL_INFERENCE_SWITCH/d
}"
template_sub "s/LINK/${LINK}/" "${PRED_SITE_MEANS_TEMPFILE}"
model_sub "/PRED_SITE_MEAN_SWITCH/{
	r ${PRED_SITE_MEANS_TEMPFILE}
	/PRED_SITE_MEAN_SWITCH/d
}"
model_sub "/PRED_SITE_NEW_OBS_SWITCH/{
	r ${GEN_QUANT_BLOCK}/y-sim/${LIKELIHOOD}
	/PRED_SITE_NEW_OBS_SWITCH/d
}"
model_sub "/PRED_SITE_OOS_SWITCH/{
	r ${PRED_SITE_OOS_SWITCH_FILE}
	/PRED_SITE_OOS_SWITCH/d
}"
template_sub "s/NEW_OBS_PREFIX/pred/g" "${MOD_FILE}"
template_sub "s/INDICES/j.pred[i], x.pred.index[i], k.pred[i]/g" "${MOD_FILE}"
# template_sub "s/SECIDNI/j.pred[i], k.pred[i], x.pred.index[i]/g" "${MOD_FILE}"
template_sub "s/K_IND/k.pred[i]/g" "${MOD_FILE}"
template_sub "s/SL_IND/j.pred[i], k.pred[i]/g" "${MOD_FILE}"
template_sub "s/TABHOLDER //g" "${MOD_FILE}"
model_sub "/TAU_SITE_SWITCH/{
	r ${GEN_QUANT_BLOCK}/hat-site-means/${TAU_SITE_FILE}
	/TAU_SITE_SWITCH/d
}"
model_sub "/SITE_EPS_SWITCH/{
	r ${GEN_QUANT_BLOCK}/epsilon/${SITE_EPS_SWITCH_FILE}
	/SITE_EPS_SWITCH/d
}"
model_sub "/HAT_STRAT_DEFL_MEAN/{
	r ${HAT_STRAT_DEFL_MEAN}
	/HAT_STRAT_DEFL_MEAN/d
}"
model_sub "/HAT_PNORM_SWITCH/{
	r ${HAT_PNORM}
	/HAT_PNORM_SWITCH/d
}"
model_sub "/PRED_PNORM_SWITCH/{
	r ${PRED_PNORM}
	/PRED_PNORM_SWITCH/d
}"
model_sub "/PRED_STRAT_DEFL_MEAN/{
	r ${PRED_STRAT_DEFL_MEAN}
	/PRED_STRAT_DEFL_MEAN/d
}"
template_sub "s/_B0_4_SITE_FROM_ZONE_/${_B0_4_SITE_FROM_ZONE_}/g" "${MOD_FILE}"
template_sub "s/_B1_4_SITE_FROM_ZONE_/${_B1_4_SITE_FROM_ZONE_}/g" "${MOD_FILE}"
template_sub "s/PRED_SITE_IDX/${PRED_SITE_IDX}/" "${MOD_FILE}"

# Model checking.
model_sub "/MODEL_CHECKING_COMMENT/{
	r ${MOD_CHECK_BLOCK}/${MODEL_CHECKING_COMMENT_FILE}
	/MODEL_CHECKING_COMMENT/d
}"
model_sub "/MODEL_CHECKING_SWITCH/{
	r ${MOD_CHECK_BLOCK}/${LIKELIHOOD}
	/MODEL_CHECKING_SWITCH/d
}"
template_sub "s/TRUNCATION_SWITCH/${TRUNCATION_STRING}/" "${MOD_FILE}"
model_sub "/BAYESIAN_P_VALUES_SWITCH/{
	r ${MOD_CHECK_BLOCK}/${BAYESIAN_P_VALUES_SWITCH_FILE}
	/BAYESIAN_P_VALUES_SWITCH/d
}"
template_sub "s/Y_VEC_SWITCH/${Y_VEC}/" "${MOD_FILE}"
model_sub "/CENSOR_SWITCH_CHECKING/{
	r ${CENSOR_SWITCH_CHECKING_FILE}
	/CENSOR_SWITCH_CHECKING/d
}"

# Extra variance.
model_sub "/HAT_STRATUM_EXTRA_VARIANCE_SWITCH/{
	r ${HAT_STRATUM_EXTRA_VARIANCE_SWITCH_FILE}
	/HAT_STRATUM_EXTRA_VARIANCE_SWITCH/d
}"
model_sub "/PRED_STRATUM_EXTRA_VARIANCE_SWITCH/{
	r ${PRED_STRATUM_EXTRA_VARIANCE_SWITCH_FILE}
	/PRED_STRATUM_EXTRA_VARIANCE_SWITCH/d
}"
model_sub "/SITE_Z_SWITCH/{
	r ${SITE_Z_SWITCH_FILE}
	/SITE_Z_SWITCH/d
}"
model_sub "/STRATUM_Z_SWITCH/{
	r ${STRATUM_Z_SWITCH_FILE}
	/STRATUM_Z_SWITCH/d
}"
model_sub "/HAT_ZI_ZONE_QUANTS/{
	r ${HAT_ZI_ZONE_QUANTS}
	/HAT_ZI_ZONE_QUANTS/d
}"
model_sub "/PRED_ZI_ZONE_QUANTS/{
	r ${PRED_ZI_ZONE_QUANTS}
	/PRED_ZI_ZONE_QUANTS/d
}"

model_sub "/DRIVER_EVAL_SWITCH/{
	r ${SLI}/${DRIVER_MARGINAL_EFFECTS}
	/DRIVER_EVAL_SWITCH/d
}"

# Final LINK replacements.
template_sub "s/LINK/${LINK}/" "${MOD_FILE}"
template_sub "s/HAT_QUANT_SUFFIX/${HAT_QUANT_SUFFIX}/g" "${MOD_FILE}"

# Variances structure.
#VAR_LEVEL_PATH="${PRIORS_BLOCK}/variances/${VAR_LEVEL}"
model_sub "/VARIANCES_SPEC_SWITCH/{
	r ${PRIORS_BLOCK}/variances/${VAR_LEVEL}-level${VAR_MOD_KIND}
	/VARIANCES_SPEC_SWITCH/d
}"
model_sub "/VARIANCES_SPEC_HYPERS_SWITCH/{
	r ${PRIORS_BLOCK}/variances/${VAR_LEVEL}-level-hypers${VAR_MOD_KIND}
	/VARIANCES_SPEC_HYPERS_SWITCH/d
}"
template_sub "s/MU_G_SHAPE/${MU_G_SHAPE}/" "${MOD_FILE}"
template_sub "s/MU_G_RATE/${MU_G_RATE}/" "${MOD_FILE}"
template_sub "s/SIGMA_G_SHAPE/${SIGMA_G_SHAPE}/" "${MOD_FILE}"
template_sub "s/SIGMA_G_RATE/${SIGMA_G_RATE}/" "${MOD_FILE}"
template_sub "s/MU_SIGMA_PRIOR_DIST/${MU_SIGMA_PRIOR_DIST}/" "${MOD_FILE}"
template_sub "s/SIGMA_SIGMA_PRIOR_DIST/${SIGMA_SIGMA_PRIOR_DIST}/" "${MOD_FILE}"
template_sub "s/SIGMA_K_UPPER/${SIGMA_K_UPPER}/" "${MOD_FILE}"
template_sub "s/SIGMA_B0_PRIOR/${SIGMA_B0_PRIOR}/" "${MOD_FILE}"
template_sub "s/SIGMA_B1_PRIOR/${SIGMA_B1_PRIOR}/" "${MOD_FILE}"
model_sub "/HAT_EPS_TILDE/{
	r ${SLI}/extra-variance/site-eps-tilde/hat-${VAR_TYPE}-${VAR_LEVEL}-var
	/HAT_EPS_TILDE/d
}"
model_sub "/PRED_EPS_TILDE/{
	r ${SLI}/extra-variance/site-eps-tilde/pred-${VAR_TYPE}-${VAR_LEVEL}-var
	/PRED_EPS_TILDE/d
}"

# Any additional quantities.
model_sub "/HAT_STRATUM_ADDITIONAL_QUANTS/{
	r ${HAT_STRATUM_ADDITIONAL_QUANTS_FILE}
	/HAT_STRATUM_ADDITIONAL_QUANTS/d
}"
model_sub "/PRED_STRATUM_ADDITIONAL_QUANTS/{
	r ${PRED_STRATUM_ADDITIONAL_QUANTS_FILE}
	/PRED_STRATUM_ADDITIONAL_QUANTS/d
}"
template_sub "s/HAT_STRAT_LIN_PRED/${HAT_STRAT_LIN_PRED}/g" "${MOD_FILE}"
template_sub "s/PRED_STRAT_LIN_PRED/${PRED_STRAT_LIN_PRED}/g" "${MOD_FILE}"
template_sub "s/DRIVER_STRAT_LIN_PRED/${DRIVER_STRAT_LIN_PRED}/g" "${MOD_FILE}"

# model_sub "/MM_DETAILS/{
# 	r ${MM_DETAILS_FILE}
# 	/MM_DETAILS/d
# }"
model_sub "/SIGMA_TILDE_DRAW/{
	r ${SIGMA_TILDE_FIILENAME}
	/SIGMA_TILDE_DRAW/d
}"

template_sub "s/WHICH_BETA/${THIS_BETA}/g" "${MOD_FILE}"
# template_sub "s/VAGUE_COEF_PRIOR_SLOPE/${COEF_PRIOR_SLOPE}/g" "${MOD_FILE}"

# template_sub "s/MM_DETAILS/${MM_DETAILS}/" "${MOD_FILE}"
# model_sub "/LIKELIHOOD_RELATED_DQS_SWITCH/{
# 	r ${LIKELIHOOD_RELATED_DQS_SWITCH_FILE}
# 	/LIKELIHOOD_RELATED_DQS_SWITCH/d
# }"
# model_sub "/GAP_SIZE_THRESH_SWITCH/{
# 	r ${GAP_SIZE_THRESH_SWITCH_FILE}
# 	/GAP_SIZE_THRESH_SWITCH/d
# }"

# # TODO: add the following.
model_sub "/HAT_SITE_VARIANCE_SWITCH/{
	r ${HAT_SITE_VARIANCE}
	/HAT_SITE_VARIANCE_SWITCH/d
}"
model_sub "/PRED_SITE_VARIANCE_SWITCH/{
	r ${HAT_SITE_VARIANCE}
	/PRED_SITE_VARIANCE_SWITCH/d
}"
model_sub "/HAT_STRATUM_VARIANCE_SWITCH/{
	r ${HAT_STRATUM_VARIANCE}
	/HAT_STRATUM_VARIANCE_SWITCH/d
}"
model_sub "/PRED_STRATUM_VARIANCE_SWITCH/{
	r ${HAT_STRATUM_VARIANCE}
	/PRED_STRATUM_VARIANCE_SWITCH/d
}"
model_sub "s/OFFSET/${EXPOSURE}/" # possibly redundant with the call above

rm "${HAT_SITE_MEANS_TEMPFILE}"
rm "${PRED_SITE_MEANS_TEMPFILE}"
rm "${GLE_PRIORS_TEMPFILE}"
rm "${ADDED_COV_PRIORS_TEMPFILE}"
echo "${MOD_FILE}" > "${MOD_DIR}/crumbs.txt"
