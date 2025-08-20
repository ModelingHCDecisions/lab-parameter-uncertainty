import hiv_model_econ_eval.input_data as data
import hiv_model_param_uncertainty.model_classes as model
import hiv_model_param_uncertainty.param_classes as param
import hiv_model_param_uncertainty.support as support

N_COHORTS = 100  # number of cohorts
POP_SIZE = 50  # population size of each cohort

# create a multi-cohort to simulate under mono therapy
multiCohortMono = model.MultiCohort(
    ids=range(N_COHORTS),
    pop_size=POP_SIZE,
    therapy=param.Therapies.MONO
)

multiCohortMono.simulate(n_time_steps=data.SIM_TIME_STEPS)

# create a multi-cohort to simulate under combo therapy
multiCohortCombo = model.MultiCohort(
    ids=range(N_COHORTS),
    pop_size=POP_SIZE,
    therapy=param.Therapies.COMBO
)

multiCohortCombo.simulate(n_time_steps=data.SIM_TIME_STEPS)

# print the estimates for the mean survival time and mean time to AIDS
support.print_outcomes(multi_cohort_outcomes=multiCohortMono.multiCohortOutcomes,
                       therapy_name=param.Therapies.MONO)
support.print_outcomes(multi_cohort_outcomes=multiCohortCombo.multiCohortOutcomes,
                       therapy_name=param.Therapies.COMBO)

# draw survival curves and histograms
support.plot_survival_curves_and_histograms(multi_cohort_outcomes_mono=multiCohortMono.multiCohortOutcomes,
                                            multi_cohort_outcomes_combo=multiCohortCombo.multiCohortOutcomes)

# print comparative outcomes
support.print_comparative_outcomes(multi_cohort_outcomes_mono=multiCohortMono.multiCohortOutcomes,
                                   multi_cohort_outcomes_combo=multiCohortCombo.multiCohortOutcomes)

# report the CEA results
support.report_CEA_CBA(multi_cohort_outcomes_mono=multiCohortMono.multiCohortOutcomes,
                       multi_cohort_outcomes_combo=multiCohortCombo.multiCohortOutcomes)