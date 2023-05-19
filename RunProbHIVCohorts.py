import deampy.plots.histogram as hist
import deampy.plots.sample_paths as path

import hiv_model_econ_eval.input_data as data
import hiv_model_param_uncertainty.model_classes as model
import hiv_model_param_uncertainty.param_classes as param
import hiv_model_param_uncertainty.support as support

POP_SIZE = 1000             # cohort population size
N_COHORTS = 50              # number of cohorts
therapy = param.Therapies.MONO  # selected therapy

# create multiple cohort
multiCohort = model.MultiCohort(
    ids=range(N_COHORTS),
    pop_size=POP_SIZE,
    therapy=therapy)

multiCohort.simulate(n_time_steps=data.SIM_TIME_STEPS)

# plot the sample paths
path.plot_sample_paths(
    sample_paths=multiCohort.multiCohortOutcomes.survivalCurves,
    title='Survival Curves',
    x_label='Time-Step (Year)',
    y_label='Number Survived',
    transparency=0.5)

# plot the histogram of average survival time
hist.plot_histogram(
    data=multiCohort.multiCohortOutcomes.meanSurvivalTimes,
    title='Histograms of Mean Survival Time',
    x_label='Survival Time (year)',
    bin_width=0.5,
    x_range=[5, 20])

# print the outcomes of this simulated cohort
support.print_outcomes(multi_cohort_outcomes=multiCohort.multiCohortOutcomes,
                       therapy_name=therapy)
