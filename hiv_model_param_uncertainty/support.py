import deampy.econ_eval as econ
import deampy.plots.histogram as hist
import deampy.plots.sample_paths as path
import deampy.statistics as stat

import hiv_model_econ_eval.input_data as data


def print_outcomes(multi_cohort_outcomes, therapy_name):
    """ prints the outcomes of a simulated cohort
    :param multi_cohort_outcomes: outcomes of a simulated multi-cohort
    :param therapy_name: the name of the selected therapy
    """
    # mean and prediction interval of patient survival time
    survival_mean_PI_text = multi_cohort_outcomes.statMeanSurvivalTime.get_formatted_mean_and_interval(
        interval_type='p', alpha=data.ALPHA, deci=2)

    # mean and prediction interval text of time to AIDS
    time_to_HIV_death_PI_text = multi_cohort_outcomes.statMeanTimeToAIDS.get_formatted_mean_and_interval(
        interval_type='p', alpha=data.ALPHA, deci=2)

    # mean and prediction interval text of discounted total cost
    cost_mean_PI_text = multi_cohort_outcomes.statMeanCost.get_formatted_mean_and_interval(
        interval_type='p', alpha=data.ALPHA, deci=2, form=',')

    # mean and prediction interval text of discounted total QALY
    utility_mean_PI_text = multi_cohort_outcomes.statMeanQALY.get_formatted_mean_and_interval(
        interval_type='p', alpha=data.ALPHA, deci=2)

    # print outcomes
    print(therapy_name)
    print("  Estimate of mean survival time and {:.{prec}%} uncertainty interval:".format(1 - data.ALPHA, prec=0),
          survival_mean_PI_text)
    print("  Estimate of mean time to AIDS and {:.{prec}%} uncertainty interval:".format(1 - data.ALPHA, prec=0),
          time_to_HIV_death_PI_text)
    print("  Estimate of mean discounted cost and {:.{prec}%} uncertainty interval:".format(1 - data.ALPHA, prec=0),
          cost_mean_PI_text)
    print("  Estimate of mean discounted utility and {:.{prec}%} uncertainty interval:".format(1 - data.ALPHA, prec=0),
          utility_mean_PI_text)
    print("")


def plot_survival_curves_and_histograms(multi_cohort_outcomes_mono, multi_cohort_outcomes_combo):
    """ plot the survival curves and the histograms of survival times
    :param multi_cohort_outcomes_mono: outcomes of a multi-cohort simulated under mono therapy
    :param multi_cohort_outcomes_combo: outcomes of a multi-cohort simulated under combination therapy
    """

    # get survival curves of both treatments
    sets_of_survival_curves = [
        multi_cohort_outcomes_mono.survivalCurves,
        multi_cohort_outcomes_combo.survivalCurves
    ]

    # graph survival curve
    path.plot_sets_of_sample_paths(
        sets_of_sample_paths=sets_of_survival_curves,
        title='Survival Curves',
        x_label='Simulation Time Step (year)',
        y_label='Number of Patients Alive',
        legends=['Mono Therapy', 'Combination Therapy'],
        transparency=0.4,
        color_codes=['green', 'blue'],
        figure_size=(6, 5),
        file_name='figs/survival_curves.png'
    )

    # histograms of survival times
    set_of_survival_times = [
        multi_cohort_outcomes_mono.meanSurvivalTimes,
        multi_cohort_outcomes_combo.meanSurvivalTimes
    ]

    # graph histograms
    hist.plot_histograms(
        data_sets=set_of_survival_times,
        title='Histograms of Mean Survival Time',
        x_label='Survival Time (year)',
        y_label='Counts',
        bin_width=0.5,
        x_range=[5, 20],
        legends=['Mono Therapy', 'Combination Therapy'],
        color_codes=['green', 'blue'],
        transparency=0.5,
        figure_size=(6, 5),
        file_name='figs/survival_times.png'
    )


def print_comparative_outcomes(multi_cohort_outcomes_mono, multi_cohort_outcomes_combo):
    """ prints average increase in survival time, discounted cost, and discounted utility
    under combination therapy compared to mono therapy
    :param multi_cohort_outcomes_mono: outcomes of a multi-cohort simulated under mono therapy
    :param multi_cohort_outcomes_combo: outcomes of a multi-cohort simulated under combination therapy
    """

    # increase in mean survival time under combination therapy with respect to mono therapy
    increase_mean_survival_time = stat.DifferenceStatPaired(
        name='Increase in mean survival time',
        x=multi_cohort_outcomes_combo.meanSurvivalTimes,
        y_ref=multi_cohort_outcomes_mono.meanSurvivalTimes)

    # estimate and PI
    estimate_PI = increase_mean_survival_time.get_formatted_mean_and_interval(
        interval_type='p', alpha=data.ALPHA, deci=2)
    print("Increase in mean survival time and {:.{prec}%} uncertainty interval:"
          .format(1 - data.ALPHA, prec=0), estimate_PI)

    # increase in mean discounted cost under combination therapy with respect to mono therapy
    increase_mean_discounted_cost = stat.DifferenceStatPaired(
        name='Increase in mean discounted cost',
        x=multi_cohort_outcomes_combo.meanCosts,
        y_ref=multi_cohort_outcomes_mono.meanCosts)

    # estimate and PI
    estimate_PI = increase_mean_discounted_cost.get_formatted_mean_and_interval(
        interval_type='p', alpha=data.ALPHA, deci=2, form=',')
    print("Increase in mean discounted cost and {:.{prec}%} uncertainty interval:"
          .format(1 - data.ALPHA, prec=0), estimate_PI)

    # increase in mean discounted QALY under combination therapy with respect to mono therapy
    increase_mean_discounted_qaly = stat.DifferenceStatPaired(
        name='Increase in mean discounted QALY',
        x=multi_cohort_outcomes_combo.meanQALYs,
        y_ref=multi_cohort_outcomes_mono.meanQALYs)

    # estimate and PI
    estimate_PI = increase_mean_discounted_qaly.get_formatted_mean_and_interval(
        interval_type='p', alpha=data.ALPHA, deci=2)
    print("Increase in mean discounted utility and {:.{prec}%} uncertainty interval:"
          .format(1 - data.ALPHA, prec=0), estimate_PI)


def report_CEA_CBA(multi_cohort_outcomes_mono, multi_cohort_outcomes_combo):
    """ performs cost-effectiveness and cost-benefit analyses
    :param multi_cohort_outcomes_mono: outcomes of a multi-cohort simulated under mono therapy
    :param multi_cohort_outcomes_combo: outcomes of a multi-cohort simulated under combination therapy
    """

    # define two strategies
    mono_therapy_strategy = econ.Strategy(
        name='Mono Therapy',
        cost_obs=multi_cohort_outcomes_mono.meanCosts,
        effect_obs=multi_cohort_outcomes_mono.meanQALYs,
        color='green'
    )
    combo_therapy_strategy = econ.Strategy(
        name='Combination Therapy',
        cost_obs=multi_cohort_outcomes_combo.meanCosts,
        effect_obs=multi_cohort_outcomes_combo.meanQALYs,
        color='blue'
    )

    # do CEA
    CEA = econ.CEA(
        strategies=[mono_therapy_strategy, combo_therapy_strategy],
        if_paired=True
    )

    # show the cost-effectiveness plane
    CEA.plot_CE_plane(
        title='Cost-Effectiveness Analysis',
        x_label='Additional Discounted QALY',
        y_label='Additional Discounted Cost',
        fig_size=(6, 5),
        add_clouds=True,
        transparency=0.2,
        file_name='figs/cea.png')

    # report the CE table
    CEA.build_CE_table(
        interval_type='p',  # uncertainty (projection) interval for cost and effect estimates but
                            # for ICER, confidence interval will be reported.
        alpha=data.ALPHA,
        cost_digits=0,
        effect_digits=2,
        icer_digits=2,
        file_name='CETable.csv')

    # CBA
    NBA = econ.CBA(
        strategies=[mono_therapy_strategy, combo_therapy_strategy],
        wtp_range=(0, 50000),
        if_paired=True
    )
    # show the net monetary benefit figure
    NBA.plot_marginal_nmb_lines(
        title='Cost-Benefit Analysis',
        x_label='Willingness-To-Pay for One Additional QALY ($)',
        y_label='Incremental Net Monetary Benefit ($)',
        interval_type='c', # show confidence interval
        show_legend=True,
        figure_size=(6, 5),
        file_name='figs/nmb.png'
    )
