import deampy.statistics as stat

from hiv_model_econ_eval.model_classes import Cohort
from hiv_model_param_uncertainty.param_classes import ParameterGenerator


class MultiCohort:
    """ simulates multiple cohorts with different parameters """

    def __init__(self, ids, pop_size, therapy):
        """
        :param ids: (list) of ids for cohorts to simulate
        :param pop_size: (int) population size of cohorts to simulate
        :param therapy: selected therapy
        """
        self.ids = ids
        self.popSize = pop_size
        self.therapy = therapy
        self.paramSets = []  # list of parameter sets each of which corresponds to a cohort
        self.multiCohortOutcomes = MultiCohortOutcomes()
        self.paramGenerator = ParameterGenerator(therapy=self.therapy)

    def simulate(self, n_time_steps):
        """ simulates all cohorts
        :param n_time_steps: number of simulation time steps
        """

        for i in range(len(self.ids)):

            # get a new set of parameter values
            param_set = self.paramGenerator.get_new_parameters(seed=i)

            # create a cohort
            cohort = Cohort(id=self.ids[i],
                            pop_size=self.popSize,
                            parameters=param_set)

            # simulate the cohort
            cohort.simulate(n_time_steps=n_time_steps)

            # extract the outcomes of this simulated cohort
            self.multiCohortOutcomes.extract_outcomes(simulated_cohort=cohort)

        # calculate the summary statistics of outcomes from all cohorts
        self.multiCohortOutcomes.calculate_summary_stats()


class MultiCohortOutcomes:
    def __init__(self):

        self.survivalCurves = []  # list of survival curves from all simulated cohorts

        self.meanSurvivalTimes = []  # list of average patient survival time from each simulated cohort
        self.meanTimeToAIDS = []     # list of average patient time until AIDS from each simulated cohort
        self.meanCosts = []          # list of average patient cost from each simulated cohort
        self.meanQALYs = []          # list of average patient QALY from each simulated cohort

        self.statMeanSurvivalTime = None    # summary statistics of average survival time
        self.statMeanTimeToAIDS = None      # summary statistics of average time until AIDS
        self.statMeanCost = None            # summary statistics of average cost
        self.statMeanQALY = None            # summary statistics of average QALY

    def extract_outcomes(self, simulated_cohort):
        """ extracts outcomes of a simulated cohort
        :param simulated_cohort: a cohort after being simulated"""

        # append the survival curve of this cohort
        self.survivalCurves.append(simulated_cohort.cohortOutcomes.nLivingPatients)

        # store mean survival time from this cohort
        self.meanSurvivalTimes.append(simulated_cohort.cohortOutcomes.statSurvivalTime.get_mean())
        # store mean time to AIDS from this cohort
        self.meanTimeToAIDS.append(simulated_cohort.cohortOutcomes.statTimeToAIDS.get_mean())
        # store mean cost from this cohort
        self.meanCosts.append(simulated_cohort.cohortOutcomes.statCost.get_mean())
        # store mean QALY from this cohort
        self.meanQALYs.append(simulated_cohort.cohortOutcomes.statUtility.get_mean())

    def calculate_summary_stats(self):
        """
        calculate the summary statistics
        """

        # summary statistics of mean survival time
        self.statMeanSurvivalTime = stat.SummaryStat(name='Average survival time',
                                                     data=self.meanSurvivalTimes)
        # summary statistics of mean time to AIDS
        self.statMeanTimeToAIDS = stat.SummaryStat(name='Average time to AIDS',
                                                   data=self.meanTimeToAIDS)
        # summary statistics of mean cost
        self.statMeanCost = stat.SummaryStat(name='Average cost',
                                             data=self.meanCosts)
        # summary statistics of mean QALY
        self.statMeanQALY = stat.SummaryStat(name='Average QALY',
                                             data=self.meanQALYs)
