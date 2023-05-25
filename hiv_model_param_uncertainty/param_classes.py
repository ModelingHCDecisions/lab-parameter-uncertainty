import math

import deampy.random_variates as rvgs
import numpy as np
import scipy.stats as stat

import hiv_model_econ_eval.input_data as data
from hiv_model_econ_eval.param_classes import Therapies, get_prob_matrix_combo


class Parameters:
    """ class to include parameter information to simulate the model """

    def __init__(self, therapy):

        self.therapy = therapy              # selected therapy
        self.initialHealthState = data.HealthStates.CD4_200to500     # initial health state
        self.annualTreatmentCost = 0        # annual treatment cost
        self.probMatrix = []                # transition probability matrix of the selected therapy
        self.annualStateCosts = []          # annual state costs
        self.annualStateUtilities = []      # annual state utilities
        self.discountRate = data.DISCOUNT   # discount rate


class ParameterGenerator:
    """ class to generate parameter values from the selected probability distributions """

    def __init__(self, therapy):

        self.therapy = therapy
        self.probMatrixRVG = []     # list of dirichlet distributions for transition probabilities
        self.lnRelativeRiskRVG = None  # normal distribution for the natural log of the treatment relative risk
        self.annualStateCostRVGs = []  # list of gamma distributions for the annual cost of states
        self.annualStateUtilityRVGs = []  # list of beta distributions for the annual utility of states
        self.annualZidovudineCostRVG = None   # gamma distribution for the cost of zidovudine
        self.annualLamivudineCostRVG = None   # gamma distribution for the cost of lamivudine

        # create Dirichlet distributions for transition probabilities
        for probs in data.TRANS_MATRIX:
            # note:  for a Dirichlet distribution all values of the argument 'a' should be non-zero.
            # setting if_ignore_0s to True allows the Dirichlet distribution to take 'a' with zero values.
            self.probMatrixRVG.append(rvgs.Dirichlet(
                a=probs, if_ignore_0s=True))

        # treatment relative risk
        rr_ci = [0.365, 0.71]  # confidence interval of the treatment relative risk

        # find the mean and st_dev of the normal distribution assumed for ln(RR)
        # sample mean ln(RR)
        mean_ln_rr = math.log(data.TREATMENT_RR)
        # sample standard deviation of ln(RR)
        std_ln_rr = \
            (math.log(rr_ci[1]) - math.log(rr_ci[0])) / (2 * stat.norm.ppf(1 - 0.05 / 2))
        # create a normal distribution for ln(RR)
        self.lnRelativeRiskRVG = rvgs.Normal(loc=mean_ln_rr,
                                             scale=std_ln_rr)

        # create gamma distributions for annual state cost
        for cost in data.ANNUAL_STATE_COST:

            # if cost is zero, add a constant 0, otherwise add a gamma distribution
            if cost == 0:
                self.annualStateCostRVGs.append(rvgs.Constant(value=0))
            else:
                # find shape and scale of the assumed gamma distribution
                # no data available to estimate the standard deviation, so we assumed st_dev=cost / 5
                fit_output = rvgs.Gamma.fit_mm(mean=cost, st_dev=cost / 5)
                # append the distribution
                self.annualStateCostRVGs.append(
                    rvgs.Gamma(a=fit_output["a"],
                               loc=0,
                               scale=fit_output["scale"]))

        # create a gamma distribution for annual treatment cost with each drug
        # first fit the gamma distribution to the cost of each drug
        fit_output_zido = rvgs.Gamma.fit_mm(mean=data.Zidovudine_COST, st_dev=data.Zidovudine_COST / 5)
        fit_output_lami = rvgs.Gamma.fit_mm(mean=data.Lamivudine_COST, st_dev=data.Lamivudine_COST / 5)
        # then create the gamma distribution for the cost of each drug
        self.annualZidovudineCostRVG = rvgs.Gamma(a=fit_output_zido["a"], loc=0, scale=fit_output_zido["scale"])
        self.annualLamivudineCostRVG = rvgs.Gamma(a=fit_output_lami["a"], loc=0, scale=fit_output_lami["scale"])

        # create beta distributions for annual state utility
        for utility in data.ANNUAL_STATE_UTILITY:
            # if utility is zero, add a constant 0, otherwise add a beta distribution
            if utility == 0:
                self.annualStateUtilityRVGs.append(rvgs.Constant(value=0))
            else:
                # find alpha and beta of the assumed beta distribution
                # no data available to estimate the standard deviation, so we assumed st_dev=cost / 4
                fit_output = rvgs.Beta.fit_mm(mean=utility, st_dev=utility / 4)
                # append the distribution
                self.annualStateUtilityRVGs.append(
                    rvgs.Beta(a=fit_output["a"], b=fit_output["b"]))

    def get_new_parameters(self, seed):
        """
        :param seed: seed for the random number generator used to a sample of parameter values
        :return: a new parameter set
        """

        rng = np.random.RandomState(seed=seed)

        # create a parameter set
        param = Parameters(therapy=self.therapy)

        # calculate transition probabilities
        prob_matrix = []    # probability matrix without background mortality added
        # for all health states
        for s in data.HealthStates:
            # if this state is not death
            if s != data.HealthStates.HIV_DEATH:
                # sample from the dirichlet distribution to find the transition probabilities between hiv states
                # fill in the transition probabilities out of this state
                prob_matrix.append(self.probMatrixRVG[s.value].sample(rng))

        # sampled relative risk
        rr = math.exp(self.lnRelativeRiskRVG.sample(rng))

        # calculate transition probabilities between hiv states
        if self.therapy == Therapies.MONO:
            # calculate transition probability matrix for the mono therapy
            param.probMatrix = prob_matrix

        elif self.therapy == Therapies.COMBO:
            # calculate transition probability matrix for the combination therapy
            param.probMatrix = get_prob_matrix_combo(
                prob_matrix_mono=prob_matrix,
                combo_rr=rr)

        # sample from gamma distributions that are assumed for annual state costs
        for dist in self.annualStateCostRVGs:
            param.annualStateCosts.append(dist.sample(rng))

        # sample from gamma distributions that are assumed for annual treatment costs
        zido_cost = self.annualZidovudineCostRVG.sample(rng)
        lami_cost = self.annualLamivudineCostRVG.sample(rng)

        # calculate the annual treatment cost
        if self.therapy == Therapies.MONO:
            param.annualTreatmentCost = zido_cost
        elif self.therapy == Therapies.COMBO:
            param.annualTreatmentCost = zido_cost + lami_cost

        # sample from beta distributions that are assumed for annual state utilities
        for dist in self.annualStateUtilityRVGs:
            param.annualStateUtilities.append(dist.sample(rng))

        # return the parameter set
        return param
