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


        # treatment relative risk
        rr_ci = [0.365, 0.71]  # confidence interval of the treatment relative risk

        # find the mean and st_dev of the normal distribution assumed for ln(RR)
        # sample mean ln(RR)

        # sample standard deviation of ln(RR)

        # create a normal distribution for ln(RR)

        # create gamma distributions for annual state cost
        for cost in data.ANNUAL_STATE_COST:

            # if cost is zero, add a constant 0, otherwise add a gamma distribution
            if cost == 0:

            else:
                # find shape and scale of the assumed gamma distribution
                # no data available to estimate the standard deviation, so we assumed st_dev=cost / 5

                # append the distribution

        # create a gamma distribution for annual treatment cost with each drug
        # first fit the gamma distribution to the cost of each drug

        # then create the gamma distribution for the cost of each drug

        # create beta distributions for annual state utility
        for utility in data.ANNUAL_STATE_UTILITY:
            # if utility is zero, add a constant 0, otherwise add a beta distribution
            if utility == 0:

            else:
                # find alpha and beta of the assumed beta distribution
                # no data available to estimate the standard deviation, so we assumed st_dev=cost / 4

                # append the distribution

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

                # sample from the dirichlet distribution to find the transition probabilities between hiv states
                # fill in the transition probabilities out of this state

        # sampled relative risk

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


        # sample from gamma distributions that are assumed for annual treatment costs


        # calculate the annual treatment cost
        if self.therapy == Therapies.MONO:
            param.annualTreatmentCost = zido_cost
        elif self.therapy == Therapies.COMBO:
            param.annualTreatmentCost = zido_cost + lami_cost

        # sample from beta distributions that are assumed for annual state utilities


        # return the parameter set
        return param
