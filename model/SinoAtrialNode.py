import numpy as np


class SinoAtrialNode:
    """
    The sinoatrial node model implementation is based on the work by Fabbri et al. [1],
    which provides a comprehensive computational analysis of the human sinus node action potential.
    The CellML model used here is available on the CellML repository.
    We thank the authors for sharing their research and model.

    This class implements the model described by Fabbri et al.
    The values for constants and initial conditions in this JSON structure were obtained from the
    CellML model available at:
    (https://models.cellml.org/e/568/HumanSAN_Fabbri_Fantini_Wilders_Severi_2017.cellml/view)

    These values were derived from the original publication:
    Fabbri, Alan, et al. "Computational analysis of the human sinus node action potential: model
    development and effects of mutations." The Journal of Physiology 595.7 (2017): 2365-2396.
    (DOI: doi.org/10.1113/JP273259)

    Attributes:
        NUM_ALGEBRAIC (int): Number of algebraic variables in the model.
        NUM_STATES (int): Number of state variables in the model.
        NUM_CONSTANTS (int): Number of constant parameters in the model.
        c (np.array): Constant parameters.
        y (np.array): State variables.

    Methods:
        __init__(constants, initial_conditions, state_descriptions, constant_descriptions, acetylcholine, noradrenaline): Initializes the SAN model.
        update(acetylcholine, noradrenaline): Updates the SAN model's initial conditions with new concentrations of ACh and Nor.
        calculate_constants(): Calculates derived constants based on model parameters.
        calculate_derivatives(t, y): Calculates derivatives of state variables for numerical integration.
        print_constants(self): Prints the model's constant parameters with descriptions.
        print_initial_conditions(self): Prints the model's initial conditions with descriptions.
        get_state_labels(self): Returns a list of labels for state variables.
        get_constant_labels(self): Returns a list of labels for constant parameters.
        info(self): Prints information about the model and its parameters.
    """

    def __init__(self, constants, initial_conditions, state_descriptions, constant_descriptions, acetylcholine=0,
                 noradrenaline=0):
        """
        Initializes the SinoAtrialNode class with the given parameters.

        Args:
            constants (np.array): Array of constant parameters for the model.
            initial_conditions (np.array): Array of initial state variable values.
            state_descriptions (list of str): Descriptions of state variables.
            constant_descriptions (list of str): Descriptions of constant parameters.
            acetylcholine (float, optional): Concentration of Acetylcholine. Defaults to 0.
            noradrenaline (float, optional): Concentration of Noradrenaline. Defaults to 0.

        Returns:
            None.
        """
        self.NUM_ALGEBRAIC = 101
        self.NUM_STATES = 33
        self.NUM_CONSTANTS = 116
        self.c = np.zeros(self.NUM_CONSTANTS)
        self.c[0:91] = constants
        self.y = np.array(initial_conditions)
        self.y0 = self.y
        self.state_descriptions = state_descriptions
        self.constant_descriptions = constant_descriptions
        # Update the model with the provided concentrations of Acetylcholine and Noradrenaline
        self.update(acetylcholine, noradrenaline)

    def update(self, acetylcholine=0, noradrenaline=0):
        """
        Updates the SinoAtrialNode model's state variables based on the given concentrations of acetylcholine and noradrenaline.

        Args:
            acetylcholine (float): Concentration of Acetylcholine to be used in the model.
            noradrenaline (float): Concentration of Noradrenaline to be used in the model.

        Returns:
            None.
        """
        self.c[9] = acetylcholine
        self.c[10] = noradrenaline
        self.calculate_constants()

    def calculate_constants(self):
        """
        Calculates derived constants based on the model's parameter values.

        Returns:
            None.
        """
        self.c[91] = (self.c[0] * self.c[1]) / self.c[2]
        self.c[92] = self.c[16] / (self.c[13] / (self.c[13] + self.c[17]))
        self.c[93] = np.piecewise(self.c[10], [self.c[10] > 0.0, self.c[9] > 0.0],
                                  [-0.25, lambda x: (0.7 * x) / (9.0e-05 + x), 0.0])
        self.c[94] = np.piecewise(self.c[10], [self.c[10] > 0.0],
                                  [lambda x: 1.2 * self.c[88], self.c[88]])
        self.c[95] = np.piecewise(self.c[9], [self.c[9] > 0.0],
                                  [lambda x: -1.0 - (9.89800 * np.power(1.0 * x, 0.618)) /
                                             (np.power(1.0 * x, 0.618) + 0.00122423), 0.0])
        self.c[96] = self.c[91] * np.log(self.c[13] / self.c[12])
        self.c[97] = self.c[92] / (self.c[18] + 1.0)
        self.c[98] = self.c[61] * (1.0 - self.c[93])
        self.c[99] = np.piecewise(self.c[10], [self.c[10] > 0.0], [7.50000, 0.0])
        self.c[100] = self.c[18] * self.c[97]
        self.c[101] = (self.c[97] * self.c[13]) / (self.c[13] + self.c[17])
        self.c[102] = np.piecewise(self.c[10], [self.c[10] > 0.0], [1.20000, 1.0])
        self.c[103] = (self.c[100] * self.c[13]) / (self.c[13] + self.c[17])
        self.c[104] = self.c[11] / (self.c[33] + self.c[11])
        self.c[105] = np.piecewise(self.c[10], [self.c[10] > 0.0], [1.23000, 1.0])
        self.c[106] = (0.310000 * self.c[9]) / (self.c[9] + 9.0e-05)
        self.c[107] = np.piecewise(self.c[10], [self.c[10] > 0.0], [-8.00000, 0.0])
        self.c[108] = np.piecewise(self.c[10], [self.c[10] > 0.0], [-27.0000, 0.0])
        self.c[109] = 1.00000e-09 * np.pi * np.power(self.c[82], 2.0) * self.c[83]
        self.c[110] = 1.00000e-09 * 2.00000 * np.pi * self.c[84] * (self.c[82] - self.c[84] / 2.0) * self.c[83]
        self.c[111] = self.c[79] * self.c[109]
        self.c[112] = self.c[80] * self.c[109] - self.c[110]
        self.c[113] = self.c[81] * self.c[109]
        self.c[114] = np.piecewise(self.c[10], [self.c[10] > 0.00000], [-14.0000, 0.00000])
        with np.errstate(divide='ignore'):
            self.c[115] = (3.59880 - 0.0256410) / \
                          (1.00000 + 1.21550e-06 / np.power(1.00000 * self.c[9], 1.69510)) + 0.0256410

    def calculate_derivatives(self, t, y):
        """
        Calculates the derivatives of the SinoAtrialNode model's state variables.

        Args:
            t (float): The current time point in the simulation.
            y (numpy.array): The current state of the model, including all state variables.

        Returns:
            dydt (numpy.array): The derivatives of the state variables with respect to time.
        """
        a = np.zeros(self.NUM_ALGEBRAIC)
        dydt = np.zeros(self.NUM_STATES)
        self.y = y
        a[6] = self.c[69] * self.c[78] * (1.0 - (self.y[22] + self.y[18])) - self.c[75] * self.y[18]
        dydt[18] = a[6]
        a[1] = self.c[47] / (self.c[47] + self.y[1])
        a[7] = (0.001 * a[1]) / self.c[46]
        dydt[8] = (a[1] - self.y[8]) / a[7]
        a[2] = self.c[51] - (self.c[51] - self.c[52]) / (1.0 + np.power(self.c[53] / self.y[15], self.c[54]))
        a[8] = self.c[55] / a[2]
        a[17] = self.c[56] * a[2]
        dydt[11] = (self.c[57] * self.y[14] - a[17] * self.y[1] * self.y[11]) - \
                   (a[8] * (np.power(self.y[1], 2.0)) * self.y[11] - self.c[58] * self.y[12])
        dydt[12] = (a[8] * (np.power(self.y[1], 2.0)) * self.y[11] - self.c[58] * self.y[12]) - \
                   (a[17] * self.y[1] * self.y[12] - self.c[57] * self.y[13])
        dydt[13] = (a[17] * self.y[1] * self.y[12] - self.c[57] * self.y[13]) - \
                   (self.c[58] * self.y[13] - a[8] * (np.power(self.y[1], 2.0)) * self.y[14])
        dydt[14] = (self.c[58] * self.y[13] - a[8] * (np.power(self.y[1], 2.0)) * self.y[14]) - \
                   (self.c[57] * self.y[14] - a[17] * self.y[1] * self.y[11])
        a[5] = np.piecewise(t, [np.greater(t, self.c[5]) & np.less(t, self.c[5] + self.c[6]),
                                ~(np.greater(t, self.c[5]) & np.less(t, self.c[5] + self.c[6]))],
                            [self.c[7], self.c[8]])
        a[9] = np.piecewise(self.c[4], [np.greater_equal(self.c[4], 1.0),
                                        ~np.greater_equal(self.c[4], 1.0)], [a[5], self.y[0]])
        a[10] = 1.0 / ((0.36 * (((a[9] + 148.8) - self.c[95]) - self.c[99])) /
                       (np.exp(0.066 * (((a[9] + 148.8) - self.c[95]) - self.c[99])) - 1.0) +
                       (0.1 * (((a[9] + 87.3) - self.c[95]) - self.c[99])) /
                       (1.0 - np.exp(-0.2 * (((a[9] + 87.3) - self.c[95]) - self.c[99])))) - 0.054
        a[29] = np.piecewise(a[9], [np.less(a[9], -(((80.0 - self.c[95]) - self.c[99]) - self.c[20])),
                                    ~np.less(a[9], -(((80.0 - self.c[95]) - self.c[99]) - self.c[20]))],
                             [0.01329 + 0.999210 / (1.0 + np.exp(
                                 ((((a[9] + 97.1340) - self.c[95]) - self.c[99]) - self.c[20]) / 8.17520)),
                              0.0002501 * np.exp(-(((a[9] - self.c[95]) - self.c[99]) - self.c[20]) / 12.8610)])
        dydt[3] = (a[29] - self.y[3]) / a[10]
        a[14] = 1.0 / (1.0 + np.exp((a[9] + 37.4 + self.c[44]) / (5.30 + self.c[45])))
        a[33] = 0.001 * (44.3 + 230.0 * np.exp(-(np.power((a[9] + 36.0) / 10.0, 2.0))))
        dydt[7] = (a[14] - self.y[7]) / a[33]
        a[15] = 1.0 / (1.0 + np.exp(-(a[9] + 38.3) / 5.5))
        a[34] = 0.001 / (1.068 * np.exp((a[9] + 38.3) / 30.0) + 1.068 * np.exp(-(a[9] + 38.3) / 30.0))
        dydt[9] = (a[15] - self.y[9]) / a[34]
        a[16] = 1.0 / (1.0 + np.exp((a[9] + 58.7) / 3.8))
        a[35] = 1.0 / (16.6700 * np.exp(-(a[9] + 75.0) / 83.3) + 16.67 * np.exp((a[9] + 75.0) / 15.38)) + self.c[49]
        dydt[10] = (a[16] - self.y[10]) / a[35]
        a[37] = 0.009 / (1.0 + np.exp((a[9] + 5.0) / 12.0)) + 0.0005
        a[19] = 1.0 / (1.0 + np.exp((a[9] + 6.0) / -8.6))
        dydt[24] = (a[19] - self.y[24]) / a[37]
        a[38] = 0.590 / (1.0 + np.exp((a[9] + 60.0) / 10.0)) + 3.05
        a[20] = 1.0 / (1.0 + np.exp((a[9] + 7.5) / 10.0))
        dydt[25] = (a[20] - self.y[25]) / a[38]
        a[21] = 1.0 / (1.0 + np.exp((a[9] + 49.0) / 13.0))
        a[39] = 0.001 * 0.6 * (65.17 / (
                0.57 * np.exp(-0.08 * (a[9] + 44.0)) + 0.065 * np.exp(0.1 * (a[9] + 45.93))) + 10.1)
        dydt[26] = (a[21] - self.y[26]) / a[39]
        a[22] = 1.0 / (1.0 + np.exp(-(a[9] - 19.3) / 15.0))
        a[40] = 0.001 * 0.66 * 1.4 * \
                (15.59 / (1.03700 * np.exp(0.09 * (a[9] + 30.61)) + 0.369 * np.exp(-0.12 * (a[9] + 23.84))) + 2.98)
        dydt[27] = (a[22] - self.y[27]) / a[40]
        a[23] = 1.0 / (1.0 + np.exp(-(a[9] + 10.0144) / 7.6607))
        a[41] = 0.846554 / (4.2 * np.exp(a[9] / 17.0) + 0.15 * np.exp(-a[9] / 21.6))
        dydt[28] = (a[23] - self.y[28]) / a[41]
        a[42] = 1.0 / (30.0 * np.exp(a[9] / 10.0) + np.exp(-a[9] / 12.0))
        dydt[29] = (a[23] - self.y[29]) / a[42]
        a[43] = 1.0 / (1.0 + np.exp((a[9] + 28.6) / 17.1))
        a[26] = 1.0 / (100.0 * np.exp(-a[9] / 54.6450) + 656.0 * np.exp(a[9] / 106.157))
        dydt[30] = (a[43] - self.y[30]) / a[26]
        a[28] = 10.0 * np.exp(0.0133 * (a[9] + 40.0))
        a[45] = self.c[115] / (self.c[115] + a[28])
        a[51] = 1.0 / (self.c[115] + a[28])
        dydt[32] = (a[45] - self.y[32]) / a[51]
        a[12] = 1.0 / (1.0 + np.exp((a[9] + 69.804) / 4.4565))
        a[31] = 20.0 * np.exp(-0.125 * (a[9] + 75.0))
        a[47] = 2000.0 / (320.0 * np.exp(-0.1 * (a[9] + 75.0)) + 1.0)
        a[53] = 1.0 / (a[31] + a[47])
        dydt[5] = (a[12] - self.y[5]) / a[53]
        a[27] = np.power(1.0 / (1.0 + np.exp(-((a[9] + 0.6383) - self.c[114]) / 10.7071)), 1.0 / 2)
        a[44] = 28.0 / (1.0 + np.exp(-((a[9] - 40.0) - self.c[114]) / 3.0))
        a[50] = 1.0 * np.exp(-((a[9] - self.c[114]) - 5.0) / 25.0)
        a[56] = 1.0 / (a[44] + a[50])
        dydt[31] = (a[27] - self.y[31]) / a[56]
        a[11] = 1.0 / (1.0 + np.exp(-(a[9] + 42.0504) / 8.3106))
        a[30] = a[9] + 41.0
        a[46] = np.piecewise(a[30], [np.less(np.fabs(a[30]), self.c[40]), ~np.less(np.fabs(a[30]), self.c[40])],
                             [2000.0, (200.0 * a[30]) / (1.0 - np.exp(-0.1 * a[30]))])
        a[52] = 8000.0 * np.exp(-0.056 * (a[9] + 66.0))
        a[57] = 1.0 / (a[46] + a[52])
        dydt[4] = (a[11] - self.y[4]) / a[57]
        a[13] = 1.0 / (1.0 + np.exp(
            -((a[9] - self.c[43]) - self.c[107]) / (self.c[42] * (1.0 + self.c[108] / 100.0))))
        a[32] = np.piecewise(a[9], [np.equal(a[9], -41.80), np.equal(a[9], 0.0), np.equal(a[9], -6.8),
                                    ~(np.equal(a[9], -41.8) | np.equal(a[9], 0.0) | np.equal(a[9], -6.8))],
                             [-41.8, 0.0, -6.80001, a[9]])
        a[48] = (-0.02839 * (a[32] + 41.8)) / (np.exp(-(a[32] + 41.8) / 2.5) - 1.0) - (
                0.0849 * (a[32] + 6.8)) / (np.exp(-(a[32] + 6.8) / 4.8) - 1.0)
        a[54] = np.piecewise(a[9], [np.equal(a[9], -1.8), ~np.equal(a[9], -1.8)], [-1.80001, a[9]])
        a[58] = (0.01143 * (a[54] + 1.8)) / (np.exp((a[54] + 1.8) / 2.5) - 1.0)
        a[60] = 0.001 / (a[48] + a[58])
        dydt[6] = (a[13] - self.y[6]) / a[60]
        a[18] = self.y[2]
        a[36] = self.c[91] * np.log(self.c[11] / a[18])
        a[61] = self.c[102] * self.c[23] * (np.power(1.0 + np.power(self.c[21] / self.c[13], 1.2), -1.0)) * \
                (np.power(1.0 + np.power(self.c[22] / a[18], 1.3), -1.0)) * \
                (np.power(1.0 + np.exp(-((a[9] - a[36]) + 110.0) / 20.0), -1.0))
        a[63] = np.exp((-self.c[26] * a[9]) / (2.0 * self.c[91]))
        a[69] = 1.0 + (self.c[14] / self.c[36]) * (1.0 + np.exp((self.c[27] * a[9]) / self.c[91])) + (
                    self.c[11] / self.c[34]) * \
                (1.0 + (self.c[11] / self.c[35]) * (1.0 + self.c[11] / self.c[33]))
        a[71] = ((((self.c[11] / self.c[34]) * self.c[11]) / self.c[35]) * (1.0 + self.c[11] / self.c[33]) *
                 np.exp((-self.c[26] * a[9]) / (2.0 * self.c[91]))) / a[69]
        a[70] = ((self.c[14] / self.c[36]) * np.exp((self.c[27] * a[9]) / self.c[91])) / a[69]
        a[67] = np.exp((self.c[26] * a[9]) / (2.0 * self.c[91]))
        a[62] = a[18] / (self.c[28] + a[18])
        a[72] = a[63] * self.c[104] * (a[71] + a[70]) + a[70] * a[67] * (a[62] + a[63])
        a[64] = 1.0 + (self.y[1] / self.c[29]) * (
                    1.0 + np.exp((-self.c[25] * a[9]) / self.c[91]) + a[18] / self.c[32]) + \
                (a[18] / self.c[30]) * (1.0 + (a[18] / self.c[31]) * (1.0 + a[18] / self.c[28]))
        a[65] = ((self.y[1] / self.c[29]) * np.exp((-self.c[25] * a[9]) / self.c[91])) / a[64]
        a[66] = ((((a[18] / self.c[30]) * a[18]) / self.c[31]) * (1.0 + a[18] / self.c[28]) *
                 np.exp((self.c[26] * a[9]) / (2.0 * self.c[91]))) / a[64]
        a[68] = a[67] * a[62] * (a[66] + a[65]) + a[63] * a[65] * (self.c[104] + a[67])
        a[73] = a[66] * a[62] * (a[71] + a[70]) + a[65] * a[71] * (a[62] + a[63])
        a[74] = a[71] * self.c[104] * (a[66] + a[65]) + a[66] * a[70] * (self.c[104] + a[67])
        a[75] = ((1.0 - self.c[37]) * self.c[24] * (a[68] * a[70] - a[72] * a[65])) / (a[72] + a[68] + a[73] + a[74])
        a[76] = self.c[91] * np.log((self.c[11] + 0.12 * self.c[13]) / (a[18] + 0.12 * self.c[12]))
        a[77] = self.c[38] * (np.power(self.y[4], 3.0)) * self.y[5] * (a[9] - a[76])
        a[78] = self.c[39] * (np.power(self.y[4], 3.0)) * (a[9] - a[76])
        a[79] = a[77] + a[78]
        a[49] = self.y[3] * self.c[103] * (a[9] - a[36]) * (1.0 - self.c[19])
        a[82] = ((1.85e-05 * self.c[41] * (a[9] - 0.0)) / (
                    self.c[91] * (1.0 - np.exp((-1.0 * (a[9] - 0.0)) / self.c[91])))) * \
                (a[18] - self.c[11] * np.exp((-1.0 * (a[9] - 0.0)) / self.c[91])) * self.y[6] * self.y[7] * self.y[8]
        dydt[2] = ((1.0 - self.c[15]) * -1.0 * (a[79] + a[49] + a[82] + 3.0 * a[61] + 3.0 * a[75])) / \
                  (1.0 * (self.c[112] + self.c[110]) * self.c[2])
        a[90] = self.c[71] * self.y[1] * (1.0 - self.y[20]) - self.c[76] * self.y[20]
        dydt[20] = a[90]
        a[84] = ((2.0 * self.c[48] * a[9]) / (self.c[91] * (1.0 - np.exp((-1.0 * a[9] * 2.0) / self.c[91])))) * \
                (self.y[1] - self.c[14] * np.exp((-2.0 * a[9]) / self.c[91])) * self.y[9] * self.y[10]
        a[80] = ((2.0 * self.c[41] * (a[9] - 0.0)) / (
                    self.c[91] * (1.0 - np.exp((-1.0 * (a[9] - 0.0) * 2.0) / self.c[91])))) * \
                (self.y[1] - self.c[14] * np.exp((-2.0 * (a[9] - 0.0)) / self.c[91])) * self.y[6] * self.y[7] * self.y[
                    8]
        a[86] = self.c[50] * self.y[12] * (self.y[15] - self.y[1])
        a[88] = (self.y[1] - self.y[17]) / self.c[59]
        dydt[1] = (a[86] * self.c[111]) / self.c[110] - (((a[80] + a[84]) - 2.0 * a[75]) /
                                                         (2.0 * self.c[2] * self.c[110]) + a[88] + self.c[66] * a[90])
        a[93] = self.c[68] * self.y[17] * (1.0 - self.y[21]) - self.c[73] * self.y[21]
        dydt[21] = a[93]
        a[91] = self.c[98] / (1.0 + np.exp((-self.y[17] + self.c[62]) / self.c[63]))
        a[94] = (self.y[16] - self.y[15]) / self.c[60]
        dydt[16] = a[91] - (a[94] * self.c[111]) / self.c[113]
        a[96] = self.c[70] * self.y[17] * (1.0 - (self.y[22] + self.y[18])) - self.c[74] * self.y[22]
        dydt[22] = a[96]
        a[97] = self.c[72] * self.y[15] * (1.0 - self.y[23]) - self.c[77] * self.y[23]
        dydt[23] = a[97]
        dydt[15] = a[94] - (a[86] + self.c[67] * a[97])
        a[99] = self.c[71] * self.y[17] * (1.0 - self.y[19]) - self.c[76] * self.y[19]
        dydt[19] = a[99]
        dydt[17] = (1.0 * (a[88] * self.c[110] - a[91] * self.c[113])) / \
                   self.c[112] - (self.c[66] * a[99] + self.c[64] * a[93] + self.c[65] * a[96])
        a[55] = self.y[3] * self.c[101] * (a[9] - self.c[96]) * (1.0 - self.c[19])
        a[59] = a[49] + a[55]
        a[89] = self.c[87] * (a[9] - self.c[96]) * (0.9 * self.y[29] + 0.1 * self.y[28]) * self.y[30]
        a[92] = self.c[91] * np.log((self.c[13] + 0.12 * self.c[11]) / (self.c[12] + 0.12 * a[18]))
        a[95] = self.c[94] * (a[9] - a[92]) * (np.power(self.y[31], 2.0))
        a[87] = self.c[86] * (a[9] - self.c[96]) * self.y[26] * self.y[27]
        a[81] = ((0.000365 * self.c[41] * (a[9] - 0.0)) /
                 (self.c[91] * (1.0 - np.exp((-1.0 * (a[9] - 0.0)) / self.c[91])))) * \
                (self.c[12] - self.c[13] * np.exp((-1.0 * (a[9] - 0.0)) / self.c[91])) * self.y[6] * self.y[7] * self.y[
                    8]
        a[83] = (a[80] + a[81] + a[82]) * (1.0 - self.c[106]) * 1.0 * self.c[105]
        a[98] = np.piecewise(self.c[9], [np.greater(self.c[9], 0.0), ~np.greater(self.c[9], 0.0)],
                             [self.c[90] * self.c[89] * (a[9] - self.c[96]) *
                              (1.0 + np.exp((a[9] + 20.0) / 20.0)) * self.y[32], 0.0])
        a[85] = self.c[85] * self.y[24] * self.y[25] * (a[9] - self.c[96])
        a[100] = a[59] + a[89] + a[95] + a[87] + a[61] + a[75] + a[79] + a[83] + a[84] + a[98] + a[85]
        dydt[0] = -a[100] / self.c[3]
        return dydt

    def print_constants(self):
        """
        Print the constant parameters of the SinoAtrialNode model.

        This method prints a list of constant parameter values that define the behavior of the SinoAtrialNode model.
        The constants include various physiological and biophysical properties used in the model.

        Args:
            None

        Returns:
            None
        """
        print("Constants:")
        for value, description in zip(self.c, self.constant_descriptions):
            parts = description.split(" in ")
            variable_name = parts[0]
            component_name = parts[1].split("[")[0].strip()
            unit = parts[1].split("[")[1].split("]")[0]
            print(f"{variable_name} in {component_name}: {value} {unit}")

    def print_initial_conditions(self):
        """
        Print the initial conditions of the SinoAtrialNode model.

        This method prints the initial values of the state variables that represent the internal states of the SinoAtrialNode model.
        The state variables describe the physiological states of the model at the beginning of the simulation.

        Args:
            None

        Returns:
            None
        """
        print("Initial Conditions:")
        for value, description in zip(self.y0, self.state_descriptions):
            label = description.split("[")[0].strip()
            unit = description.split("[")[1][:-1]  # Extract the unit within square brackets
            print(f"{label}: {value} {unit}")

    def get_state_labels(self):
        """
        Get the labels of the state variables.

        This method returns a list of labels for the state variables. These labels describe the different internal physiological states
        that the SinoAtrialNode model simulates.

        Args:
            None

        Returns:
            List[str]: A list of strings containing the labels of the state variables.
        """
        return self.state_descriptions

    def get_constant_labels(self):
        """
        Get the labels of the constant parameters.

        This method returns a list of labels for the constant parameters used in the SinoAtrialNode model.
        These labels provide information about the different physiological and biophysical properties that affect the model's behavior.

        Args:
            None

        Returns:
            List[str]: A list of strings containing the labels of the constant parameters.
        """
        return self.constant_descriptions

    def info(self):
        """
       Print information about the SinoAtrialNode model.

       This method provides general information about the SinoAtrialNode model, including its source, references,
       and the source of constant values. It also calls the print_constants and print_initial_conditions methods to provide
       detailed information about the model's constants and initial conditions.

       Args:
           None

       Returns:
           None
       """
        print("Model Information:\n")
        print("This class implements the model described by Fabbri et al. [1]")
        print("The values for constants and initial conditions in this JSON structure were obtained from the")
        print(
            "CellML model available at: https://models.cellml.org/e/568/HumanSAN_Fabbri_Fantini_Wilders_Severi_2017.cellml")
        print("These values were derived from the original of Fabbri et al. [1].\n")
        print("[1] Fabbri, Alan, et al. \"Computational analysis of the human sinus node action potential:")
        print("model development and effects of mutations. The Journal of Physiology 595.7 (2017): 2365-2396.")
        print("DOI: https://doi.org/10.1113/JP273259\n\n")
        self.print_constants()
        self.print_initial_conditions()
