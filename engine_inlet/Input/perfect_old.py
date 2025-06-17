"""
Gas models for perfect gases.

This module facilitates the calculation of thermodynamic properties for any
non-reacting perfect gas (both thermally perfect and calorically perfect). A
number of gases can be specified using a convenience function SetGasType() or
the gas properties can be set explicitly using setters.
"""

import numpy as np


class GeneralGasBase:
    """
    Base class for implementing a general gas model.

    This class provides the most general implementations of the thermodynamic relations
    needed for a gas model. By overriding the functions specific to the gas model of
    interest a gas model can easily be created. This model will not be efficient.
    """

    pass


_R_UNIV = np.float64(8.31446261815324e3)

_GAS_PROP_DATA = {'air': ('Air',               # string name
                          np.float64(28.966),  # molecular weight [kg/kmol]
                          np.float64(1.4),     # CP ratio of specific heats [-]
                          ),
                  'N2':  ('Nitrogen Gas',       # string name
                          np.float64(28.0134),  # molecular weight [kg/kmol]
                          np.float64(1.4),      # CP ratio of specific heats [-]
                          ),
                  'He':  ('Helium',              # string name
                          np.float64(4.002602),  # molecular weight [kg/kmol]
                          np.float64(5/3),       # CP ratio of specific heats [-]
                          ),
                  'Ne':  ('Neon',              # string name
                          np.float64(20.179),  # molecular weight [kg/kmol]
                          np.float64(5/3),     # CP ratio of specific heats [-]
                          ),
                  'Ar':  ('Argon',             # string name
                          np.float64(39.948),  # molecular weight [kg/kmol]
                          np.float64(5/3),     # CP ratio of specific heats [-]
                          ),
                  'H2':  ('Hydrogen Gas',     # string name
                          np.float64(2.016),  # molecular weight [kg/kmol]
                          np.float64(1.405),  # CP ratio of specific heats [-]
                          ),
                  'O2':  ('Oxygen Gas',         # string name
                          np.float64(31.9988),  # molecular weight [kg/kmol]
                          np.float64(1.4),      # CP ratio of specific heats [-]
                          ),
                  'CO2': ('Carbon Dioxide',   # string name
                          np.float64(44.01),  # molecular weight [kg/kmol]
                          np.float64(1.289),  # CP ratio of specific heats [-]
                          ),
                  }


class CaloricallyPerfectGas:
    """
    Calorically perfect gas model.

    This class represents a simple calorically perfect gas model. It uses all
    of the standard CPG relations to calculate the properties of a gas. For all
    calculations two states need to be supplied, and they need to be supplied as
    keyword arguments. The typical choices are pressure, density, and temperature.

    Attributes
    ----------
        R (np.float64): Gas constant
        gamma (np.float64): Ratio of specific heats
        forumla (string): String representation of gas forumula
        name (string): Common name of gas
    """

    @staticmethod
    def get_gas_type_list():
        """
        Return list of gases that this class has configuration information.

        Returns
        -------
            List of gase names.
        """
        return [*_GAS_PROP_DATA]

    # Need another constructor that takes a string
    @classmethod
    def from_database(cls, gas):
        """
        Create new class of type 'gas'.

        The gas must be one of the gases in the classes internal database.

        Args
        ----
            gas (string): String name of gas

        Raises
        ------
            ValueError: Gas name is not in database.
        """
        g = _GAS_PROP_DATA[gas]
        rtn = cls(R=_R_UNIV/g[1], gamma=g[2],
                  formula=gas, name=g[0])

        # manually set the private custom flag to False
        rtn._CaloricallyPerfectGas__custom = False
        return rtn

    def _complete_input_state(self, p, rho, T):
        """
        Compute the missing property from the ones given.

        There must be one argument sent in as None. That property will be calculated
        using the equation of state. All three properties for the state will be returned.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return all three states. One of which was calculated using the equation of
            state.

        Throws
        ------
            Throws a value error if either zero or two variables were passed in as None.
        """
        need_p = (p is None)
        need_rho = (rho is None)
        need_T = (T is None)

        # TODO: Fix raise syntax
        # TODO: Test this exception
        # check to make sure only have 1
        if (need_p + need_rho + need_T) != 1:
            raise(ValueError,
                  "Need to pass exactly two thermodynamic properties to this function.")

        # figure out which property needs to be calculated
        if need_p:
            p = self.p(rho=rho, T=T)
        elif need_rho:
            rho = self.rho(p=p, T=T)
        else:
            T = self.T(p=p, rho=rho)

        return p, rho, T

    def __init__(self, R, gamma, formula='Custom', name=None):
        """
        Initialize the state of gas.

        Args
        ----
            R (np.float64): Gas constant
            gamma (np.float64): Ratio of specific heats
            formula (string): String representing chemical formula. The default is
                              'Custom'
            name (string): Common name for gas. The default is None.

        Throws
        ------
            ValueError when an invalid input value is given
        """
        # Check for special case of hex string being sent in
        if isinstance(R, str):
            R = np.float64.fromhex(R)
        if isinstance(gamma, str):
            gamma = np.float64.fromhex(gamma)
        self.R = np.float64(R)
        self.gamma = np.float64(gamma)
        self.formula = formula
        self.name = name
        self.__custom = True

    def __repr__(self):
        """
        Return a string representation of how this instance can be constructed.

        Returns
        -------
            String representation of instance construction.

        """
        strout = "%s" % (self.__class__.__name__)
        if (self.__custom):
            strout += "(R=%r, gamma=%r, " % (self.R, self.gamma)
            strout += "formula=%r, name=%r)" % (self.formula, self.name)
        else:
            strout += ".from_database(%r)" % (self.formula)
        return strout

    def __str__(self):
        """
        Return a readable presentation of instance.

        Returns
        -------
            Readable string representation of instance.
        """
        strout = "%s:" % (self.__class__.__name__)
        if (not self.__custom):
            strout += " (from database)\n"
        else:
            strout += "\n"
        strout += "    Formula: %s\n" % (self.formula)
        strout += "    Name: %s\n" % (self.name)
        strout += "    R: %r\n" % (self.R)
        strout += "    gamma: %r\n" % (self.gamma)

        return strout

    def __eq__(self, other):
        """
        Compare the values of this class with other class to see if they are the same.

        Args
        ----
            other (self.__class__): class to compare against.

        Returns
        -------
            Return true if all values are same between two classes.
        """
        if self.__class__ == other.__class__:
            return self.__dict__ == other.__dict__
        else:
            # TODO: Fix raise syntax
            # TODO: Test this exception
            raise(TypeError, "Invalid type in comparison")

    # Create alias-like property so user can access molecular
    #   weight, but class does not have to store that and R.
    def _set_MW(self, MWin):
        self.R = _R_UNIV/MWin

    def _get_MW(self):
        """Molecular weight of gas in kg/kmol."""
        return _R_UNIV/self.R

    MW = property(_get_MW, _set_MW)

    @property
    def R(self):
        """
        Gas constant for this gas.

        For this class the gas constant must be in units of J/kg/K.
        """
        return self.__R

    @R.setter
    def R(self, R):
        """
        Set the gas constant for this gas.

        Args
        ----
            R (np.float64): Value of gas constant

        Raises
        ------
            ValueError: If value for R passed in is less than or equal to zero
        """
        if (R <= 0):
            # TODO: Fix raise syntax
            # TODO: Test this exception
            raise(ValueError, "R must be positive")
        self.__custom = True
        self.__R = R

    @property
    def gamma(self):
        """Ratio of specific heats for this gas."""
        return self.__gamma

    @gamma.setter
    def gamma(self, gamma):
        """
        Set the ratio of specific heats for this gas.

        Args
        ----
            R (np.float64): Value of gas constant

        Raises
        ------
            ValueError: If value for gamma passed in is less than or equal to one
        """
        if (gamma <= 1):
            # TODO: Fix raise syntax
            # TODO: Test this exception
            raise(ValueError, "gamma must be greater than one")
        self.__custom = True
        self.__gamma = gamma

    @property
    def formula(self):
        """Checmical formula of the gas."""
        return self.__formula

    @formula.setter
    def formula(self, fin):
        self.__custom = True
        self.__formula = fin

    @property
    def name(self):
        """Human readable name to be used for this gas model."""
        return self.__name

    @name.setter
    def name(self, nin):
        self.__custom = True
        if nin is None:
            self.__name = self.__formula
        else:
            self.__name = nin

    def p(self, *, rho, T):
        """
        Calculate the pressure given the density and temperature.

        Return the pressure that corresponds to the given density and temperature. The
        arguments must be passed to this method in key-value form.

        Args
        ----
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Pressure.
        """
        return rho*self.R*T

    def rho(self, *, p, T):
        """
        Calculate the density given the pressure and temperature.

        Return the density that corresponds to the given pressure and temperature. The
        arguments must be passed to this method in key-value form.

        Args
        ----
            p (np.float64): Pressure.
            T (np.float64): Temperature.

        Returns
        -------
            Density.
        """
        return p/(self.R*T)

    def T(self, *, p, rho):
        """
        Calculate the temperature given the pressure and density.

        Return the temperature that corresponds to the given pressure and density. The
        arguments must be passed to this method in key-value form.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.

        Returns
        -------
            Temperature.
        """
        return p/(rho*self.R)

    def Z(self, *, p=None, rho=None, T=None):
        """
        Calculate the compressibility factor at specified state.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the compressibility factor.

        Throws
        ------
            Throws a value error if either zero or two variables were passed in as None.
        """
        p, rho, T = self._complete_input_state(p, rho, T)

        return p/(rho*self.R*T)

    def k_T(self, *, p=None, rho=None, T=None):
        """
        Calculate the isothermal compressibility.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the isothermal compressibility.

        Throws
        ------
            Throws a value error if not enough state information passed in.
        """
        if p is None:
            if (rho is None) or (T is None):
                # TODO: Fix raise syntax
                # TODO: Test this exception
                raise(ValueError, "Need to pass pressure or temperature and density"
                      + " to this function.")
            p = self.p(T=T, rho=rho)

        return 1/p

    def k_T_dp_T(self, *, p=None, rho=None, T=None):
        """
        Calculate a derivative of the isothermal compressibility.

        Calculate the partial derivative of the isothermal compressibility with respect
        to pressure holding temperature constant.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the derivative of the isothermal compressibility.

        Throws
        ------
            Throws a value error if not enough state information passed in.
        """
        if p is None:
            if (rho is None) or (T is None):
                # TODO: Fix raise syntax
                # TODO: Test this exception
                raise(ValueError, "Need to pass pressure or temperature and density"
                      + " to this function.")
            p = self.p(T=T, rho=rho)

        return -1/p**2

    def k_T_dp_rho(self, *, p=None, rho=None, T=None):
        """
        Calculate a derivative of the isothermal compressibility.

        Calculate the partial derivative of the isothermal compressibility with respect
        to pressure holding density constant.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the derivative of the isothermal compressibility.

        Throws
        ------
            Throws a value error if not enough state information passed in.
        """
        if p is None:
            if (rho is None) or (T is None):
                # TODO: Fix raise syntax
                # TODO: Test this exception
                raise(ValueError, "Need to pass pressure or temperature and density"
                      + " to this function.")
            p = self.p(T=T, rho=rho)

        return -1/p**2

    def k_T_dT_p(self, *, p=None, rho=None, T=None):
        """
        Calculate a derivative of the isothermal compressibility.

        Calculate the partial derivative of the isothermal compressibility with respect
        to temperature holding pressure constant.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the derivative of the isothermal compressibility.

        Throws
        ------
            Throws a value error if not enough state information passed in.
        """
        return np.float64(0)

    def k_T_dT_rho(self, *, p=None, rho=None, T=None):
        """
        Calculate a derivative of the isothermal compressibility.

        Calculate the partial derivative of the isothermal compressibility with respect
        to temperature holding density constant.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the derivative of the isothermal compressibility.

        Throws
        ------
            Throws a value error if not enough state information passed in.
        """
        if rho is not None:
            if (T is None) and (p is None):
                # TODO: Fix raise syntax
                # TODO: Test this exception
                raise(ValueError, "Need to pass two thermodynamic properties"
                      + " to this function.")
            if p is not None:
                if T is not None:
                    # TODO: Fix raise syntax
                    # TODO: Test this exception
                    raise(ValueError, "Can only pass two thermodynamic properties"
                          + " to this function.")
                T = self.T(p=p, rho=rho)
            else:
                p = self.p(T=T, rho=rho)

        return -1/(p*T)

    def k_T_drho_T(self, *, p=None, rho=None, T=None):
        """
        Calculate a derivative of the isothermal compressibility.

        Calculate the partial derivative of the isothermal compressibility with respect
        to density holding temperature constant.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the derivative of the isothermal compressibility.

        Throws
        ------
            Throws a value error if not enough state information passed in.
        """
        if T is not None:
            if (rho is None) and (p is None):
                # TODO: Fix raise syntax
                # TODO: Test this exception
                raise(ValueError, "Need to pass two thermodynamic properties"
                      + " to this function.")
            if p is not None:
                if rho is not None:
                    raise(ValueError, "Can only pass two thermodynamic properties"
                          + " to this function.")
                rho = self.rho(p=p, T=T)
            else:
                p = self.p(T=T, rho=rho)

        return -1/(rho*p)

    def k_T_drho_p(self, *, p=None, rho=None, T=None):
        """
        Calculate a derivative of the isothermal compressibility.

        Calculate the partial derivative of the isothermal compressibility with respect
        to density holding pressure constant.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the derivative of the isothermal compressibility.

        Throws
        ------
            Throws a value error if not enough state information passed in.
        """
        return np.float64(0)

    def alpha_v(self, *, p=None, rho=None, T=None):
        """
        Calculate the volumetric expansion coefficient.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the volumetric expansion coefficient.

        Throws
        ------
            Throws a value error if not enough state information passed in.
        """
        if T is None:
            if (rho is None) or (p is None):
                # TODO: Fix raise syntax
                # TODO: Test this exception
                raise(ValueError, "Need to pass temperature or pressure and density"
                      + " to this function.")
            T = self.T(p=p, rho=rho)

        return 1/T

    def alpha_v_dp_T(self, *, p=None, rho=None, T=None):
        """
        Calculate a derivative of the volumetric expansion coefficient.

        Calculate the partial derivative of the volumetric expansion coefficient with
        respect to pressure holding temperature constant.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the derivative of the volumetric expansion coefficient.

        Throws
        ------
            Throws a value error if not enough state information passed in.
        """
        return np.float64(0)

    def alpha_v_dp_rho(self, *, p=None, rho=None, T=None):
        """
        Calculate a derivative of the volumetric expansion coefficient.

        Calculate the partial derivative of the volumetric expansion coefficient with
        respect to pressure holding density constant.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the derivative of the volumetric expansion coefficient.

        Throws
        ------
            Throws a value error if not enough state information passed in.
        """
        if rho is not None:
            if (T is None) and (p is None):
                # TODO: Fix raise syntax
                # TODO: Test this exception
                raise(ValueError, "Need to pass two thermodynamic properties"
                      + " to this function.")
            if p is not None:
                if T is not None:
                    # TODO: Fix raise syntax
                    # TODO: Test this exception
                    raise(ValueError, "Can only pass two thermodynamic properties"
                          + " to this function.")
                T = self.T(p=p, rho=rho)
            else:
                p = self.p(T=T, rho=rho)

        return -1/(p*T)

    def alpha_v_dT_p(self, *, p=None, rho=None, T=None):
        """
        Calculate a derivative of the volumetric expansion coefficient.

        Calculate the partial derivative of the volumetric expansion coefficient with
        respect to temperature holding pressure constant.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the derivative of the volumetric expansion coefficient.

        Throws
        ------
            Throws a value error if not enough state information passed in.
        """
        if T is None:
            if (rho is None) or (p is None):
                # TODO: Fix raise syntax
                # TODO: Test this exception
                raise(ValueError, "Need to pass temperature pressure or and density"
                      + " to this function.")
            T = self.T(p=p, rho=rho)

        return -1/T**2

    def alpha_v_dT_rho(self, *, p=None, rho=None, T=None):
        """
        Calculate a derivative of the volumetric expansion coefficient.

        Calculate the partial derivative of the volumetric expansion coefficient with
        respect to temperature holding density constant.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the derivative of the volumetric expansion coefficient.

        Throws
        ------
            Throws a value error if not enough state information passed in.
        """
        if T is None:
            if (rho is None) or (p is None):
                # TODO: Fix raise syntax
                # TODO: Test this exception
                raise(ValueError, "Need to pass temperature pressure or and density"
                      + " to this function.")
            T = self.T(p=p, rho=rho)

        return -1/T**2

    def alpha_v_drho_T(self, *, p=None, rho=None, T=None):
        """
        Calculate a derivative of the volumetric expansion coefficient.

        Calculate the partial derivative of the volumetric expansion coefficient with
        respect to density holding temperature constant.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the derivative of the volumetric expansion coefficient.

        Throws
        ------
            Throws a value error if not enough state information passed in.
        """
        return np.float64(0)

    def alpha_v_drho_p(self, *, p=None, rho=None, T=None):
        """
        Calculate a derivative of the volumetric expansion coefficient.

        Calculate the partial derivative of the volumetric expansion coefficient with
        respect to density holding pressure constant.

        Args
        ----
            p (np.float64): Pressure.
            rho (np.float64): Density.
            T (np.float64): Temperature.

        Returns
        -------
            Return the derivative of the volumetric expansion coefficient.

        Throws
        ------
            Throws a value error if not enough state information passed in.
        """
        if p is not None:
            if (rho is None) and (T is None):
                # TODO: Fix raise syntax
                # TODO: Test this exception
                raise(ValueError, "Need to pass two thermodynamic properties"
                      + " to this function.")
            if T is not None:
                if rho is not None:
                    # TODO: Fix raise syntax
                    # TODO: Test this exception
                    raise(ValueError, "Can only pass two thermodynamic properties"
                          + " to this function.")
                rho = self.rho(p=p, T=T)
            else:
                T = self.T(p=p, rho=rho)

        return 1/(rho*T)
