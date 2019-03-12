"""Class representations of heatsinks."""

import math
from scipy import constants as const

from materials import Aluminium_6063 as aluminium


class Heatsink:
    """
    A Heatsink.

    Extended by form factor subclasses
    """

    def __init__(self, material, configuration):
        """Init material and configuration variables."""
        self.material = material
        self.configuration = configuration


class CylindricalAnnularFin(Heatsink):
    """Extend base heatsink class with a cylindrical annular fin heatsink."""

    def __init__(self, material, finSpacing, finRadius,
                 finThickness, cylinderDiameter, numberOfFins,
                 ambAirTemp, maxJunctionTemp, maxSurfaceTemp):
        """
        Init remainder of class variables.

        NOTE: all models are based off of the finSpacing variable
        NOTE: using the simplified model for calculation efficiency.

        finSpacing      : gap between adjacent fins
        finRadius       : radius of fin minus central support cylinder
                          (alternatively, fin depth)
        finThickness    : thickness of individual fin
        cylinderDiameter: diameter of support cylinder
        heatsinkLength  : overall axial length of heatsink
        overall diameter: outside diameter of heatsink including fins.
        """
        self.finSpacing = finSpacing  # in meters
        self.finRadius = finRadius  # in meters
        self.finThickness = finThickness  # in meters
        self.cylinderDiameter = cylinderDiameter  # in meters
        self.numberOfFins = numberofFins
        self.heatsinkLength = ((self.finThickness * self.numberOfFins)
                               + ((self.numberOfFins - 1) * self.finSpacing))
        self.overallDiameter = self.cylinderDiameter + (2 * finRadius)
        self.ambAirTemp = ambAirTemp  # degrees kelvin
        self.maxJunctionTemp = maxJunctionTemp
        self.maxSurfaceTemp = maxSurfaceTemp

        """
        NOTE: in order to prevent ridiculously long variable names, all
        Nusselt Numbers are abbreviated as follows:
        nn      = Nusselt Number
        nn0     = Nusselt Number 0 (Diffusive Limit)
        nnOut   = Nusselt Number for outer surfaces
        nnIn    = Nusselt Number for inner surfaces
        nnInT   = Nusselt Number for the thin boundry layer of inner surface
        nnInFD  = Nusselt Number for fully developed regime inner surface
        """
        # thermal diffusivity of air at atmospheric pressure at 25C
        alpha = 22.39 * 10**(-6)  # (meters^2) / seconds
        # Volumetric coefficient of thermal expansion
        beta = aluminium.expansionCoefficient  # 1/kelvin
        heatsinkSurfaceTemp =  # TODO kelvin
        # at atmospheric pressure at 25C
        kinematicViscosity = 15.52 * 10**(-6)  # meter^2/second
        deltaT = heatsinkSurfaceTemp - ambAirTemp  # kelvin
        hLoD = self.heatsinkLength / self.overallDiameter
        cDoD = self.cylinderDiameter / self.overallDiameter
        oneChannelArea = (math.pi * (((self.overallDiameter**2
                                      - self.cylinderDiameter**2) / 2)
                                     + (self.cylinderDiameter
                                        * self.finSpacing)))
        # area of circumscribed cylinder
        areaCC = (math.pi * (((self.overallDiameter**2) / 2)
                  + self.overallDiameter * self.heatsinkLength))  # meter^2
        # inner surface area of heatsink
        areaIn = (self.numberOfFins - 1) * oneChannelArea  # meter^2
        # outer surface area of heatsink
        areaOut = (math.pi * (((self.overallDiameter**2) / 2)
                   + (self.numberOfFins
                      * self.overallDiameter
                      * self.finThickness)))  # meter^2
        # overall area of heatsink
        areaHS = areaIn + areaOut  # meter^2
        RayleighNbrFinSpacing = ((const.g
                                  * beta
                                  * deltaT
                                  * self.finSpacing**4)
                                 / (kinematicViscosity
                                    * alpha
                                    * self.overallDiameter))
        RayleighNbrOverallDiameter = ((const.g
                                       * beta
                                       * deltaT
                                       * self.overallDiameter**3)
                                      / (kinematicViscosity * alpha))
        if 0.1 <= hLoD <= 8:
            self.nn0 = ((3.36 + (0.087 * hLoD))
                        * math.sqrt(areaCC)
                        * (self.finSpacing / areaHS)
                        )

        if 0.1 <= (self.finThickness
                   * self.numberOfFins
                   / self.overallDiameter) <= 8:
            self.nnOut = ((0.499 - (0.026 * math.log(self.finThickness
                                                     * self.numberOfFins
                                                     / self.overallDiameter)))
                          * math.pow(RayleighNbrFinSpacing, 0.25)
                          * (areaOut/areaHS)
                          )
        if (0.1 <= cdoD <= 8) and (2.9 * 10**4
                                   <= RayleighNbrOverallDiameter
                                   <= 2.3 * 10**5):
            nnInT = ((0.573-(0.184 * cdoD) + (0.0388 * cdoD**2))
                     * math.pow(RayleighNbrFinSpacing, 0.25))
            nnInFD = (((0.0323
                       - (0.0517 * cdoD)
                       + (0.11 * cdoD**2))
                       * math.pow(RayleighNbrFinSpacing, 0.25))
                      + (0.0516 + (0.0154 * cdoD)
                      - (0.0433 * cdoD**2)
                      + (0.0792 * cdoD**3)) * RayleighNbrFinSpacing)
        n = 1
        self.nnIn = (math.pow(math.pow(nnInT, -n)
                     + math.pow(nnInFD, -n), (-1/n)
                              )
                     * (areaIn/areaHS)
                     )
        self.nn = (self.nnIn + self.nnOut + self.nn0)

        super(Child, self).__init__(material, self.__name__)

        """
        Nusselt number = (Qconv * b) / (Ahs deltaT k)
        Qconv = heat flow rate by convection (Watts)
        b = finSpacing (meters)
        Ahs = Area of heatsink (meter^2)
        deltaT = temperature difference between surface temp of
                heatsink and ambient air temp.
        k = thermal conductivity of material (Watts / (meter kelvin))
        """
