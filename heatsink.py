"""Python Script to Calculate parameters of various heatsinks."""

import math


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
                 finThickness, cylinderDiameter, heatsinkLength):
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
        """
        self.finSpacing = finSpacing
        self.finRadius = finRadius
        self.finThickness = finThickness
        self.cylinderDiameter = cylinderDiameter
        self.heatsinkLength = heatsinkLength

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
        if 0.1 <= L/D <= 8:
            self.nn0 = ((3.36 + (0.087 * (L/D)))
                        * math.sqrt(acc)
                        * (finSpacing / ahs)
                        )

        if 0.1 <= tNf/D <= 8:
            self.nnOut = ((0.499 - (0.026 * math.log(tNf/D)))
                          * math.pow(Rab, -0.25)
                          * (Aout/Ahs)
                          )
        if (0.1 <= d/D <= 8) and (2.9 * 10**4 <= Rad <= 2.3 * 10**5):
            nnInT = (0.573-(0.184 * (d/D)) +( 0.0388 * (d/D)**2))
        n = 1
        self.nnIn = (math.pow(math.pow(nnInT, -n)
                     + math.pow(nnInFD, -n), (-1/n)
                              )
                     * (Ain/ahs)
                     )
        self.nn = (self.nnIn + self.nnOut + self.nn0)
        super(Child, self).__init__(material, self.__name__)
