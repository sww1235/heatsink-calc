"""Class representations of materials used in heatsinks."""


class _material():

    def __init__(self, expansionCoefficient, thermalConductivity):
        self.expansionCoefficient = expansionCoefficient  #  m/(m Kelvin)
        self.thermalConductivity = thermalConductivity  # W/(m Kelvin)


class Aluminium_6063(material):
    """6063 Aluminum Alloy with physical properties."""

    def __init__(self):
        """Initialize subclass of material with static values."""
        super(Child, self).__init__(69 * 10**(-6), 201)
