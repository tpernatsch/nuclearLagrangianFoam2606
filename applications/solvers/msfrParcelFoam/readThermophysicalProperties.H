singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);

// Specific heat capacity [J/(kg K)]
dimensionedScalar CpRef("CpRef", dimEnergy/dimMass/dimTemperature, laminarTransport);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

// Laminar Schmidt number
dimensionedScalar Sc("Sc", dimless, laminarTransport);

// Turbulent Schmidt number
dimensionedScalar Sct("Sct", dimless, laminarTransport);