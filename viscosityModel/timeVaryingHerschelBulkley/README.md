# Time-Varying Herschel-Bulkley Viscosity Model for OpenFOAM

[![OpenFOAM v2406](https://img.shields.io/badge/OpenFOAM-v2406-blue)](https://www.openfoam.com/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


## Overview

This repository contains an implementation of a time-varying Herschel-Bulkley viscosity model for OpenFOAM v2406. The model extends the classical Herschel-Bulkley model by introducing time-dependent consistency coefficient, enabling simulation of complex rheological behaviors in time-evolving non-Newtonian fluids.

## Mathematical Formulation

The classical Herschel-Bulkley model is given by:

$$\tau = \tau_0 + k\dot{\gamma}^n \quad \text{for} \quad \tau > \tau_0$$

where:
- $\tau$ is the shear stress
- $\tau_0$ is the yield stress
- $k$ is the consistency coefficient
- $\dot{\gamma}$ is the shear rate
- $n$ is the flow behavior index

This implementation introduces time-dependency in the consistency coefficient $k(t)$ through two models:

1. **Power-law variation**: $k(t) = A \cdot t^B$
2. **Exponential variation**: $k(t) = A \cdot e^{B \cdot t}$

## Features

- ✅ Fully compatible with OpenFOAM v2406
- ✅ Two time-variation models (power-law and exponential)
- ✅ Robust numerical implementation with bounded viscosity
- ✅ Easy integration with existing OpenFOAM solvers
- ✅ Comprehensive parameter control through dictionary files

## Installation

### Prerequisites
- OpenFOAM v2406
- GNU C++ Compiler
- Make build system

### Build Instructions

1. Clone this repository:
```bash
git clone https://github.com/searchinvaild/timeVaryingHerschelBulkley.git
cd timeVaryingHerschelBulkley
```

2. Set OpenFOAM environment:
```bash
source /opt/openfoam-v2406/etc/bashrc
```

3. Compile the library:
```bash
wmake libso
```

4. The compiled library will be available in `$FOAM_USER_LIBBIN`

## Usage

### Basic Configuration

In your case's `constant/transportProperties` file:

```cpp
transportModel  timeVaryingHerschelBulkley;

timeVaryingHerschelBulkleyCoeffs
{
    // Herschel-Bulkley parameters
    k0          0.01;               // Initial consistency coefficient [Pa.s^n]
    n           0.8;                // Flow behavior index [-]
    tau0        10;                 // Yield stress [Pa]
    nuMin       1e-06;              // Minimum kinematic viscosity [m²/s]
    nuMax       1e-01;              // Maximum kinematic viscosity [m²/s]
    
    // Time variation parameters
    timeVariationType   "power";    // Options: "power" or "exponential"
    A           0.01;               // Coefficient A [Pa.s^n]
    B           0.5;                // Exponent B [-] for power, [1/s] for exponential
}
```

### Example Applications

#### 1. Cement Grout Injection
```cpp
// Simulating cement hydration with increasing consistency
timeVariationType   "power";
A           0.001;          // Initial low consistency
B           1.2;            // Superlinear growth
```

#### 2. Polymer Degradation
```cpp
// Simulating polymer breakdown with decreasing viscosity
timeVariationType   "exponential";
A           0.1;            // Initial high consistency
B           -0.05;          // Decay rate
```



## Academic Use and Citation

This model has been developed as part of academic research in computational rheology. If you use this implementation in your research, please contant us.
cite:

```bibtex
@software{timeVaryingHB2024,
  author       = {Your Name},
  title        = {Time-Varying Herschel-Bulkley Viscosity Model for OpenFOAM},
  year         = {2024},
  publisher    = {GitHub},
  version      = {v1.0},
  doi          = none,
  url          = {https://github.com/searchinvaild/timeVaryingHerschelBulkley}
}
```

### Collaboration and Academic Integrity

As this is an active research tool, I encourage:
- 🤝 **Collaboration**: Please contact me if you're using this model for research purposes
- 📧 **Communication**: Reach out for potential collaborations or if you need assistance
- 📊 **Data Sharing**: Share your validation cases to improve the model

**Contact**: [2023010110@mail.hfut.edu.cn](mailto:2023010110@mail.hfut.edu.cn)

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- OpenFOAM Foundation for the excellent CFD framework
- [Your Institution/Lab Name] for supporting this research
- Contributors and users who provided feedback

## Related Publications



---

### Note on Academic Use

While this software is released under GPL v3.0 and can be freely used, I kindly request that researchers using this model for academic publications **consider reaching out for potential collaboration**. This helps:
- Ensure proper implementation and usage
- Foster academic collaboration
- Track the impact of this research tool
- Potentially co-author methodological papers

This is not a legal requirement but a professional courtesy in the academic community.

---

**Last Updated**: July 2025  
