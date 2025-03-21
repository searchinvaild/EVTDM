
# timeSlurryPower Non-Newtonian Viscosity Model for OpenFOAM 

[![OpenFOAM-v2406](https://img.shields.io/badge/OpenFOAM-2406-blue)](https://www.openfoam.com/news/main-news/openfoam-v24-06)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

![图片1](https://github.com/user-attachments/assets/988ac58f-a74f-41aa-ba20-ef9bb8aa9afc)




An enhanced viscosity model for simulating time-dependent non-Newtonian slurry transport in pipeline systems, particularly optimized for grouting process simulations with inflatable bag sealing.

## Key Features
✅ **Time-dependent rheology**  
`τ = τ₀ + k·t^α·(γ̇)^n`  
Combines time evolution (e.g., cement hydration) with shear-thinning/thickening behavior

✅ **Industrial applications**  
- Grout migration in pipe sealing  
- TBM shield tail void filling  
- Geotechnical permeation grouting

✅ **Validation-ready**  
- Built-in verification cases (Hagen-Poiseuille, Backward Facing Step)  
- Experimental comparison module

## Installation
### Requirements
- OpenFOAM v2012 or later
- GNU Make 4.0+

### Compilation
```bash
# Clone repository
git clone https://github.com/yourusername/timeSlurryPower.git
cd timeSlurryPower

# Set OpenFOAM environment
source /opt/openfoam/etc/bashrc

# Compile 
wmake
```

## Usage
### Model Activation
Add to `constant/transportProperties`:
```cpp
transportModel  timeSlurryPower;

timeSlurryPowerCoeffs
{
    k        0.05 [ 0 2 -1 0 0 0 0 ];  // Consistency index
    n        0.4;                       // Power-law index 
    tau0     10 [ 0 2 -2 0 0 0 0 ];     // Yield stress
    timeCoeff 0.2;                       // Time exponent
    nuMax    100 [ 0 2 -1 0 0 0 0 ];     // Max viscosity
}
```

### Key Parameters
| Parameter  | Units       | Physical Meaning                     |
|------------|-------------|---------------------------------------|
| k          | Pa·sⁿ       | Consistency index                    |
| n          | -           | Shear-rate sensitivity               |  
| τ₀         | Pa          | Yield stress                         |
| α          | -           | Time evolution rate                  |
| ν_max      | m²/s        | Viscosity upper limit                |

## Validation Cases
### 1. Transient Pipe Flow
```bash
# Run verification case
cd validation/pipeFlow
./Allrun
```


### 2. Backward Facing Step 
```bash 
cd validation/backwardStep
./Allrun -parallel 4
```


## Contributing
We welcome contributions through:
- [Bug reports](https://github.com/yourusername/timeSlurryPower/issues)
- Feature requests
- Validation case additions  
See [CONTRIBUTING.md](CONTRIBUTING.md) for development guidelines.

## Citation
If used in academic work, please cite:
```bibtex
@software{YourName_2024_timeSlurryPower,
  author = {Your Name},
  title = {timeSlurryPower: A Time-dependent Non-Newtonian Model for Grout Migration},
  year = {2024},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/yourusername/timeSlurryPower}}
}
```

## License
GNU GPLv3 © 2024 Your Name  
This project is not affiliated with the OpenFOAM Foundation.
```

---

**Suggested Repository Structure**

timeSlurryPower/
├── Make/
├── timeSlurryPower.C
├── timeSlurryPower.H
├── validation/
│   ├── pipeFlow/
│   │   ├── 0/ 
│   │   ├── system/
│   │   └── Allrun
│   └── backwardStep/
│       ├── 0.orig/
│       ├── system/
│       └── Allrun
├── docs/
│   ├── velocity_profile.gif
│   ├── pipe_comparison.png
│   └── bfs_streamlines.png
├── CONTRIBUTING.md
└── README.md

Let me know if you need help creating the validation case files!
