# XYZ Protonation State Analyzer

A Python tool for analyzing molecular dynamics trajectories in XYZ format, specifically designed to identify and visualize different protonation states of oxygen atoms (OH<sup>-</sup>, H<sub>2</sub>O, H<sub>3</sub>O<sup>+</sup>) through color coding for [OVITO](https://www.ovito.org/) visualization.

## Overview

This tool processes XYZ trajectory files issued from molecular dynamics simulations and automatically identifies the protonation state of oxygen atoms based on their hydrogen coordination. Each oxygen atom is assigned a specific RGB color based on its protonation state, enabling easy visualization and analysis in [OVITO](https://www.ovito.org/).

## Features

- **Automatic Protonation State Detection**: Identifies OH<sup>-</sup>, H<sub>2</sub>O, and H<sub>3</sub>O<sup>+</sup> species based on O-H bond distances
- **Periodic Boundary Conditions**: Supports periodic boundary conditions for distance calculations
- **Color Coding**: Assigns distinct RGB colors to different protonation states
- **OVITO Compatible**: Outputs XYZ files with embedded color information readable by [OVITO](https://www.ovito.org/)
- **Progress Tracking**: Reports statistics for each frame processed
- **Flexible Parameters**: Customizable O-H bond distance threshold

## Color Scheme

| Species | Description | RGB Color | Visual |
|---------|-------------|-----------|---------|
| OH<sup>-</sup> | Hydroxide ion (1 hydrogen) | Blue (0.0, 0.0, 1.0) | 🔵 |
| H<sub>2</sub>O | Water molecule (2 hydrogens) | Red (1.0, 0.0, 0.0) | 🔴 |
| H<sub>3</sub>O<sup>+</sup> | Hydronium ion (3 hydrogens) | Orange (1.0, 0.5, 0.0) | 🟠 |
| H | Hydrogen atoms | White (1.0, 1.0, 1.0) | ⚪ |
| Xe | Xenon atoms | Grey (0.5, 0.5, 0.5) | ⚫ |

## Installation

### Prerequisites

- Python 3.6 or higher
- NumPy

### Setup

1. Clone this repository:
```bash
git clone https://github.com/ozakary/xyz-proton-tracker.git
cd xyz-proton-tracker
```

2. Install required dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

```bash
python protonation_analyzer.py input_trajectory.xyz output_colored.xyz
```

### Example

```bash
python protonation_analyzer.py water_simulation_in.xyz water_simulation_out.xyz
```

This will:
1. Read the input XYZ trajectory file
2. Analyze each frame to identify protonation states
3. Output a new XYZ file with color information embedded
4. Display progress and statistics during processing

### Expected Output

```
Processing water_simulation.xyz to water_colored.xyz...
Frame 0 stats: OH-=5, H2O=128, H3O+=3, other=0
Frame 10 stats: OH-=7, H2O=126, H3O+=3, other=0
Frame 20 stats: OH-=4, H2O=130, H3O+=2, other=0
...
Successfully processed 100 frames
Color scheme used:
  OH- (1 hydrogen): Blue [0.0, 0.0, 1.0]
  H2O (2 hydrogens): Red [1.0, 0.0, 0.0]
  H3O+ (3 hydrogens): Orange [1.0, 0.5, 0.0]
  Hydrogen: White [1.0, 1.0, 1.0]
  Xenon: Grey [0.5, 0.5, 0.5]
```

## Input File Format

The script expects XYZ trajectory files with the following format:

```
number_of_atoms
Properties=species:S:1:pos:R:3:id:I:1 Lattice="Lx 0.0 0.0 0.0 Ly 0.0 0.0 0.0 Lz"
species x y z atom_id
species x y z atom_id
...
```

Where:
- `species`: Atom type (O, H, Xe, etc.)
- `x y z`: Cartesian coordinates
- `atom_id`: Unique atom identifier
- `Lattice`: Periodic boundary conditions (optional)

## Output File Format

The output XYZ file includes color information:

```
number_of_atoms
Properties=species:S:1:pos:R:3:id:I:1:color:R:3 Lattice="Lx 0.0 0.0 0.0 Ly 0.0 0.0 0.0 Lz"
species x y z atom_id r g b
species x y z atom_id r g b
...
```

## Visualization in [OVITO](https://www.ovito.org/)

1. Open [OVITO](https://www.ovito.org/)
2. Load the output colored XYZ file
3. The color information will be automatically recognized
4. Use the "Color coding" modifier if needed to apply the embedded colors
5. Different protonation states will be clearly visible with distinct colors

## Algorithm Details

### Protonation State Identification

1. **Distance Calculation**: For each oxygen atom, calculate distances to all hydrogen atoms
2. **Bond Detection**: Identify O-H bonds using a cutoff distance (default: 1.3 Å)
3. **Coordination Counting**: Count the number of bonded hydrogens per oxygen
4. **Classification**:
   - 1 hydrogen --> OH<sup>-</sup> (hydroxide)
   - 2 hydrogens --> H<sub>2</sub>O (water)
   - 3 hydrogens --> H<sub>3</sub>O<sup>+</sup> (hydronium)

### Periodic Boundary Conditions

When cell dimensions are available, the script uses the minimum image convention:
- Calculates the shortest distance between atoms considering periodic images
- Essential for accurate bond detection in periodic systems

## Parameters

### Customizable Parameters

- **O-H Bond Distance**: Default 1.3 Å (can be modified in the script)
- **Color Scheme**: RGB values can be customized in the `colors` dictionary

### File Structure

```
xyz-proton-tracker/
├── protonation_analyzer.py    # Main script
├── requirements.txt           # Python dependencies
├── README.md                 # This file
├── examples/                 # Example input/output files
│   ├── sample_input.xyz
│   └── sample_output.xyz
└── tests/                    # Test cases
    └── test_analyzer.py
```

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-feature`)
3. Commit your changes (`git commit -am 'Add new feature'`)
4. Push to the branch (`git push origin feature/new-feature`)
5. Create a Pull Request

### Common Issues

1. **"No cell dimensions found"**: The script will still work but won't use periodic boundary conditions
2. **Memory issues with large trajectories**: Process in chunks if needed
3. **Incorrect bond detection**: Adjust the O-H distance threshold

### Contact

For questions, issues, or contributions, please open an issue on GitHub or contact [ouail.zakary@oulu.fi].
