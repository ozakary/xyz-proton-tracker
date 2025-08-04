#!/usr/bin/env python
import sys
import os
import numpy as np
import re

def extract_cell_dimensions(header_line):
    """
    Extract cell dimensions from the header line of an XYZ file.
    
    Args:
        header_line: The header line containing lattice information
        
    Returns:
        Array of [Lx, Ly, Lz] or None if not found
    """
    # Look for Lattice="Lx 0.0 0.0 0.0 Ly 0.0 0.0 0.0 Lz" pattern
    lattice_match = re.search(r'Lattice="([0-9.]+)\s+0\.0\s+0\.0\s+0\.0\s+([0-9.]+)\s+0\.0\s+0\.0\s+0\.0\s+([0-9.]+)"', header_line)
    
    if lattice_match:
        lx = float(lattice_match.group(1))
        ly = float(lattice_match.group(2))
        lz = float(lattice_match.group(3))
        return [lx, ly, lz]
    
    return None

def calculate_distance(pos1, pos2, cell_dimensions=None):
    """
    Calculate distance between two points, optionally accounting for periodic boundary conditions.
    
    Args:
        pos1, pos2: Coordinate arrays [x, y, z]
        cell_dimensions: Cell dimensions [Lx, Ly, Lz] for periodic boundary conditions
    
    Returns:
        Minimum distance considering periodic boundaries if cell_dimensions is provided
    """
    if cell_dimensions is None:
        # Standard Euclidean distance
        return np.sqrt(np.sum((np.array(pos1) - np.array(pos2))**2))
    
    # Convert to numpy arrays
    pos1 = np.array(pos1)
    pos2 = np.array(pos2)
    cell = np.array(cell_dimensions)
    
    # Calculate direct distance
    delta = pos1 - pos2
    
    # Apply minimum image convention for periodic boundaries
    for i in range(3):
        if abs(delta[i]) > 0.5 * cell[i]:
            delta[i] = delta[i] - np.sign(delta[i]) * cell[i]
    
    return np.sqrt(np.sum(delta**2))

def identify_species_in_frame(atoms, cell_dimensions, oh_bond_max_distance=1.3):
    """
    Identify OH-, H2O, and H3O+ species in a single frame.
    
    Args:
        atoms: List of atom dictionaries with 'species', 'pos', 'id'
        cell_dimensions: Cell dimensions for periodic boundary conditions
        oh_bond_max_distance: Maximum O-H bond distance to consider
        
    Returns:
        Dictionary with classification of oxygen atoms
    """
    # Separate atoms by species
    oxygens = [atom for atom in atoms if atom['species'] == 'O']
    hydrogens = [atom for atom in atoms if atom['species'] == 'H']
    
    # For each oxygen, count bonded hydrogens
    oxygen_coordination = {}
    
    for oxygen in oxygens:
        o_id = oxygen['id']
        o_pos = oxygen['pos']
        bonded_h = []
        
        for hydrogen in hydrogens:
            h_pos = hydrogen['pos']
            h_id = hydrogen['id']
            
            dist = calculate_distance(o_pos, h_pos, cell_dimensions)
            if dist <= oh_bond_max_distance:
                bonded_h.append(h_id)
        
        # Classify based on number of hydrogens
        coordination = len(bonded_h)
        oxygen_coordination[o_id] = {
            'pos': o_pos,
            'coordination': coordination,
            'bonded_h': bonded_h
        }
    
    # Count species types
    oh_minus_count = sum(1 for info in oxygen_coordination.values() if info['coordination'] == 1)
    h2o_count = sum(1 for info in oxygen_coordination.values() if info['coordination'] == 2)
    h3o_plus_count = sum(1 for info in oxygen_coordination.values() if info['coordination'] == 3)
    other_count = sum(1 for info in oxygen_coordination.values() if info['coordination'] not in [1, 2, 3])
    
    return {
        'oxygen_coordination': oxygen_coordination,
        'stats': {
            'OH-': oh_minus_count,
            'H2O': h2o_count,
            'H3O+': h3o_plus_count,
            'other': other_count
        }
    }

def process_xyz_file(input_file, output_file, oh_bond_distance=1.3):
    """
    Process an XYZ trajectory file, coloring oxygens based on protonation state.
    
    Args:
        input_file: Path to input XYZ trajectory file
        output_file: Path to output XYZ trajectory file with color embedding
        oh_bond_distance: Maximum O-H bond distance to consider
    """
    # Define colors for different species
    colors = {
        'Xe': [0.5, 0.5, 0.5],   # Grey for Xenon
        'H': [1.0, 1.0, 1.0],    # White for Hydrogen
        'OH-': [0.0, 0.0, 1.0],  # Blue for OH- oxygen
        'H2O': [1.0, 0.0, 0.0],  # Red for H2O oxygen
        'H3O+': [1.0, 0.5, 0.0]  # Orange for H3O+ oxygen
    }
    
    # Read input file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    output_lines = []
    line_idx = 0
    frame_count = 0
    
    print(f"Processing {input_file} to {output_file}...")
    
    while line_idx < len(lines):
        try:
            # Parse number of atoms
            natoms = int(lines[line_idx].strip())
            line_idx += 1
            
            # Parse header
            header = lines[line_idx].strip()
            line_idx += 1
            
            # Extract cell dimensions
            cell_dimensions = extract_cell_dimensions(header)
            if cell_dimensions is None:
                print(f"Warning: No cell dimensions found in frame {frame_count}, not using PBC")
            
            # Ensure the header has color:R:3 property
            if "color:R:3" not in header and "Properties=" in header:
                header = header.replace("id:I:1", "id:I:1:color:R:3")
            
            # Parse atoms in this frame
            atoms = []
            for i in range(natoms):
                if line_idx + i >= len(lines):
                    break
                
                parts = lines[line_idx + i].strip().split()
                if len(parts) < 5:  # Need at least species, pos, id
                    continue
                    
                species = parts[0]
                pos = [float(parts[1]), float(parts[2]), float(parts[3])]
                atom_id = int(parts[4])
                
                atoms.append({
                    'species': species,
                    'pos': pos,
                    'id': atom_id
                })
            
            # Identify species in this frame
            frame_info = identify_species_in_frame(atoms, cell_dimensions, oh_bond_distance)
            
            # Print stats for first frame and every 10th frame
            if frame_count == 0 or frame_count % 10 == 0:
                stats = frame_info['stats']
                print(f"Frame {frame_count} stats: OH-={stats['OH-']}, H2O={stats['H2O']}, H3O+={stats['H3O+']}, other={stats['other']}")
            
            # Write frame to output
            output_lines.append(f"{natoms}\n")
            output_lines.append(f"{header}\n")
            
            for atom in atoms:
                species = atom['species']
                pos = atom['pos']
                atom_id = atom['id']
                
                # Determine color
                if species == 'O':
                    coordination = frame_info['oxygen_coordination'][atom_id]['coordination']
                    if coordination == 1:
                        color = colors['OH-']
                    elif coordination == 2:
                        color = colors['H2O']
                    elif coordination == 3:
                        color = colors['H3O+']
                    else:
                        # Unusual coordination - use a distinctive color
                        color = [1.0, 0.0, 1.0]  # Magenta
                else:
                    color = colors.get(species, [0.5, 0.5, 0.5])  # Default grey
                
                # Format the output line
                color_str = ' '.join(map(str, color))
                output_line = f"{species} {pos[0]} {pos[1]} {pos[2]} {atom_id} {color_str}\n"
                output_lines.append(output_line)
            
            # Move to next frame
            line_idx += natoms
            frame_count += 1
            
        except Exception as e:
            print(f"Error processing frame at line {line_idx}: {e}")
            break
    
    # Write output file
    with open(output_file, 'w') as f:
        f.writelines(output_lines)
    
    print(f"Successfully processed {frame_count} frames")
    print("Color scheme used:")
    print(f"  OH- (1 hydrogen): Blue {colors['OH-']}")
    print(f"  H2O (2 hydrogens): Red {colors['H2O']}")
    print(f"  H3O+ (3 hydrogens): Orange {colors['H3O+']}")
    print(f"  Hydrogen: White {colors['H']}")
    print(f"  Xenon: Grey {colors['Xe']}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py input_xyz_file output_xyz_file")
        print("  This script colors oxygen atoms based on their protonation state:")
        print("  - OH- (1 hydrogen): Blue")
        print("  - H2O (2 hydrogens): Red")
        print("  - H3O+ (3 hydrogens): Orange")
        return
        
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} not found")
        return
    
    try:
        process_xyz_file(input_file, output_file)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
