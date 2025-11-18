"""
I/O utilities for reading and writing 4D array data
"""

import numpy as np


def read4D(filename, data, dim1, dim2, dim3, dim4):
    """
    Read 4D array data from a file.
    
    Args:
        filename: Path to the input file
        data: numpy array to fill with data
        dim1, dim2, dim3, dim4: Dimensions of the array
        
    Returns:
        -1: File open error
         0: Parsing error
         1: Success
    """
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except IOError:
        print(f"Error opening file: {filename}")
        return -1
    
    if dim1 == 0 or dim2 == 0 or dim3 == 0 or dim4 == 0:
        print("Error: Dimensions must be greater than zero.")
        return 0
    
    # Read dimension line
    first_line = lines[0].strip().split()
    try:
        file_dim2, file_dim3, file_dim4 = int(first_line[0]), int(first_line[1]), int(first_line[2])
    except (ValueError, IndexError):
        print(f"Error: First line of {filename} does not contain three dimension integers.")
        return 0
    
    if file_dim2 != dim2 or file_dim3 != dim3 or file_dim4 != dim4:
        print(f"Error: Dimension mismatch in {filename}. Expected ({dim2}, {dim3}, {dim4}), got ({file_dim2}, {file_dim3}, {file_dim4}).")
        return 0
    
    # Read data lines
    for line in lines[1:]:
        tokens = line.strip().split()
        if len(tokens) < 3:
            continue
            
        try:
            i, j, k = int(tokens[0]) - 1, int(tokens[1]) - 1, int(tokens[2]) - 1
            
            col = 0
            for value_str in tokens[3:]:
                if col >= dim1:
                    break
                    
                value = float(value_str)
                idx = ((col * dim2 + i) * dim3 + j) * dim4 + k
                
                if idx >= data.size:
                    print(f"Index out of bounds while reading {filename}: {idx}({col}, {i}, {j}, {k}) for data size {data.size}")
                    return 0
                
                data.flat[idx] = value
                col += 1
                
        except (ValueError, IndexError) as e:
            print(f"Error parsing line: {line}")
            print(f"Exception: {e}")
            return 0
    
    return 1


def write4D(filename, data, dim1, dim2, dim3, dim4, vtk=False):
    """
    Write 4D array data to a file.
    
    Args:
        filename: Path to the output file
        data: numpy array containing data
        dim1, dim2, dim3, dim4: Dimensions of the array
        vtk: If True, write in VTK format (not implemented yet)
    """
    print(f"Writing data to {filename}")
    
    if vtk:
        print("VTK format is not implemented yet.")
        return
    
    with open(filename, 'w') as f:
        # Write dimensions in the first line
        f.write(f"{dim2:6d} {dim3:5d} {dim4:5d}\n")
        
        # Write data in the requested format
        for i in range(dim2):
            for j in range(dim3):
                for k in range(dim4):
                    f.write(f"{i+1:6d} {j+1:5d} {k+1:5d} ")
                    for p in range(dim1):
                        index = ((p * dim2 + i) * dim3 + j) * dim4 + k
                        f.write(f"{data.flat[index]:14.7E} ")
                    f.write("\n")
