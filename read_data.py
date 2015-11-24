import vtk
import string

# pdb file:
# chain: [21]
# atom name: [13:16]
# resnr: [23:26]

# pdb cols:
# 0 ATOM
# 1 atomnr
# 2 atomname
# 3 resname
# 4 chain
# 5 resnr
# 6 x
# 7 y
# 8 z

def read_points(fileName):
    '''Read chain A Calpha coordinates from a PDB file'''
    points = vtk.vtkPoints()
    with open(fileName) as f :
        for line in f :
            data = line.split()
            if not data or data[0] != 'ATOM' or data[2] != 'CA' :
                continue
            #or data[4] != 'A' ::
            points.InsertNextPoint(float(data[6]), float(data[7]), float(data[8]))
    return points

def read_scalars(fileName):
    '''Read scores'''
    scalars = list()
    text_file = open(fileName)
    line = text_file.readline()
    while line:
        data = string.split(line)
        if data and data[0] != '#':
            x = float(data[0])
            scalars.append(x)
        line = text_file.readline()
    text_file.close()
    return scalars

def read_connections(filename):
    """Reads molecule connections from an ASCII file."""
    connections = vtk.vtkCellArray()
    text_file = open(filename)
    line = text_file.readline()
    while line:
        data = string.split(line)
        if data and data[0] != '#':
            a, b = int(data[0]), int(data[1])
            connections.InsertNextCell(2)
            connections.InsertCellPoint(a)
            connections.InsertCellPoint(b)
        line = text_file.readline()
    text_file.close()
    return connections
