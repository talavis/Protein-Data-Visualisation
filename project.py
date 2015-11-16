#!/usr/bin/env python2

import sys

# check input parameters, generate INDATA list
# code placed here to avoid unneccesary vtk loading
if len(sys.argv) == 1 :
    sys.stderr.write('Using standard data set, to use other data:\n'.format(sys.argv[0]))
    sys.stderr.write('Usage: {0}  <pdb file> <connection file> |<protname> <data file>| [repeat || for multiple proteins]\n'.format(sys.argv[0]))
    INDATA = ('1U3W.pdb', 'adh_connections.txt',
              ('ADH1', 'adh1_conservation.txt'),
              ('ADH2', 'adh2_conservation.txt'),
              ('ADH3', 'adh3_conservation.txt'),
              ('ADH4', 'adh4_conservation.txt'),
              ('ADH5', 'adh5_conservation.txt'))
elif len(sys.argv) > 3 and len(sys.argv) % 2 != 1 :
    sys.stderr.write('Usage: {0} <pdb file> <connection file> |<protname> <data file>| [repeat || for multiple proteins]\n'.format(sys.argv[0]))
    sys.exit()
else :
    INDATA = list()
    INDATA.append(sys.argv[1])
    INDATA.append(sys.argv[2])
    for i in range(3, len(sys.argv), 2) :
        INDATA.append((sys.argv[i], sys.argv[i+1]))

import Tkinter
import vtk
from vtk.tk.vtkTkRenderWidget import vtkTkRenderWidget

import read_data

class Protein :
    '''Stores information about a protein'''
    def __init__(self, name, data, aAtoms, aBonds) :
        self._name = name
        self._data = data
        self._atoms = aAtoms
        self._bonds = aBonds

    @property
    def name(self) :
        '''Get the protein name'''
        return self._name
    @name.setter
    def name(self, newName) :
        '''Set the protein name'''
        self._name = newName

    @property
    def data(self) :
        '''Get the data set'''
        return self._data

    @property
    def atoms(self) :
        '''Get the residue actor'''
        return self._atoms

    @property
    def bonds(self) :
        '''Get the bond actor'''
        return self._bonds
    
class Interface :
    def __init__(self, indata) :
        self._proteins = list()
        self._protStruct = self.initVtkProtein(indata[0])
        self._protCoord = read_data.read_points(indata[0])
        self._protConnections = read_data.read_connections(indata[1])

        self.initVtkColors()
        for i in range(2, len(indata)) :
            self.addProtein(*indata[i])
        self.initUI()

    def addProtein(self, protName, dataFile) :
        '''Add a protein'''
        data = self.readData(dataFile)
        aAtoms = self.initVtkAtoms(data)
        aBonds = self.initVtkBonds(data)
        self._proteins.append(Protein(protName, data, aAtoms, aBonds))
        
    def doMenuFileExit(self):
        sys.exit(0)

    def getScoreRange(self) :
        '''Get the score range for the protein data sets that are currently visible'''
        maxScore = 0
        minScore = 0
        for i in range(len(self._toggleDataVars)) :
            if self._toggleDataVars[i].get() > 0 :
                low, high = self._proteins[i].data.GetScalarRange()
                if minScore == 0 or low < minScore :
                    minScore = low
                if high > maxScore :
                    maxScore = high
        return minScore, maxScore
        
    def initMenu(self) :
        menuBar = Tkinter.Menu(self.root)
        menuFile = Tkinter.Menu(menuBar, tearoff=0)
        menuFile.add_command(label='Exit', command=self.doMenuFileExit)
        menuBar.add_cascade(label='File', menu=menuFile)

        return menuBar
    
    def initUI(self) :
        # initialise tkinter
        self.root = Tkinter.Tk()
        self.root.config(menu = self.initMenu())
        self.root.title('ADH residue uniqueness')
        
        # vtk
        self.renderWidget = vtkTkRenderWidget(self.root, width=800, height=800)
        self.renderWidget.pack(expand='true', fill='both')        
        wMain = self.renderWidget.GetRenderWindow()
        wMain.AddRenderer(self.initVtk())

        # tkinter ui
        # protein toggle
        self._proteinStructureVisible = Tkinter.IntVar()
        self._proteinStructureVisible.set(1)
        self._toggleProteinStructureBox = Tkinter.Checkbutton(self.root, text = 'Show protein structure',
                                                              command = self.toggleProteinStructure, state = 'active',
                                                              var = self._proteinStructureVisible)
        self._toggleProteinStructureBox.pack()
        self._toggleDataVars = list()
        self._toggleDataBoxes = list()
        for i in range(len(self._proteins)) :
            varData = Tkinter.IntVar()
            varData.set(0)
            # I've been unable to find a way to only update current protein
            # so the function will check the status of all proteins (datasets)
            boxData = Tkinter.Checkbutton(self.root, text = 'Show data of {protname}'.format(protname = self._proteins[i].name),
                                          command = self.toggleProteinData, state = 'active',
                                          var = varData)
            boxData.pack()
            self._toggleDataVars.append(varData)
            self._toggleDataBoxes.append(boxData)
            # make sure the checkboxes match the data shown from start
            self.toggleProteinData()
        
        self.root.mainloop()
    
    def initVtk(self) :
        main = vtk.vtkRenderer()
        main.SetBackground(0.2, 0.2, 0.2)
        main.AddActor(self._protStruct)
        for i in range(len(self._proteins)) :
            main.AddActor(self._proteins[i].atoms)
            main.AddActor(self._proteins[i].bonds)
            self.updateVtkColors(*self._proteins[i].data.GetScalarRange())
#        self._aSBar = self.initVtkBar()
#        main.AddActor(self._aSBar)
        main.ResetCamera()

        return main
    
    def initVtkAtoms(self, data) :
        sAtom = vtk.vtkSphereSource()
        sAtom.SetRadius(0.5)
        sAtom.SetThetaResolution(15)
        sAtom.SetPhiResolution(15)
        
        gAtom = vtk.vtkGlyph3D()
        gAtom.SetInputData(data)
        gAtom.SetSourceConnection(sAtom.GetOutputPort())
        gAtom.SetScaleModeToDataScalingOff()
    
        mAtom = vtk.vtkPolyDataMapper()
        mAtom.SetInputConnection(gAtom.GetOutputPort())
        mAtom.SetLookupTable(self._lut)
        
        aAtom = vtk.vtkActor()
        aAtom.SetMapper(mAtom)

        return aAtom

    def initVtkBar(self) :
        aSBar = vtk.vtkScalarBarActor()
        aSBar.SetOrientationToVertical()
        aSBar.SetLookupTable(self._lut)
        aSBar.SetTitle("Residue score")
        aSBar.GetLabelTextProperty().SetColor(0.8, 0.8, 0.8)
        aSBar.GetTitleTextProperty().SetColor(0.8, 0.8, 0.8)
        aSBar.SetWidth(0.1)
        aSBar.SetHeight(0.9)
        spc = aSBar.GetPositionCoordinate()
        spc.SetCoordinateSystemToNormalizedViewport()
        spc.SetValue(0.05, 0.05)

        return aSBar
    
    def initVtkBonds(self, data) :
        bond = vtk.vtkTubeFilter()
        bond.SetNumberOfSides(6)
        bond.SetInputData(data)
        bond.SetRadius(0.15)
        bond.SetVaryRadiusToVaryRadiusOff()

        mBond = vtk.vtkPolyDataMapper()
        mBond.SetLookupTable(self._lut)
        mBond.SetInputConnection(bond.GetOutputPort())

        aBond = vtk.vtkActor()
        aBond.SetMapper(mBond)
        aBond.GetProperty().SetSpecular(0.1)
        
        return aBond
    
    def initVtkColors(self) :
        self._lut = vtk.vtkLookupTable()
        self._lut.SetHueRange(0.667, 0.0)
        self._lut.SetValueRange(1.0, 1.0)
        self._lut.SetSaturationRange(1.0, 1.0)
        self._lut.SetTableRange(0.0, 1.0)
        self._lut.SetScaleToLog10()
        
    def initVtkProtein(self, pdbFile) :
        reader = vtk.vtkPDBReader()
        reader.SetFileName(pdbFile)
        
        protein = vtk.vtkProteinRibbonFilter()
        protein.SetInputConnection(reader.GetOutputPort())
        protein.SetDrawSmallMoleculesAsSpheres(False)
        protein.Update()
        mProtein = vtk.vtkPolyDataMapper()
        mProtein.SetInputData(protein.GetOutput())
        mProtein.Update()
        aProtein = vtk.vtkActor()
        aProtein.SetMapper(mProtein)
        
        return aProtein

    def readData(self, scoreFile) :
        data = vtk.vtkPolyData()
        data.SetPoints(self._protCoord)
        data.GetPointData().SetScalars(read_data.read_scalars(scoreFile))
        data.SetLines(self._protConnections)

        return data

    def toggleProteinData(self) :
        '''Check which proteins are active and make them visible, and request color table update'''
        # upgrade visibility
        for i in range(len(self._toggleDataVars)) :
            self._proteins[i].atoms.SetVisibility(self._toggleDataVars[i].get())
            self._proteins[i].bonds.SetVisibility(self._toggleDataVars[i].get())

        # also update the color bar
#        self.updateVtkColors(*self.getScoreRange())

        self.renderWidget.GetRenderWindow().Render()

    def toggleProteinStructure(self) :
        '''Toggle the protein structure'''
        self._protStruct.SetVisibility(self._proteinStructureVisible.get())
        self.renderWidget.GetRenderWindow().Render()
        
    def updateVtkColors(self, minValue, maxValue) :
        '''Update the value range of the color table'''
        self._lut.SetTableRange(minValue, maxValue)
    
if __name__ == '__main__' :
    interface = Interface(INDATA)
