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

    # alternative data set
    # INDATA = ('1U3W.pdb', 'adh_connections.txt',
    #           ('ADH1', 'adh1_scores.txt'),
    #           ('ADH2', 'adh2_scores.txt'),
    #           ('ADH3', 'adh3_scores.txt'),
    #           ('ADH4', 'adh4_scores.txt'),
    #           ('ADH5', 'adh5_scores.txt'))


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
        low, high = self._proteins[self._currentData1.get()].data.GetScalarRange()
        return low, high
        
    def initMenu(self) :
        menuBar = Tkinter.Menu(self._root)
        menuFile = Tkinter.Menu(menuBar, tearoff=0)
        menuFile.add_command(label='Exit', command=self.doMenuFileExit)
        menuBar.add_cascade(label='File', menu=menuFile)

        return menuBar
    
    def initUI(self) :
        '''Initialize the UI'''
        # initialise tkinter
        self._root = Tkinter.Tk()
        self._root.config(menu = self.initMenu())
        self._root.title('ADH residue uniqueness')
        
        # tkinter ui
        # vtk
        self.renderWidget = vtkTkRenderWidget(self._root, width=800, height=800)
        self.renderWidget.pack(expand='true', fill='both', side = Tkinter.RIGHT)
        wMain = self.renderWidget.GetRenderWindow()
        wMain.AddRenderer(self.initVtk())

        # toggle protein structure
        settingsManager = Tkinter.Frame(self._root)
        settingsManager.pack(side = Tkinter.LEFT)
        
        self._proteinStructureVisible = Tkinter.IntVar()
        self._proteinStructureVisible.set(1)
        self._toggleProteinStructureBox = Tkinter.Checkbutton(settingsManager, text = 'Show protein structure',
                                                              command = self.toggleProteinStructure, state = 'active',
                                                              var = self._proteinStructureVisible)
        self._toggleProteinStructureBox.pack()

        # toggle current data set
        dataManager = Tkinter.Frame(settingsManager)
        dataManager.pack()
        
        groupData1 = Tkinter.LabelFrame(dataManager, text='Data', padx = 5, pady = 5)
        groupData1.pack(padx = 10, pady = 10, side = Tkinter.LEFT, anchor = Tkinter.N)

        self._currentData1 = Tkinter.IntVar()
        self._currentData1.set(0)

        for i in range(len(self._proteins)) :
            Tkinter.Radiobutton(groupData1, text = self._proteins[i].name,
                                command = self.toggleProteinData,
                                var = self._currentData1,
                                value = i).pack(anchor = Tkinter.W)
            # make sure the correct data is shown from start
            self.toggleProteinData()
            
        groupData2 = Tkinter.LabelFrame(dataManager, text='Compare with', padx = 5, pady = 5)
        groupData2.pack(padx = 10, pady = 10, side = Tkinter.RIGHT)

        self._currentData2 = Tkinter.IntVar()
        self._currentData2.set(len(self._proteins) + 1)
        
        for i in range(len(self._proteins)) :
            Tkinter.Radiobutton(groupData2, text = self._proteins[i].name,
                                command = self.toggleProteinData,
                                var = self._currentData1,
                                value = i).pack(anchor = Tkinter.W)
            # make sure the correct data is shown from start
            self.toggleProteinData()
        Tkinter.Radiobutton(groupData2, text = 'None',
                            command = self.toggleProteinData,
                            var = self._currentData1,
                            value = len(self._proteins) + 1).pack(anchor = Tkinter.W)
            
        self._root.mainloop()
    
    def initVtk(self) :
        '''Initialize the VTK renderer'''
        main = vtk.vtkRenderer()
        main.SetBackground(0.2, 0.2, 0.2)
        main.AddActor(self._protStruct)
        for i in range(len(self._proteins)) :
            main.AddActor(self._proteins[i].atoms)
            main.AddActor(self._proteins[i].bonds)
            self.updateVtkColors(*self._proteins[i].data.GetScalarRange())
        self._aSBar = self.initVtkBar()
        main.AddActor(self._aSBar)
        main.ResetCamera()

        return main
    
    def initVtkAtoms(self, data) :
        '''Initialize the residue representations'''
        sAtom = vtk.vtkSphereSource()
        sAtom.SetRadius(0.5)# + len(self._proteins)/10.0)
        sAtom.SetThetaResolution(15)
        sAtom.SetPhiResolution(15)
        
        gAtom = vtk.vtkGlyph3D()
        gAtom.SetInputData(data)
        gAtom.SetSourceConnection(sAtom.GetOutputPort())
        gAtom.SetScaleModeToDataScalingOff()
    
        mAtom = vtk.vtkPolyDataMapper()
        mAtom.SetInputConnection(gAtom.GetOutputPort())
        mAtom.SetLookupTable(self._lut)
        mAtom.SetScalarRange(*data.GetScalarRange())
        
        aAtom = vtk.vtkActor()
#        aAtom.GetProperty().SetSpecular(0.1)
#        aAtom.GetProperty().SetOpacity(0.5)
        aAtom.SetMapper(mAtom)

        return aAtom

    def initVtkBar(self) :
        '''Initialize the bond connectors'''
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
        bond.SetRadius(0.15) # + len(self._proteins)/20.0)
        bond.SetVaryRadiusToVaryRadiusOff()

        mBond = vtk.vtkPolyDataMapper()
        mBond.SetLookupTable(self._lut)
        mBond.SetInputConnection(bond.GetOutputPort())
        mBond.SetScalarRange(*data.GetScalarRange())

        aBond = vtk.vtkActor()
        aBond.SetMapper(mBond)
#        aBond.GetProperty().SetSpecular(0.1)
#        aBond.GetProperty().SetOpacity(0.5)
        
        return aBond
    
    def initVtkColors(self) :
        self._lut = vtk.vtkLookupTable()
        self._lut.SetHueRange(0.5, 0.0)
        self._lut.SetValueRange(1.0, 1.0)
        self._lut.SetSaturationRange(1.0, 1.0)
        self._lut.SetTableRange(0.0, 1.0)
#        self._lut.SetScaleToLog10()
        
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

    def joinData(self, set1, set2) :
        scal1 = data1.GetPoints().GetScalars()
        scal2 = data2.GetPoints().GetScalars()
        newScal = abs(scal1 - scal2)
        newData = vtk.vtkPolyData()
        newData.SetPoints(scal1.GetPoints())
        newData.GetPoints().SetScalars(newScal)
        return newData
    
    def readData(self, scoreFile) :
        data = vtk.vtkPolyData()
        data.SetPoints(self._protCoord)
        data.GetPointData().SetScalars(read_data.read_scalars(scoreFile))
        data.SetLines(self._protConnections)

        return data

    def toggleProteinData(self) :
        '''Check which proteins are active and make them visible, and request color table update'''
        # set color scale            
        # upgrade visibility
        for i in range(len(self._proteins)) :
            if self._currentData1.get() == i :
                self._proteins[i].atoms.SetVisibility(1)
                self._proteins[i].bonds.SetVisibility(1)
            else :
                self._proteins[i].atoms.SetVisibility(0)
                self._proteins[i].bonds.SetVisibility(0)

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
