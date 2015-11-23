#!/usr/bin/env python2

import sys

# check input parameters, generate INDATA list
# code placed here to avoid unneccesary vtk loading
if len(sys.argv) == 1 :
    sys.stderr.write('Using standard data set, to use other data:\n'.format(sys.argv[0]))
    sys.stderr.write('Usage: {0}  <pdb file> <connection file> |<protname> <data file>| [repeat || for multiple proteins]\n'.format(sys.argv[0]))
    PREFIX = 'data/'
    INDATA = (PREFIX + '1U3W.pdb', PREFIX + 'adh_connections.txt',
              (PREFIX + 'ADH1', PREFIX + 'adh1_conservation.txt'),
              (PREFIX + 'ADH2', PREFIX + 'adh2_conservation.txt'),
              (PREFIX + 'ADH3', PREFIX + 'adh3_conservation.txt'),
              (PREFIX + 'ADH4', PREFIX + 'adh4_conservation.txt'),
              (PREFIX + 'ADH5', PREFIX + 'adh5_conservation.txt'))

    # alternative data set
    # INDATA = (PREFIX + '1U3W.pdb', PREFIX + 'adh_connections.txt',
    #           (PREFIX + 'ADH1', PREFIX + 'adh1_scores.txt'),
    #           (PREFIX + 'ADH2', PREFIX + 'adh2_scores.txt'),
    #           (PREFIX + 'ADH3', PREFIX + 'adh3_scores.txt'),
    #           (PREFIX + 'ADH4', PREFIX + 'adh4_scores.txt'),
    #           (PREFIX + 'ADH5', PREFIX + 'adh5_scores.txt'))


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
    def __init__(self, name, data, scalars) :
        self._name = name
        self._data = data
        self._scalars = scalars

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
    def scalars(self) :
        '''Get a list of the scalars'''
        return self._scalars

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
        data, scalars = self.readData(dataFile)
        self._proteins.append(Protein(protName, data, scalars))
        
    def doMenuFileExit(self):
        sys.exit(0)

    def getScoreRange(self) :
        '''Get the score range for the protein data sets that are currently visible'''
        maxScore = 0
        # minScore will be 0
        minScore = 0
        for i in range(len(self._proteins)) :
            low, high = self._proteins[self._currentData1.get()].data.GetScalarRange()
            if maxScore == 0 or high > maxScore :
                maxScore = high
        return minScore, maxScore
    
    def initUI(self) :
        '''Initialize the UI'''
        # initialise tkinter
        self._root = Tkinter.Tk()
        self._root.title('Protein data visualisation')
        
        # tkinter ui
        # vtk
        self.renderWidget = vtkTkRenderWidget(self._root, width=1000, height=800)
        self.renderWidget.pack(expand='true', fill='both', side = Tkinter.RIGHT)
        wMain = self.renderWidget.GetRenderWindow()
        wMain.AddRenderer(self.initVtk())

        # toggle protein structure
        settingsManager = Tkinter.Frame(self._root)
        settingsManager.pack(side = Tkinter.LEFT)
        
        self._proteinStructureVisible = Tkinter.IntVar()
        self._proteinStructureVisible.set(1)
        self._toggleProteinStructureBox = Tkinter.Checkbutton(settingsManager, text = 'Show protein structure',
                                                              command = self.toggleProteinStructure,
                                                              var = self._proteinStructureVisible)
        self._toggleProteinStructureBox.pack()

        # toggle current data set
        dataManager = Tkinter.Frame(settingsManager)
        dataManager.pack()
        
        groupData1 = Tkinter.LabelFrame(dataManager, text = 'Data', padx = 5, pady = 5)
        groupData1.pack(padx = 10, pady = 10, side = Tkinter.LEFT, anchor = Tkinter.N)

        self._currentData1 = Tkinter.IntVar()
        self._currentData1.set(0)

        for i in range(len(self._proteins)) :
            Tkinter.Radiobutton(groupData1, text = self._proteins[i].name,
                                command = self.toggleProteinData,
                                var = self._currentData1,
                                value = i).pack(anchor = Tkinter.W)
            
        groupData2 = Tkinter.LabelFrame(dataManager, text='Compare with', padx = 5, pady = 5)
        groupData2.pack(padx = 10, pady = 10, side = Tkinter.RIGHT)

        self._currentData2 = Tkinter.IntVar()
        self._currentData2.set(0)
        
        for i in range(len(self._proteins)) :
            Tkinter.Radiobutton(groupData2, text = self._proteins[i].name,
                                command = self.toggleProteinData,
                                var = self._currentData2,
                                value = i).pack(anchor = Tkinter.W)

        # make sure the correct data set is shown
        self.toggleProteinData()

        # color scaling
        colorManager = Tkinter.LabelFrame(settingsManager, text = 'Color scaling', padx = 5, pady = 5)
        colorManager.pack()

        colorManagerLow = Tkinter.LabelFrame(colorManager, text = 'Lower limit', padx = 5, pady = 5)
        colorManagerLow.pack(side = Tkinter.LEFT)
        self._colorLow = Tkinter.StringVar()
        self._colorLow.set("0")
        self._colorScalerLow = Tkinter.Spinbox(colorManagerLow, from_ = 0, to = self.getScoreRange()[1],
                                               textvariable = self._colorLow, width = 8,
                                               command = self.updateColorScale, increment = 0.1)
        self._colorScalerLow.pack()

        colorManagerHigh = Tkinter.LabelFrame(colorManager, text = 'Upper limit', padx = 5, pady = 5)
        colorManagerHigh.pack(side = Tkinter.RIGHT)
        self._colorHigh = Tkinter.StringVar()
        self._colorHigh.set(str(self.getScoreRange()[1]))
        self._colorScalerHigh = Tkinter.Spinbox(colorManagerHigh, from_ = 0, to = self.getScoreRange()[1],
                                                textvariable = self._colorHigh, width = 8,
                                                command = self.updateColorScale, increment = 0.1)
        self._colorScalerHigh.pack(side = Tkinter.RIGHT)
        
        
        self._root.mainloop()

    def initVtk(self) :
        '''Initialize the VTK renderer'''
        main = vtk.vtkRenderer()
        main.SetBackground(0.2, 0.2, 0.2)
        main.AddActor(self._protStruct)
        self._dataVisualiserAtoms = self.initVtkAtoms(self._proteins[0].data)
        self._dataVisualiserBonds = self.initVtkBonds(self._proteins[0].data)
        # [0] : actor, [1] : datamanager
        main.AddActor(self._dataVisualiserAtoms[0])
        main.AddActor(self._dataVisualiserBonds[0])
        self._aSBar = self.initVtkBar()
        main.AddActor(self._aSBar)
        main.ResetCamera()

        return main
    
    def initVtkAtoms(self, data) :
        '''Initialize the residue representations'''
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
        mAtom.SetScalarRange(*data.GetScalarRange())
        
        aAtom = vtk.vtkActor()
        aAtom.SetMapper(mAtom)

        return aAtom, gAtom

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
        bond.SetRadius(0.15)
        bond.SetVaryRadiusToVaryRadiusOff()

        mBond = vtk.vtkPolyDataMapper()
        mBond.SetLookupTable(self._lut)
        mBond.SetInputConnection(bond.GetOutputPort())
        mBond.SetScalarRange(*data.GetScalarRange())

        aBond = vtk.vtkActor()
        aBond.SetMapper(mBond)
        
        return aBond, bond
    
    def initVtkColors(self) :
        self._lut = vtk.vtkLookupTable()
        self._lut.SetHueRange(0.5, 0.0)
        self._lut.SetValueRange(1.0, 1.0)
        self._lut.SetSaturationRange(1.0, 1.0)
        self._lut.SetTableRange(0.0, 1.0)
        
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

    def joinData(self, prot1, prot2) :
        newScal = list()
        for i in range(len(self._proteins[prot1].scalars)) :
            newScal.append(abs(self._proteins[prot1].scalars[i]-self._proteins[prot2].scalars[i]))

        newData = vtk.vtkPolyData()
        newData.SetPoints(self._proteins[prot1].data.GetPoints())
        newData.GetPointData().SetScalars(self.packScalars(newScal))
        newData.SetLines(self._proteins[prot1].data.GetLines())
        return newData

    def packScalars(self, scalars) :
        outScalars = vtk.vtkFloatArray()
        for i in range(len(scalars)) :
            outScalars.InsertNextValue(scalars[i])
        return outScalars
    
    def readData(self, scoreFile) :
        data = vtk.vtkPolyData()
        data.SetPoints(self._protCoord)
        scalars = read_data.read_scalars(scoreFile)
        data.GetPointData().SetScalars(self.packScalars(scalars))
        data.SetLines(self._protConnections)

        return data, scalars

    def toggleProteinData(self) :
        '''Check which proteins are active and make them visible'''
        if self._currentData1.get() == self._currentData2.get() :
            newData = self._proteins[self._currentData1.get()].data
        else :
            newData = self.joinData(self._currentData1.get(),
                                    self._currentData2.get())

        self._dataVisualiserAtoms[1].SetInputData(newData)
        self._dataVisualiserAtoms[1].Modified()
        self._dataVisualiserAtoms[1].Update()
        self._dataVisualiserBonds[1].SetInputData(newData)
        self._dataVisualiserBonds[1].Modified()
        self._dataVisualiserBonds[1].Update()

        self.renderWidget.GetRenderWindow().Render()

    def toggleProteinStructure(self) :
        '''Toggle the protein structure'''
        self._protStruct.SetVisibility(self._proteinStructureVisible.get())
        self.renderWidget.GetRenderWindow().Render()
        
    def updateVtkColors(self, minValue, maxValue) :
        '''Update the value range of the color table'''
        self._lut.SetTableRange(minValue, maxValue)
        self._dataVisualiserBonds[0].GetMapper().SetScalarRange(minValue, maxValue)
        self._dataVisualiserBonds[0].GetMapper().Modified()
        self._dataVisualiserBonds[0].GetMapper().Update()
        self._dataVisualiserAtoms[0].GetMapper().SetScalarRange(minValue, maxValue)
        self._dataVisualiserAtoms[0].GetMapper().Modified()
        self._dataVisualiserAtoms[0].GetMapper().Update()
        
    def updateColorScale(self) :
        self._colorScalerHigh.config(from_ = self._colorLow.get())
        self._colorScalerLow.config(to = self._colorHigh.get())
        if float(self._colorHigh.get()) < float(self._colorLow.get()) :
            self._colorHigh = self._colorLow.get()
        self.updateVtkColors(float(self._colorLow.get()), float(self._colorHigh.get()))
        self.renderWidget.GetRenderWindow().Render()
            
if __name__ == '__main__' :
    interface = Interface(INDATA)
