import vtk

def test_file_writer_output():
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName('tests/test_file_writer.vtp')
    reader.Update()

    pdata = reader.GetOutput()
    assert pdata.GetNumberOfCells() == 10
    assert pdata.GetNumberOfPoints() == 10
    assert pdata.GetPoint(5) == (5.0, 5.0, 0.0)
    assert (pdata.GetPointData().GetArray('Velocity').GetTuple(5)
            == (5.0, -5.0, 0.0))
    assert (pdata.GetPointData().GetArray('Pressure').GetValue(5)
            == 5.0)
    
