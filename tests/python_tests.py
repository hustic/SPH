import vtk
import pandas as pd

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


def test_read_keys():
        '''
        Test if the post.py is producing .h5 file correctly
        '''
        hdf = pd.HDFStore('tests/output_test.h5', mode = 'r')
        df = hdf.get('/s00001')
        assert list(df.keys()) == ['x', 'y', 'z', 'vel_x', 'vel_y', 'vel_z', 'pressure', 
                                'types','density'], "keys assert failed"
        hdf.close()
    
