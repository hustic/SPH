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
        print("okay")
        hdf = pd.HDFStore('tests/output_test.h5', mode='r')
        key = hdf.keys()[0]
        df = hdf.get(key)
        assert list(df.keys()) == ['x', 'y', 'z', 'vel_x', 'vel_y', 'vel_z', 'pressure', 
                                'types','density'], "keys assert failed"
        hdf.close()


# hard coded value (or we can read the input file here)
# minx -0.6, min y -0.6
def test_particle_behavior():
        hdf = pd.HDFStore('tests/output_test.h5', mode='r')
        keylist = sorted(hdf.keys())
        for i in keylist:
                df = hdf.get(i)
                # leaking 
                assert len(df.loc[df.x<-0.6])==0, ("particles leaking out from left outter boundary, key = " + i)
                assert len(df.loc[df.x>20.6])==0, ("particles leaking out from right outter boundary, key = " + i)
                assert len(df.loc[df.y<-0.6])==0, ("particles leaking out from bottom outter boundary, key = " + i)
                assert len(df.loc[df.y>10.6])==0, ("particles leaking out from top outter boundary, key = " + i)
                # liquid leaking from the inner boundary

                liquid = df.loc[df.types=='0']
                #print("This part would not pass the test! leaking SSSSSSSsligtly into boundary, check TMR")
                #print("Also Check the type thing!")
                assert len(liquid.loc[liquid.x<0])==0, ("liquid leaking out from left inner boundary, key = " + i)
                assert len(liquid.loc[liquid.x>20.])==0, ("liquid leaking out from right inner boundary, key = " + i)
                assert len(liquid.loc[liquid.y<0])==0, ("liquid leaking out from bottom inner boundary, key = " + i)
                assert len(liquid.loc[liquid.y>10.])==0, ("liquid leaking out from top inner boundary, key = " + i)

                boundary = df.loc[df.types=='1']
                assert len(boundary.loc[boundary.vel_x!=0.])==0, ("Moving Boundary, key = " + i)
                assert len(boundary.loc[boundary.vel_y!=0.])==0, ("Moving Boundary, key = " + i)
                #assert len(boundary.loc[boundary.y<0.])==0, ("liquid leaking out from bottom inner boundary, key = " + i)
                #assert len(boundary.loc[boundary.y>10.])==0, ("liquid leaking out from top inner boundary, key = " + i)

        hdf.close()
        
