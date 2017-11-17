import numpy as np
import shapefile


def create_array_from_mflist(mflist, kper, no_data=1e-9):
    """
    Method to create a numpy array of values from mflist object
    Parmaeters:
        mflist: mflist object from package (stress_period_data) object
        kper: stress period, if negative it returns all stress periods in 4d array
        no_data: fill value for nan
    Returns:
        np.array
    """
    outarrays = {}
    fm_array = mflist.array
    for key, value in fm_array.items():
        if kper >= 0:
            arr = value[kper]
        else:
            arr = value
        arr[np.isnan(arr)] = no_data
        outarrays[key] = arr

    return outarrays


def create_shapefile_from_array(shp_name, array, nlay, xgrid, ygrid, field_name=None,
                                no_data=1e-9):
    """
    Creates a shapefile from a 3d array
    Parameters:
        shp_name: name of shapefile.shp
        array: 3d numpy array, or dict of 3d np arrays
        xgrid: model x vertices
        ygrid: model y vertices
        field_name: name for shapefile field, not used if dict is supplied
        no_data: filter the nan value out
    """
    w = shapefile.Writer(shapeType=5)
    w.field('k', 'N')
    w.field('i', 'N')
    w.field('j', 'N')

    if isinstance(array, dict):
        field_name = []
        for key in array:
            if key not in ("date", "observation_name", "soil_class"):
                w.field(key, "N", decimal=8)
            else:
                w.field(key, "C")

            field_name.append(key)

    else:
        if field_name is not None:
            w.field(field_name, "N", decimal=8)
        else:
            raise AssertionError('Field name cannot be blank for non-dict arrays')

    for k in range(nlay):
        for i in range(1, xgrid.shape[0]):
            for j in range(1, xgrid.shape[1]):
                write = False
                if isinstance(field_name, list):
                    for key in field_name:
                        try:
                            if array[key][k, i-1, j-1] != no_data:
                                write = True
                        except IndexError:
                            write = False
                else:
                    try:
                        if array[k, i-1, j-1] != no_data:
                            write = True
                    except IndexError:
                        write = False

                if write:
                    w.poly([[[xgrid[i-1, j-1], ygrid[i-1, j-1]],
                             [xgrid[i, j-1], ygrid[i, j-1]],
                             [xgrid[i, j], ygrid[i, j]],
                             [xgrid[i-1, j], ygrid[i-1, j]],
                             [xgrid[i-1, j-1], ygrid[i-1, j-1]]]])
                    if isinstance(field_name, list):
                        t = []
                        for key in field_name:
                            t.append(array[key][k, i-1, j-1])

                        w.record(k+1, i, j, *t)

                    else:
                        w.record(k+1, i, j, array[k, i-1, j-1])

    w.save(shp_name)


def create_shapefile_from_transient_array(shp_name, array, nper, nlay, xgrid,
                                          ygrid, field_name=None, no_data=1e-9):
    """
    Creates a shapefile from a 4d array

    Parameters:
        shp_name: name of shapefile.shp
        array: 3d numpy array, or dict of 3d np arrays
        nper: number of model stress periods
        xgrid: model x vertices
        ygrid: model y vertices
        field_name: name for shapefile field, not used if dict is supplied
        no_data: filter the nan value out
    """

    w = shapefile.Writer(shapeType=5)
    w.field('kper', 'N')
    w.field('k', 'N')
    w.field('i', 'N')
    w.field('j', 'N')

    if isinstance(array, dict):
        field_name = []
        for key in array:
            w.field(key, "N", decimal=8)
            field_name.append(key)

    else:
        if field_name is not None:
            w.field(field_name, "N", decimal=8)
        else:
            raise AssertionError(
                'Field name cannot be blank for non-dict arrays')
    for per in range(nper):
        for k in range(nlay):
            for i in range(1, xgrid.shape[0]):
                for j in range(1, xgrid.shape[1]):
                    write = False
                    if isinstance(field_name, list):
                        for key in field_name:
                            try:
                                if array[key][per, k, i - 1, j - 1] != no_data:
                                    write = True
                            except IndexError:
                                write = False
                    else:
                        try:
                            if array[per, k, i - 1, j - 1] != no_data:
                                write = True
                        except IndexError:
                            write = False

                    if write:
                        w.poly([[[xgrid[i - 1, j - 1], ygrid[i - 1, j - 1]],
                                 [xgrid[i, j - 1], ygrid[i, j - 1]],
                                 [xgrid[i, j], ygrid[i, j]],
                                 [xgrid[i - 1, j], ygrid[i - 1, j]],
                                 [xgrid[i - 1, j - 1], ygrid[i - 1, j - 1]]]])
                        if isinstance(field_name, list):
                            t = []
                            for key in field_name:
                                t.append(array[key][per, k, i - 1, j - 1])

                            w.record(per + 1, k + 1, i, j, *t)

                        else:
                            w.record(per + 1, k + 1, i, j, array[k, i - 1, j - 1])

    w.save(shp_name)


if __name__ == "__main__":
    import flopy as fp
    ws = r'C:\Users\jlarsen\Desktop\Lucerne\Lucerne_OWHM\V0_initial_from_MODOPTIM'
    nam = 'Lucerne.nam'
    kper = 0

    ml = fp.modflow.Modflow.load(nam, model_ws=ws, check=False)
    # ml.dis.sr._reset()
    ml.dis.sr.xul = 493902.43
    ml.dis.sr.yul = 3837294.46

    nlay = ml.dis.nlay
    nper = ml.nper
    dis = ml.get_package("DIS")
    rch = ml.get_package("RCH")
    drn = ml.get_package("DRN")
    wel = ml.get_package("WEL")

    drn_array = create_array_from_mflist(drn.stress_period_data, kper)
    rch_array = rch.rech.array[kper]
    wel_array = create_array_from_mflist(wel.stress_period_data, -1)

    print('break')
    shp_name = r'C:\Users\jlarsen\Desktop\Lucerne\GIS\V0_from_MODOPTIM\wel.shp'
    create_shapefile_from_transient_array(shp_name, wel_array, nper, nlay,
                                          dis.sr.xgrid, dis.sr.ygrid)

    # shp_name = r'C:\Users\jlarsen\Desktop\Lucerne\GIS\V0_from_MODOPTIM\rch.shp'
    # create_shapefile_from_array(shp_name, rch_array, nlay, dis.sr.xgrid,
    #                            dis.sr.ygrid, field_name='rch', no_data=0)

    # shp_name = r'C:\Users\jlarsen\Desktop\Lucerne\GIS\V0_from_MODOPTIM\drn.shp'
    # create_shapefile_from_array(shp_name, drn_array, nlay, dis.sr.xgrid,
    #                            dis.sr.ygrid)



