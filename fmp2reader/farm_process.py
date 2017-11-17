"""
__author__ = "Joshua D Larsen"
__version__ = "0.1"

Utility to read in FMP files from modflow OWHM 2.
import into your own script as:

>>> from fmp2reader import Fmp
>>>
>>> # Read in fmp data and set to Fmp as an overridden dictionary
>>> fmp = Fmp('lucerne.fmp', fmp_ws=r"C:/mydirectory", nrow=50, ncol=64,
>>>           crop_lut={1:'alfalfa', 2:'apples'}, cell_area=4000000)
>>> # Build fmp data into static and transient arrays
>>> fmp.to_arrays()
>>>

Please do not edit code after if __name__ == "__main__" this is being used
for debugging and development purposes
"""

import ArrayToShapeFile as ATS
import numpy as np
import os


class FmpLinefeed(object):
    """
    Class to read in mnw2 linefeed objects.

    Args:
        linefeed: (str, list) linefeed file name or list of linefeed file names
        ws: (str)
    """
    def __init__(self, linefeed, ws=""):

        if isinstance(linefeed, list):
            self.num = len(linefeed)
        else:
            self.num = 1
            linefeed = [linefeed]

        self.__ws = ws
        self.name = linefeed
        self.file_array = {}
        self.combined_array = {}
        self.__reader()

    def __reader(self):
        """
        Reader method for linefeed files
        """
        for name in self.name:
            fpath = os.path.join(self.__ws, name)

            with open(fpath) as lnf:
                wells = []
                read_data = False
                arr = []
                for idx, line in enumerate(lnf):
                    if line.startswith('#'):
                        pass
                    elif 'temporal' in line.lower():
                        read_data = True
                    elif read_data:
                        arr.append([int(i) for i in line.strip().split()])
                    else:
                        t = line.strip().split()
                        well = t[0]
                        l = int(t[1])
                        r = int(t[2])
                        c = int(t[3])
                        qmax = float(t[4])
                        wells.append([well, l, r, c, qmax])

            self.file_array[name] = {'wells': wells,
                                     'data': np.array(arr, dtype=np.int)}
        self.__create_single_array()

    def __create_single_array(self):
        """
        Method to combine the file_array dictionary to create a single
        array dictionary
        """
        wells = []
        for idx, value in sorted(self.file_array.items()):
            wells += [wl for wl in value['wells']]
            nper = len(value['data'])

        data = np.zeros((nper, len(wells)))

        ix = 0
        for idx, value in sorted(self.file_array.items()):
            ix1 = ix + len(value['wells'])
            data[:, ix:ix1] = value['data']
            ix += len(value['wells'])

        self.combined_array = {"wells": wells, "data": data}


class Fmp(dict):
    """
    Method to read modflow-owhm2 FMP packages. Useful for generatating arrays
    and creating shapefiles by leveraging ATS methods.

    Parameters:
        fmp_name: (str) name of fmp file
        ws: (str) base directory of fmp file
        nrow: (int) optional, but required to create crop arrays
        ncol: (int) optional, but required to create crop arrays
        nfarm: (int) optional, read from fmp file on load
        nirrigate: (int) optional, read from fmp file on load
        load: (bool) True is to load a fmp file, False for future capability
            which is not implemented. False is where a user should be able to
            build a new fmp instance using external data
        cell_area: (float or np.ndarray) optional. Used to calculate crop area arrays
            numpy array of cell_area must be the shape (nrow, ncol)
        crop_lut: (dict) optional. LUT of crop types tied to crop number. Used
            in creation of crop arrays.
    """
    # todo: add to block names or valid_parameters to extend this module!
    block_names = ('output', 'options', 'global dimension', 'wbs', 'fwell',
                   'soil', 'climate', 'crop', 'surfacewater', 'swo')

    valid_parameters = {'options': ('noprint', 'wellfield', 'recomp_q_bd',
                                    'mnwclose'),
                        'output': ('wel_cbc', 'fnr_cbc', 'isdpfl', 'fds',
                                   'ifbpfl', 'fb_compact.out', 'fb_details',
                                   'et_array', 'et_list.out',
                                   'routing_information'),
                         'global dimension': ('nfarm', 'ncrop', 'nsoil',
                                              'nirrigate', 'nrow', 'ncol',
                                              'surface_elevation', 'nrd_types'),
                         'wbs': ('location', 'efficiency',
                                 'deficiency_scenario', 'prorate_deficiency',
                                 'bare_runoff_frac'),
                         'fwell': ('qmaxreset', 'nocirnoq', 'print',
                                   'linefeed'),
                         'soil': ('soil_id', 'capillary_fringe'),
                         'climate': ('reference_et', 'precipitation'),
                         'crop': ('print', 'fraction', 'crop_coeficient',
                                  'irrigation', 'root_depth', 'pond_depth',
                                  'transpiration_fraction',
                                  'evaporation_irrigation_fraction',
                                  'surfacewater_loss_fraction_precipitation',
                                  'surfacewater_loss_fraction_irrigation',
                                  'groundwater_root_interaction'),
                        'surfacewater': ('non_routed_delivery',
                                         'semi_routed_return')
                        }

    def __init__(self, fmp_name, fmp_ws="", nrow=0, ncol=0, nfarm=0, ncrop=0,
                 nirrigate=0, load=True, cell_area=None, crop_lut={}):
        super(Fmp, self).__init__()
        self.name = fmp_name
        self.ws = fmp_ws
        self.nrow = nrow
        self.ncol = ncol
        self.nfarm = nfarm
        self.ncrop = ncrop
        self.nirrigate = nirrigate
        self.crop_lut = crop_lut

        self.transient_arrays = {}
        self.boundary_arrays = {}
        self.crop_arrays = {}
        self.crop_areas = {}

        if isinstance(cell_area, np.ndarray):
            if cell_area.shape == (self.nrow, self.ncol):
                self.cell_area = cell_area
            else:
               raise AssertionError("cell area array must be the same "
                                    "shape as nrow, ncol or a scalar")
            pass
        else:
            if cell_area is None:
                self.cell_area = None
            else:
                self.cell_area = np.ones((self.nrow, self.ncol)) * cell_area

        if load:
            self.__block_indexes = {}
            self.__input = []
            self.__read_fmp_file()

    def __prepare_line(self, line):
        """
        Strips all characters right of a #

        Args:
            line: (str) line from fmp file

        Returns:
            t: (str) line that has been stripped
        """
        t = ""
        for i in line.strip():
            if i == '#':
                return t
            else:
                t += i
        return t

    def __find_input_blocks(self):
        """
        Finds the line index of each of the parameter blocks and sets
        them to the param block dictionary
        """
        block = False
        fmp_file_list = []
        with open(os.path.join(self.ws, self.name)) as fmp:
            for line in fmp:
                if line.strip().startswith('#') or \
                                line.strip() == "":
                    pass
                else:
                    line = self.__prepare_line(line)
                    fmp_file_list.append(line)

        for idx, line in enumerate(fmp_file_list):
            if block:
                if 'end' in line.lower():
                    self.__block_indexes[block].append(idx)
                    block = False
            else:
                if line.lower().split()[-1] in Fmp.block_names:
                    block = line.lower().split()[-1]
                    self.__block_indexes[block] = [idx + 1]

                elif " ".join(line.lower().split()[-2:]) in Fmp.block_names:
                    block = " ".join(line.lower().split()[-2:])
                    self.__block_indexes[block] = [idx + 1]

        self.__input = fmp_file_list

    def __read_fmp_file(self):
        """
        Meat and potatoes fmp reader
        """
        self.__find_input_blocks()

        for key, value in self.__block_indexes.items():
            self[key] = {}
            block_data = self.__input[value[0]:value[1]]
            linefeed = False
            skip = -1

            for ix, record in enumerate(block_data):
                t = record.split()

                if linefeed:
                    data = self.__linefeed_reader(t[-1])
                    linefeed = False
                    self[key][p_name] = data

                elif t[0].lower() not in Fmp.valid_parameters[key] and ix > skip:
                    raise AssertionError('Keyword not in valid_parameters')

                elif ix <= skip:
                    pass

                else:
                    data = None
                    p_name = t[0].lower()
                    if p_name == 'print':
                        p_name = t[0].lower() + "\t" + t[1].lower()
                        data = "\t".join(t[2:])

                    elif len(t) == 1:
                        data = True
                        dtype = 'bool'

                    elif t[1].lower() == "transient" and\
                            'internal' not in record.lower():

                        if t[2].lower() == 'list':

                            if t[3].lower() == 'datafile':
                                data = self.__transient_datafile_reader(t[-1])
                                dtype = 'transient datafile'

                            else:
                                data = self.__read_tfr_file(t[-1], list_fm=True)
                                dtype = 'transient list'

                        else:
                            data = self.__read_tfr_file(t[-1])
                            dtype = 'transient array'

                    elif t[1].lower() == 'static' and\
                            'internal' not in record.lower():

                        if t[2].lower() == 'list':

                            if t[3].lower() == 'datafile':
                                data = self.__datafile_reader(t[-1])
                                dtype = "static datafile"

                            else:
                                data = self.__list_reader(t[-1])
                                dtype = "static list"

                        else:
                            # todo: implement a static array instance maybe?
                            pass

                    elif t[1].lower() == "constant":
                        data = float(t[-1])
                        dtype = "constant"

                    elif "internal" in record.lower():
                        data, skip = self.__internal_reader(ix, block_data, key)
                        dtype = "{} {}".format(t[1].lower(), t[2].lower())

                    elif t[1].lower() == "open/close":
                        data = self.__array_reader(t[-1])
                        dtype = 'static array'

                    elif t[0].lower() == "linefeed":
                        linefeed = True
                        p_name = t[0] + ' ' + t[1]

                    else:
                        try:
                            data = int(t[-1])
                            dtype = "integer"
                        except ValueError:
                            data = t[-1]
                            dtype = "string"

                    if not linefeed:
                        self[key][p_name] = {'data': data, 'dtype': dtype}

    def __read_tfr_file(self, tfr_file, list_fm=False):
        """
        Parses the TFR file
        Args
            tfr_file: (str) tfr file name
            list_fm: (bool) True = list format;
                            False = array format
        Returns
            tfr_data: (np.ndarray)
        """
        tfr_data = []
        with open(os.path.join(self.ws, tfr_file)) as tfr:
            for line in tfr:
                if line.strip().startswith('#') or \
                                line.strip() == "":
                    pass
                else:
                    line = self.__prepare_line(line)
                    t = line.split()
                    sf = 1.
                    for idx, entry in enumerate(t):
                        if entry in ('SF', 'sf', 'Sf', 'sF'):
                            sf *= float(t[idx + 1])

                    if t[0].lower() == "open/close":
                        if list_fm:
                            tfr_data.append(self.__array_reader(t[1]))
                        else:
                            tfr_data.append(self.__array_reader(t[1]) * sf)

                    elif t[0].lower() == 'repeat':
                        tfr_data.append(tfr_data[-1])

                    else:
                        pass

        return np.array(tfr_data)

    def __array_reader(self, array, internal=False):
        """
        Read UR2dARRAY instances
        think about making this a custom loader!!!!
        """
        if internal:
            # todo: create an internal array data reader
            pass
        else:
            return np.genfromtxt(os.path.join(self.ws, array))

    def __transient_array_reader(self):
        """
        Read transient array records from file
        """
        pass

    def __datafile_reader(self, datafile):
        """
        Read static datafiles

        Args:
            datafile: (str) datafile name

        Returns:
            data: (list)
        """
        data = []
        with open(os.path.join(self.ws, datafile)) as df:
            sf = 1.
            for line in df:
                line = self.__prepare_line(line)
                t = line.split()
                if line:
                    if 'sfac' == t[0].lower():
                        sf = float(t[-1])

                    else:
                        data.append([int(t[0]), float(t[1]) * sf])

        return data

    def __transient_datafile_reader(self, datafile):
        """
        Read transient datafiles

        Args:
            datafile: (str) datafile name

        Returns:
            data: (list)
        """
        data = []
        with open(os.path.join(self.ws, datafile)) as df:
            sf = 1.
            temp = []
            for line in df:
                line = self.__prepare_line(line)
                t = line.split()
                if line:
                    if 'sfac' == t[0].lower():
                        sf = float(t[-1])
                        if temp:
                            data.append(temp)
                            temp = []
                    else:
                        temp.append([int(t[0]), float(t[1]) * sf])

        return data

    def __internal_reader(self, idx, block, block_name):
        """
        Performs internal read on fmp file

        Args:
            idx: (int) __input index
            block (list) block of parameters

        Returns:
            data: ?
        """
        fmp_data = block[idx:]
        t = []
        skip = 0
        for ix, line in enumerate(fmp_data):
            if ix == 0:
                heading = line.split()

            elif line.split()[0].lower() in Fmp.valid_parameters[block_name]:
                skip = ix + idx - 1
                break

            else:
                t.append(line)
                skip = ix + idx

        if heading[1].lower() == 'static':
            if heading[2].lower() == 'list':
                arr = self.__list_reader(t, internal=True)

                return arr, skip
        else:
            pass

    def __list_reader(self, lst, internal=False):
        """
        Read list records from file
        """
        if internal:
            arr = []
            for line in lst:
                arr.append(line.split())

            arr = np.array(arr, dtype=np.float)
        else:
            with open(os.path.join(self.ws, lst)) as lstf:
                sf = 1.
                arr = []
                for line in lstf:
                    line = self.__prepare_line(line)
                    t = line.split()
                    if line:
                        if 'sfac' == t[0].lower():
                            sf = float(t[-1])

                        else:
                            t0 = [int(t[0])]
                            t1 = [float(i) * sf for i in t[1:]]
                            arr.append(t0 + t1)

        return arr

    def __linefeed_reader(self, linefeed):
        """
        Read linefeed from file, ie. fmp wells input
        """
        temp = FmpLinefeed(linefeed, ws=self.ws)
        return temp.combined_array

    def to_arrays(self):
        """
        method to create row x col numpy arrays for all data within the
        farm processes.
        """
        self.nfarm = self['global dimension']['nfarm']
        self.ncrop = self['global dimension']['ncrop']
        self.nsoil = self['global dimension']['nsoil']
        self.nirrigate = self['global dimension']['nirrigate']

        self.boundary_arrays['soil_id'] = self['soil']['soil_id']['data']

        t = self['wbs']['location']['data']
        t.shape = (t.shape[0], 1, t.shape[1], t.shape[2])
        self.transient_arrays['location'] = t  # self['wbs']['location']['data']

        for key, item in self.items():
            if key in ('soil', 'climate', 'crop',
                       'fwell', 'wbs', 'surfacewater'):

                for param, dset in item.items():
                    if param.lower() == "linefeed wbs":
                        pass
                    else:
                        if dset['dtype'] == 'transient array':
                            if param == 'fraction':
                                s1 = dset['data'].shape[0]
                                dset['data'].shape = (s1, -1,
                                                      self.nrow,
                                                      self.ncol)
                                if self.crop_lut:
                                    kper = dset['data'].shape[0]
                                    for key, crop in self.crop_lut.items():
                                        arr = np.zeros((kper,
                                                        1,
                                                        self.nrow,
                                                        self.ncol))
                                        for per, krec in enumerate(dset['data']):
                                            arr[per] = np.array([krec[key - 1]])

                                        if self.cell_area is not None:
                                            area_arr = arr * np.array([self.cell_area])
                                            self.crop_areas[crop] = area_arr

                                        self.crop_arrays[crop] = arr


                                else:
                                    self.transient_arrays[param] = dset['data']

                            elif param == 'location':
                                pass

                            else:
                                t = dset['data']
                                t.shape = (t.shape[0], 1,
                                           t.shape[1], t.shape[2])
                                self.transient_arrays[param] = dset['data']

                        elif dset['dtype'] == 'static array':
                            # todo: maybe check the key and if its pos. transient. repr. as trans.
                            self.boundary_arrays[param] = dset['data']

                        elif dset['dtype'] == 'transient list':
                            if param in ('non_routed_delievery',
                                         'semi_routed_return'):
                                pass

                        elif dset['dtype'] == 'static list':
                            #todo: redo chunks of this to link crop to location!
                            if param in ("crop_coeficient", "irrigation",
                                         "root_depth", "pond_depth",
                                         "surfacewater_loss_fraction_precipitation",
                                         "surfacewater_loss_fraction_irrigation",
                                         "groundwater_root_interaction"):
                                pass
                                """
                                array_set = []
                                for arr in self.transient_arrays['location']:
                                    arr = np.copy(arr)
                                    for line in dset['data']:
                                        arr[arr == line[0]] = line[1]
                                    array_set.append(arr)

                                self.transient_arrays[key] = np.array(array_set)
                                """

                            elif param in ("capillary_fringe",):
                                # todo: link to soil array
                                pass

                        elif dset['dtype'] == 'transient datafile':
                            if param in ("crop_coeficient", ):
                                # todo: link to the crop array
                                pass


                        elif dset['dtype'] == 'static datafile':
                            if key in ('transpiration_fraction',
                                       'evaporation_irrigation_fraction'):
                                # todo: figure out how to link this one?
                                # todo: these are by farm--by crop
                                pass


if __name__ == "__main__":


    import flopy as fp

    crop_lut = {1: 'alfalfa', 2:'apples', 3: "grain",
                4: 'grapes', 5: 'grass', 6: 'jujube',
                7: 'lake', 8: 'landscaping', 9: 'mixed',
                10: 'orchard', 11: 'pasture', 12: 'pistachios',
                13: 'row_crop', 14: 'desert'}

    ws = r"C:\Users\jlarsen\Desktop\Lucerne\Lucerne_OWHM\V2_combined_OWHM"
    fmp_name = os.path.join(ws, "Lucerne.fmp")
    dis_name = os.path.join(ws, "Lucerne.dis")
    # shp_name = r"C:\Users\jlarsen\Desktop\Lucerne\GIS\V2_NAD83\Sfr_v2.shp"

    ml = fp.modflow.Modflow('Lucerne', model_ws=ws)

    dis = fp.modflow.ModflowDis.load(dis_name, ml,
                                     ext_unit_dict={},
                                     check=False)

    ml.dis.sr.xul = 493822.36 * 3.28084
    ml.dis.sr.yul = 3837490.09 * 3.28084

    xgrid = dis.sr.xgrid / 3.28084
    ygrid = dis.sr.ygrid / 3.28084
    area = 2000.*2000.

    fmp = Fmp(fmp_name, fmp_ws=ws, nrow=dis.nrow,
              ncol=dis.ncol, cell_area=area, crop_lut=crop_lut)
    fmp.to_arrays()
    print('break')

    ws = r'C:\Users\jlarsen\Desktop\Lucerne\GIS\V2_NAD83'
    shp_name = "Crop_arrays.shp"
    ATS.create_shapefile_from_transient_array(os.path.join(ws, shp_name),
                                              fmp.crop_arrays, dis.nper, 1,
                                              xgrid, ygrid, no_data=0.)

    shp_name = "Fmp_locations.shp"
    ATS.create_shapefile_from_transient_array(os.path.join(ws, shp_name),
                                              fmp.transient_arrays, dis.nper, 1,
                                              xgrid, ygrid, no_data=0.)
