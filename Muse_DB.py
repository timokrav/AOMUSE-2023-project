import logging
from pathlib import Path
import time
from datetime import datetime
import os
from astropy.io import fits
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from tqdm import tqdm
import numpy as np
import json
from pony.orm import *
import pymysql
import warnings

provider = "mysql"
host = "scdevdb1.sc.eso.org"
user = "aomuse"
passwd = "a0mu53!"
database_name = "aomuse"

"""
provider = "mysql"
host = "127.0.0.1"
user = "aomuse"
passwd = "#aomuse2020"
database_name = "timotest"
"""

##############################
# Actual script operations start from here
##############################

db = Database()


# Target class describes an entry to the Target table in the database
# The classes have to inherit db.Entity from Pony
class Target(db.Entity):
    #   ----- Attributes -----

    target_name = Required(str, unique=True)  # Required: Cannot be None

    #   ----- Relations -----

    exposures = Set('Exposure')  # One target contains a set of exposures
    processed_exposure = Optional('Processed_Exposure')

# Exposure table class
class Exposure(db.Entity):
    #   ----- Attributes -----

    observation_time = Required(datetime, unique=True)
    obs_id = Required(int, size=32, unsigned=True)
    insMode = Required(str)
    datacube_header = Optional(Json)
    raw_exposure_header = Optional(Json)
    raw_exposure_data = Optional(Json)
    raw_exposure_filename = Optional(str)
    prm_filename = Optional(str)
    pampelmuse_params = Optional(Json)
    sources = Optional(Json)
    pampelmuse_catalog = Optional(Json)
    raman_image_header = Optional(Json)
    maoppy_data = Optional(Json)

    #   ----- Sky parameters -----
    sky_condition_start_time = Optional(float)
    sky_condition_start = Optional(LongStr)
    sky_comment_start = Optional(LongStr)
    sky_condition_end_time = Optional(float)
    sky_condition_end = Optional(LongStr)
    sky_comment_end = Optional(LongStr)

    #   ----- Relations -----

    target = Required('Target')  # One exposure belongs to a target
    processed_exposure = Optional('Processed_Exposure')

class Processed_Exposure(db.Entity):
    observation_time = Required(datetime, unique=True)
    obs_id = Required(int, size=32, unsigned=True)
    insMode = Required(str)
    raw_filename = Optional(str)
    ngs_flux = Optional(float)
    ttfree = Optional(bool)
    degraded = Optional(bool)
    glf = Optional(float)
    seeing = Optional(float)
    seeing_los = Optional(float)
    airMass = Optional(float)
    tau0 = Optional(float)
    # --------------------------------------------------------------
    num_sources = Optional(int, unsigned=True)
    sgs_data = Optional(Json) # sgs_data extension
    ag_data = Optional(Json) # ag_data extension
    sparta_cn2 = Optional(Json) # sparta_cn2 extension
    sparta_atm =Optional(Json) # sparta_atm extension
    psf_params = Optional(Json)
    sparta_iq_data = Optional(Json)
    sparta_iq_los_500nm = Optional(float)
    sparta_iq_los_500nm_nogain = Optional(float)
    
    # Relations
    
    target = Required('Target')  # One exposure belongs to a target
    exposure = Required('Exposure')


# This simplifies some header value fetches
def fetch_data(header, keyword):
    try:
        return header[keyword]
    except KeyError:
        print('Keyword ' + keyword + ' not found.')


# This method converts an np.int64 to python native int, because the json library
# at the moment cannot deal with numpy int64
# N.B. I think numpy numbers are now deprecated, so this might not do anything anymore. - Timo 18/01/23
def convert_npint64_to_int(o):
    if isinstance(o, np.int64):
        return int(o)
    raise TypeError


@db_session  # Decorator from Pony to make the function into a db session, thus the database can be modified
def muse_script():
    ##############################
    # This function read 5 files (prm, psf, reduced, raman and raw) to store the information of each exposure in a database.
    #
    # For the prm file reads the PSFPARS extension and stores it in the 'psfParams' field of the exposure.
    # For the psf file reads the 'PM_X' extensions, where X is the id of the source, and stores it in the
    # 'sources' field of the exposure.
    # With other files the scripts takes the headers and saves them as json-like LongStr
    # Also, from the primary header it stores the target and instrument mode of the exposure in their own fields
    # to make it easier the analysis from differents targets and instrument modes.
    ##############################

    new_entries = 0
    modified_entries = 0
    not_modified = 0
    rootDir = input("Enter the root directory: ")
    singleDir = rootDir + '/single/' # Reduced files
    rawDir = rootDir + '/raw/' # Raw files
    analysisDir = rootDir + '/analysis/' # Prm and psf files
    catalogDir = rootDir + '/cats/' # SExtractor files

    # Store the directories to find it later by the namefile quickly, since each file is inside a date format folder

    single_files = {}
    raw_files = {}
    catalog_files = {}
    nightlog_files = {}
    raman_files = {}

    single_count = 0
    raw_count = 0
    catalog_count = 0
    nightlog_count = 0
    raman_count = 0
    prm_count = 0

    start = time.time()
    logging.info("Listing single files directories...")
    print("Listing single files directories...")
    time.sleep(0.1)  # Makes the terminal output pretty
    for filepath in tqdm(list(Path(singleDir).rglob('*.fits'))):
        single_files[filepath.name] = filepath
        single_count += 1
    logging.info(f"{single_count} single files found.")
    print(f"{single_count} single files found.\n")

    logging.info("Listing raw files and nightlogs directories...")
    print("Listing raw files and nightlogs directories...")
    time.sleep(0.1)  # Makes the terminal output pretty
    for filepath in tqdm(list(Path(rawDir).rglob('*.fits.fz'))):
        raw_files[filepath.name] = filepath
        raw_count += 1
        nightlog_path = Path(str(filepath)[:-7] + "NL.txt")
        if Path.exists(nightlog_path):
            nightlog_files[filepath.name] = nightlog_path
            nightlog_count += 1
    logging.info(f"{raw_count} raw files found.")
    print(f"{raw_count} raw files found.")
    logging.info(f"{nightlog_count} nightlog files found.")
    print(f"{nightlog_count} nightlog files found.\n")

    logging.info("Listing catalog files directories...")
    print("Listing catalog files directories...")
    time.sleep(0.1)  # Makes the terminal output pretty
    for filepath in tqdm(list(Path(rootDir).rglob("*.detections.cat"))):
        catalog_files[filepath.name.split("_")[0] + '.fits'] = filepath
        catalog_count += 1
    logging.info(f"{catalog_count} catalog files found.")
    print(f"{catalog_count} catalog files found.\n")

    logging.info("Listing raman files directories...")
    print("Listing raman files directories...")
    time.sleep(0.1)  # Makes the terminal output pretty
    for filepath in tqdm(list(Path(rootDir).rglob("*.fits"))):
        break # REMOVE ME
        try:
            header = fits.getheader(filepath)
            if 'HIERARCH ESO PRO CATG' in header and header['HIERARCH ESO PRO CATG'] == "RAMAN_IMAGES":
                raman_files[header['DATE-OBS']] = filepath
                raman_count += 1
        except:
            pass
    logging.info(f"{raman_count} raman files found.")
    print(f"{raman_count} raman files found.\n")

    folders = os.listdir(analysisDir)
    analysisDir_moffat = ""
    analysisDir_maoppy = ""
    for folder in folders:
        if("moffat" in folder):
            analysisDir_moffat = analysisDir + folder + "/"
        elif("maoppy" in folder):
            analysisDir_maoppy = analysisDir + folder + "/"
    NFM = True
    if(analysisDir_moffat == "" or analysisDir_maoppy == ""):
        NFM = False
    prm_files = {}
    if(NFM):
        logging.info("NFM exposures found")
        print("NFM exposures found")
        logging.info("Listing prm files directories...")
        print("Listing prm files directories...")
        time.sleep(0.1)  # Makes the terminal output pretty
        for filepath in tqdm(list(Path(analysisDir_maoppy).rglob('*.prm.fits'))):
            prm_files[filepath.name] = filepath
            prm_count += 1
        logging.info(f"{single_count} prm files found.")
        print(f"{single_count} prm files found.\n")

    end = time.time()
    print('{:.3f}'.format(end-start) + " seconds to finish file search\n")
    # ----- Analysis Folder -----

    # The script starts with the analysis files (psf and prm files)    
    # FITS file scan and identification
    start = time.time()
    logging.info("Iterating analysis folder for prm files...")
    print("Creating/updating and storing exposure objects into the database")
    
        
    if(NFM):
        prm_fits = list(Path(analysisDir_moffat).rglob('*prm*.fits'))
    else:
        prm_fits = list(Path(analysisDir).rglob('*prm*.fits'))
    warnings = 0
    for prm_file in tqdm(prm_fits):

        logging.info(f"Manipulating file {prm_file}")

        observation_dictionary = {}
        prm_filename = prm_file.name
        expected_cube_name = prm_filename.replace('.prm', '')

        # ----- Single Folder -----

        # Reduced File

        try:
            cube_parameters = {}
            logging.info(f"Attempting to guess the cube name {single_files[expected_cube_name]}")
            header = fits.getheader(single_files[expected_cube_name])
            logging.info(f"Found a corresponding cube and header.")
            
            try:
                target = header['OBJECT']  # Obtain the target name
            except KeyError:
                logging.warning(str(expected_cube_name) + " does not have a header key OBJECT, skipping file.")
                warnings += 1
                print("The header OBJECT does not exist")
                continue
            logging.info(f"Testing target {target} existence...")
            target_exists = db.exists("select id from target where target_name=$target")
            logging.info(f"Test complete.")
            if target_exists:
                logging.info(f"Fetching the target...")
                cube_parameters['target'] = Target.get_by_sql(f"select * from target where target_name = '{target}'")
            else:
                logging.info(f"Generating a new target...")
                cube_parameters['target'] = Target(target_name=target)  # Else, create the target
            logging.info(f"Successful test.")
            
            try:
                cube_parameters['observation_time'] = header['DATE-OBS']
            except KeyError:
                logging.warning(str(expected_cube_name) + " the header DATE-OBS does not exist")
                warnings += 1
                print("The header DATE-OBS does not exist")
                continue
            except Exception as e:
                logging.warning(str(expected_cube_name) + " " + str(e))
                warnings += 1
                continue

            try:
                cube_parameters["obs_id"] = header["HIERARCH ESO OBS ID"]
            except KeyError:
                logging.warning(str(expected_cube_name) + " the header HIERARCH ESO OBS ID does not exist")
                warnings += 1
                print("The header HIERARCH ESO OBS ID does not exist")
                continue
            except Exception as e:
                logging.warning(str(expected_cube_name) + " " + str(e))
                warnings += 1
                continue

            try:
                raw_filename = header['PROV1']+'.fz'
            except KeyError:
                logging.warning(f"{expected_cube_name} header PROV1 does not exist. Can't find raw filename")
                warnings += 1
                print(f"{expected_cube_name} header PROV1 does not exist. Can't find raw filename")

            cube_parameters['instrument_mode'] = fetch_data(header, 'HIERARCH ESO INS MODE') # CHECK: what happens with insmode key
            cube_parameters['header'] = dict(header)
            try:
                del cube_parameters['header']['COMMENT']  # and then delete the COMMENT key, because it does not have the JSON format and
            except KeyError:  # throws an error when saving it in the database
                pass
        except Exception as e:
            logging.warning(f"Error trying to get header of reduced file {expected_cube_name}. Error: {e}")
            warnings += 1
            continue

        observation_dictionary['target'] = cube_parameters['target'] 
        observation_dictionary['observation_time'] = cube_parameters['observation_time'][:-4]
        observation_dictionary["obs_id"] = cube_parameters["obs_id"]
        observation_dictionary['insMode'] = cube_parameters['instrument_mode']
        observation_dictionary['datacube_header'] = dict(header)
        try:
            del observation_dictionary['datacube_header']['COMMENT']  # and then delete the COMMENT key, because it does not have the JSON format and
        except:  # throws an error when saving it in the database
            pass
        try:
            del observation_dictionary['datacube_header']['']  # and then delete the '' key, because it does not have the JSON format and
        except:  # throws an error when saving it in the database
            pass

        observation_dictionary['raw_exposure_filename'] = raw_filename

		# These flags decide what file extractions will be done
        new_expo = False
        read_raw = False
        read_prm = False
        read_psf = False
        read_catalog = False
        read_nightlog = False
        read_raman = False

        # Check if exposure exists, create a new exposure if not and modify existing one if needed
        logging.info(f"Finished with the header, starting with the exposure.")
        
        logging.info(f"Testing existence of {observation_dictionary['observation_time']}")
        logging.info(prm_filename)
        exposure_exists = db.exists("select id from exposure where observation_time = $observation_dictionary['observation_time']")
        
        if not exposure_exists:
            logging.info("Didn't find the exposure in database")
            new_expo = True
        else:
            logging.info("Did find the exposure in database, retrieving...")
            database_entry = Exposure.get_by_sql("select * from exposure where observation_time = $observation_dictionary['observation_time']")
            raw_attr = ["raw_exposure_header", "raw_exposure_data", "raw_exposure_filename"]
            prm_attr = ["prm_filename", "pampelmuse_params"]
            nightlog_attr = [
                "sky_condition_start_time",
                "sky_condition_start",
                "sky_comment_start",
                "sky_condition_end_time",
                "sky_condition_end",
                "sky_comment_end"]

            # Cycles through columns associated with different files. If they have empty entries, they should be reread.
            for attr in raw_attr:
                database_value = getattr(database_entry, attr)
                if(database_value == 'null' or database_value == None or database_value == ""):
                    read_raw = True
                    break
            for attr in prm_attr:
                database_value = getattr(database_entry, attr)
                if(database_value == 'null' or database_value == None or database_value == ""):
                    read_prm = True
                    break
            database_value = getattr(database_entry, "sources")
            if(database_value == 'null' or database_value == None or database_value == ""):
                read_psf = True
            database_value = getattr(database_entry, "pampelmuse_catalog")
            if(database_value == 'null' or database_value == None or database_value == ""):
                read_catalog = True
            database_value = getattr(database_entry, "raman_image_header")
            if(database_value == 'null' or database_value == None or database_value == ""):
                read_raman = True
            for attr in nightlog_attr:
                database_value = getattr(database_entry, attr)
                if(database_value == 'null' or database_value == None or database_value == ""):
                    read_nightlog = True
                    break

        # Nothing needs to be modified so we can move on to the next file
        if(not new_expo and not read_raw and not read_prm and not read_catalog and not read_nightlog and not read_raman):
            not_modified += 1
            continue
            
        # PRM File
        # Notice that single file has the same root name that prm file
        if(new_expo or read_prm):
            psfParams = {}  # Dictionary to store the extensions of the prm files
            try:
                with fits.open(prm_file) as hduList:

                    try:
                        params = hduList['PSFPARS'].data  # If the prm file does not have the PSFPARS extension, skip
                    except KeyError:
                        logging.warning(prm_filename + " Could not find PSFPARS in header. Possibly a broken prm file, skipping file.")
                        warnings += 1
                        psfParams = None
                    if(psfParams != None):
                        try:
                            for i in range(len(params['name'])):  # The column 'name' is a list with the parameter names
                                parameter = params['name'][i]
                                value = params['polyfit'][i].tolist()  # The column 'polyfit' is a list of lists
                                psfParams[parameter] = value  # Store the polyfit of the PSF parameters
                            prmStep = hduList['SPECTRA'].header['CDELT3']
                            prmRestw = hduList['SPECTRA'].header['CRVAL3']
                            prmData = np.arange(hduList['SPECTRA'].header['NAXIS3'])  # Calculate the wavelength
                            prmWavelength = (prmRestw + (prmData * prmStep)) * 10 ** 9
                            psfParams['wavelength'] = prmWavelength.tolist()
                        except Exception as e:
                            print(str(prm_file) + " " + str(e))
                            logging.warning(str(prm_file) + " " + str(e))
                            warnings += 1
                            psfParams = None
            except FileNotFoundError:
                logging.warning("The file " + str(prm_filename) + " does not exist.")
                warnings += 1
                print(f"The file {prm_filename} does not exist\n")  # If the prm file does not exists, skip
                psfParams = None
                prm_filename = None
            except Exception as e:
                logging.warning(str(prm_filename) + " " + str(e))
                warnings += 1
                psfParams = None
                prm_filename = None

            observation_dictionary['pampelmuse_params'] = psfParams
            observation_dictionary['prm_filename'] = prm_filename

        # PSF File
        if(new_expo or read_psf):
            sources = {}  # Dictionary to store the extensions (each source) of the psf files
            expected_psf_file = str(prm_file).replace('.prm', '.psf')
            try:
                with fits.open(expected_psf_file) as hduList:
                    for tupla in hduList.info(False):  # Iterate in the list of the extensions
                        if ('PM_' in tupla[1]):  # Looking for the 'PM_X' extensions, where X is the id of the source
                            sourceID = tupla[1]
                            sourceData = {}  # Dictionary that will store the data of the source
                            table = hduList[sourceID].data  # Access to the source extension
                            i = 0
                            for column in table.columns.names:  # Iterate in the list of column names
                                sourceData[column] = table[column].tolist()  # Store the column with its list of values
                                i += 1
                            sources[sourceID] = sourceData  # Store the source data
                    if sources == {}:
                        print("Possibly a bad psf file, no sources.")
                        logging.warning(expected_psf_file.name + " Possibly a bad psf file, PM_* layers are missing")
                        warnings += 1
            except FileNotFoundError:
                logging.warning("The file " + str(expected_psf_file.name) + " does not exist.")
                warnings += 1
                print(f"The file {expected_psf_file} does not exist\n")  # If the psf file does not exists, skip
                sources = None
            except Exception as e:
                logging.warning(f"{expected_psf_file} {str(e)}")
                warnings += 1
                sources = None

            observation_dictionary['sources'] = sources

        if(NFM):
            # Maoppy PRM File
            psfParams = {}  # Dictionary to store the extensions of the prm files
            observation_dictionary["maoppy_data"] = {}
            try:
                maoppy_filepath = prm_files[prm_file.name]
                with fits.open(maoppy_filepath) as hduList:

                    try:
                        params = hduList['PSFPARS'].data  # If the prm file does not have the PSFPARS extension, skip
                    except KeyError:
                        logging.warning("Maoppy " + prm_filename + " Could not find PSFPARS in header. Possibly a broken prm file, skipping file.")
                        warnings += 1
                        psfParams = None
                    if(psfParams != None):
                        try:
                            for i in range(len(params['name'])):  # The column 'name' is a list with the parameter names
                                parameter = params['name'][i]
                                value = params['polyfit'][i].tolist()  # The column 'polyfit' is a list of lists
                                psfParams[parameter] = value  # Store the polyfit of the PSF parameters
                            prmStep = hduList['SPECTRA'].header['CDELT3']
                            prmRestw = hduList['SPECTRA'].header['CRVAL3']
                            prmData = np.arange(hduList['SPECTRA'].header['NAXIS3'])  # Calculate the wavelength
                            prmWavelength = (prmRestw + (prmData * prmStep)) * 10 ** 9
                            psfParams['wavelength'] = prmWavelength.tolist()
                        except Exception as e:
                            print("Maoppy " + str(prm_file) + " " + str(e))
                            logging.warning("Maoppy " + str(prm_file) + " " + str(e))
                            warnings += 1
                            psfParams = None
            except FileNotFoundError:
                logging.warning("The file " + "Maoppy " + str(prm_filename) + " does not exist.")
                warnings += 1
                print(f"The file Maoppy {prm_filename} does not exist\n")  # If the prm file does not exists, skip
                psfParams = None
                prm_filename = None
            except Exception as e:
                logging.warning("Maoppy " + str(prm_filename) + " " + str(e))
                warnings += 1
                psfParams = None
                prm_filename = None

            observation_dictionary["maoppy_data"]['pampelmuse_params'] = psfParams
            observation_dictionary["maoppy_data"]['prm_filename'] = prm_filename

            # Maoppy PSF File
        
            sources = {}  # Dictionary to store the extensions (each source) of the psf files
            expected_psf_file = str(maoppy_filepath).replace('.prm', '.psf')
            try:
                with fits.open(expected_psf_file) as hduList:
                    for tupla in hduList.info(False):  # Iterate in the list of the extensions
                        if ('PM_' in tupla[1]):  # Looking for the 'PM_X' extensions, where X is the id of the source
                            sourceID = tupla[1]
                            sourceData = {}  # Dictionary that will store the data of the source
                            table = hduList[sourceID].data  # Access to the source extension
                            i = 0
                            for column in table.columns.names:  # Iterate in the list of column names
                                sourceData[column] = table[column].tolist()  # Store the column with its list of values
                                i += 1
                            sources[sourceID] = sourceData  # Store the source data
                    if sources == {}:
                        print("Possibly a bad Maoppy psf file, no sources.")
                        logging.warning("Maoppy " + expected_psf_file.name + " Possibly a bad psf file, PM_* layers are missing")
                        warnings += 1
            except FileNotFoundError:
                logging.warning("The file " + "Maoppy " + str(expected_psf_file.name) + " does not exist.")
                warnings += 1
                print(f"The file Maoppy {expected_psf_file.name} does not exist\n")  # If the psf file does not exists, skip
                sources = None
            except Exception as e:
                logging.warning("Maoppy " + str(expected_psf_file.name) + " " + str(e))
                warnings += 1
                sources = None

            observation_dictionary["maoppy_data"]['sources'] = sources

    
        if(new_expo or read_raw):
            #######################################
            # Raw exposure extraction starts here
            #######################################

            raw_parameters = {}

            logging.info("Raw exposure start.")
            data = {}
                
            # SGS and SPARTA header extraction    
            try:
                #data = {}
                flag = False # Flag not relevant, just to print the \n before a missed extension
                with fits.open(raw_files[raw_filename]) as hduList:
                    observation_dictionary['raw_exposure_filename'] = raw_filename 
                    try:
                        header = hduList[0].header
                    except Exception as e:
                        logging.warning(f"Error getting header from {raw_filename}")
                        warnings += 1
                        print(f"Error getting header from {raw_filename}")
                    try:
                        for tupla in hduList.info(False): # Iterate in the list of the extensions
                            if('CHAN' in tupla[1]): # Looking for the CHAN extensions
                                data[tupla[1]] = dict(hduList[tupla[1]].header)
                    except:
                        logging.warning("Channel header not found " + raw_filename)   
                        warnings += 1           
                    try:
                        sgsData = {} # Dictionary that will store the data
                        table = hduList['SGS_DATA'].data 
                        for column in table.columns.names: # Iterate in the list of column names
                            sgsData[column] = table[column].tolist() # Store the column with its list of values
                        data['SGS_DATA'] = sgsData
                    except KeyError:
                        data['SGS_DATA'] = None # If the extension does not exists, stores it anyway as None
                        flag = True # Flag, not important, for the \n print
                        logging.warning("SGS_DATA not found in " + raw_filename)
                        warnings += 1
                    try:
                        agData = {}
                        table = hduList['AG_DATA'].data
                        for column in table.columns.names: # Iterate in the list of column names
                            agData[column] = table[column].tolist() # Store the column with its list of values
                        data['AG_DATA'] = agData    
                    except KeyError:
                        data['AG_DATA'] = None # If the extension does not exists, stores it anyway as None
                        flag = True # Flag, not important, for the \n print
                        logging.warning("AG_DATA not found in " + raw_filename)
                        warnings += 1
                    try:
                        asmData = {}
                        table = hduList['ASM_DATA'].data
                        for column in table.columns.names: # Iterate in the list of column names
                            asmData[column] = table[column].tolist() # Store the column with its list of values
                        data['ASM_DATA'] = asmData      
                    except KeyError:
                        data['ASM_DATA'] = None # If the extension does not exists, stores it anyway as None
                        flag = True # Flag, not important, for the \n print
                        logging.warning("ASM_DATA not found in " + raw_filename)  
                        warnings += 1           
                    try:
                        spartaAtmData = {}
                        table = hduList['SPARTA_ATM_DATA'].data
                        for column in table.columns.names: # Iterate in the list of column names
                            spartaAtmData[column] = table[column].tolist() # Store the column with its list of values
                        data['SPARTA_ATM_DATA'] = spartaAtmData
                    except KeyError:
                        data['SPARTA_ATM_DATA'] = None # If the extension does not exists, stores it anyway as None
                        flag = True # Flag, not important, for the \n print
                        logging.warning("SPARTA_ATM_DATA not found in " + raw_filename)
                        warnings += 1
                    try:
                        spartaCn2Data = {}
                        table = hduList['SPARTA_CN2_DATA'].data
                        for column in table.columns.names: # Iterate in the list of column names
                            spartaCn2Data[column] = table[column].tolist() # Store the column with its list of values
                        data['SPARTA_CN2_DATA'] = spartaCn2Data
                    except KeyError:
                        data['SPARTA_CN2_DATA'] = None # If the extension does not exists, stores it anyway as None
                        flag = True # Flag, not important, for the \n print
                        logging.warning("SPARTA_CN2_DATA not found in " + raw_filename)  
                        warnings += 1      
                    if(flag):
                        logging.warning("") # The \n print
            except FileNotFoundError:
                logging.warning("The file {raw_filename} does not exist\n") # If the raw file does not exists, skip
                warnings += 1
                header = None 
                data = {"error": "FileNotFound"}
                raw_filename = None
            except Exception as e:
                logging.warning(f"Error at reading raw file {raw_filename}")
                warnings += 1
                header = None 
                data = {"error": f"{e}"}
                raw_filename = None
            if(header != None):
                observation_dictionary['raw_exposure_header']  = dict(header)
            logging.info("Dumping variable 'data'...")
            logging.info(type(data))
            logging.info(data)
            observation_dictionary['raw_exposure_data']  = data
            try:
                del observation_dictionary['raw_exposure_header']['COMMENT']  # and then delete the COMMENT key, because it does not have the JSON format and
            except:  # throws an error when saving it in the database
                pass
            try:
                del observation_dictionary['raw_exposure_header']['']  # and then delete the '' key, because it does not have the JSON format and
            except:  # throws an error when saving it in the database
                pass


        ##############################
        # Catalog parsing starts here
        ##############################

        # detections.cat file scan and identification
        # pampelmuse source catalogs
        if(new_expo or read_catalog):
            catalog_filename = expected_cube_name
            try:
                with open(catalog_files[catalog_filename], "r") as cat_file:
                    catalog_text = cat_file.read()

                cat_dict = {'catalog': catalog_text.replace("\'", "").replace("\n", "\\n")}
                
                """
                catalog_table = Table.read(catalog_files[catalog_filename], format="ascii.ecsv")
                catalog_array = np.asarray(catalog_table)
                list_of_cat_entries = []
                for i in range(catalog_array.shape[0]):
                    entry_dict = {k: v for k, v in zip(catalog_array[i].dtype.names, catalog_array[i])}
                    list_of_cat_entries.append(entry_dict)
                """
            except Exception as e:
                logging.warning(f"{catalog_filename} cannot open {e}")
                warnings += 1
                list_of_cat_entries = None
            observation_dictionary['pampelmuse_catalog'] = cat_dict

        ##############################
        # Nightlog starts here
        ##############################
        if(new_expo or read_nightlog):
            skip = False
            try:
                with open(nightlog_files[raw_filename]) as f:
                    read_lines = f.readlines()
            except (FileExistsError, KeyError):
                logging.warning("Nightlog file for " + expected_cube_name + " not found.")
                warnings += 1
                start_weather_time = -1
                observation_start_condition = "None"
                observation_start_comment = "None"
                end_weather_time = -1
                observation_end_condition = "None"
                observation_end_comment = "None"
                skip = True

            if(not skip):
                try:
                    local_start_time = float(observation_dictionary['raw_exposure_header']["UTC"])
                    integration_time = float(observation_dictionary['raw_exposure_header']["EXPTIME"])
                except Exception as e:
                    logging.warning(raw_filename + " " + str(e))
                    warnings += 1
                    start_weather_time = -1
                    observation_start_condition = "None"
                    observation_start_comment = "None"
                    end_weather_time = -1
                    observation_end_condition = "None"
                    observation_end_comment = "None"
                    skip = True

                if(not skip):
                    local_end_time = (local_start_time + integration_time) % 86400.0
                    if local_start_time > 43200.0:
                        local_start_time -= 86400.0

                    line_extraction = []

                    # Read through the nightlog, starting from the weather report
                    # Take the hh:mm time and transform into UTC seconds
                    for line in read_lines[10:]:
                        if line == "---------------------------------------------------\n" or line == "\n":
                            break
                        line_split = line[:-1].split("\t")
                        try:
                            weather_time = float(line_split[0].split(":")[0]) * 3600.0 + float(line_split[0].split(":")[1]) * 60.0
                        except ValueError:
                            continue
                        except Exception as e:
                            logging.warning("Nightlog for " + expected_cube_name + " " + str(e))
                            warnings += 1
                            continue

                        if weather_time > 43200.0:
                            weather_time -= 86400.0
                        #print(weather_time)
                        weather_line = [weather_time,
                                        line_split[1],
                                        line_split[2]]
                        line_extraction.append(weather_line)

                    start_weather_time = 0.0
                    end_weather_time = 0.0
                    observation_start_condition = "None"
                    observation_end_condition = "None"
                    observation_start_comment = "None"
                    observation_end_comment = "None"
                    i = -1
                    # Loop through the weather comments
                    try:
                        for line in line_extraction:
                            weather_comment_time = line[0]
                            i += 1
                            if weather_comment_time > 43200:
                                weather_comment_time -= 86400
                            # Check we've reached past the exposure start
                            if weather_comment_time > local_start_time:
                                break
                            start_weather_time = line[0]
                            end_weather_time = line[0]
                            observation_start_condition = line[1]
                            observation_end_condition = line[1]
                            observation_start_comment = line[2]
                            observation_end_comment = line[2]

                        # Loop through the weather comments, but start where we left off previously
                        for line in line_extraction[i:]:
                            weather_comment_time = line[0]
                            if line[0] > 43200:
                                weather_comment_time -= 86400
                            # Check we've reached past the exposure end
                            if weather_comment_time > local_end_time:
                                break
                            end_weather_time = line[0]
                            observation_end_condition = line[1]
                            observation_end_comment = line[2]
                    except Exception as e:
                        logging.warning("Nightlog for " + expected_cube_name + " " + str(e))
                        warnings += 1
                        start_weather_time = -1
                        observation_start_condition = "None"
                        observation_start_comment = "None"
                        end_weather_time = -1
                        observation_end_condition = "None"
                        observation_end_comment = "None"

            observation_dictionary["sky_condition_start_time"] = start_weather_time
            observation_dictionary["sky_condition_start"] = observation_start_condition
            observation_dictionary["sky_comment_start"] = observation_start_comment
            observation_dictionary["sky_condition_end_time"] = end_weather_time
            observation_dictionary["sky_condition_end"] = observation_end_condition
            observation_dictionary["sky_comment_end"] = observation_end_comment


        # Raman file extraction
        if(new_expo or read_raman):
            raman_parameters = {}
            try:
                header = fits.getheader(raman_files[observation_dictionary['observation_time']])
                observation_dictionary['raman_image_header'] = dict(header)
                try:
                    del raman_parameters['header']['COMMENT']  # and then delete the COMMENT key, because it does not have the JSON format and
                except KeyError:  # throws an error when saving it in the database
                    pass
            except KeyError:
                observation_dictionary['raman_image_header'] = None
                pass
            except Exception as e:
                print(f"Error with raman {observation_dictionary['observation_time']} obs-date.\n{e}")
                logging.warning(f"Error with raman {observation_dictionary['observation_time']} obs-date.\n{e}")
                warnings += 1
                observation_dictionary['raman_image_header'] = None
            # Modify headers and tables into more JSON-like format
            if 'pampelmuse_catalog' in observation_dictionary:
                observation_dictionary['pampelmuse_catalog'] = json.dumps(observation_dictionary['pampelmuse_catalog'])#, default=convert_npint64_to_int)
            if 'pampelmuse_params' in observation_dictionary:
                observation_dictionary['pampelmuse_params'] = json.dumps(observation_dictionary['pampelmuse_params'])
            if 'sources' in observation_dictionary:
                observation_dictionary['sources'] = json.dumps(observation_dictionary['sources'])
            if 'raman_image_header' in observation_dictionary:
                observation_dictionary['raman_image_header'] = json.dumps(observation_dictionary['raman_image_header'])
            if 'raw_exposure_header' in observation_dictionary:
                observation_dictionary['raw_exposure_header'] = json.dumps(observation_dictionary['raw_exposure_header'])
            if 'raw_exposure_data' in observation_dictionary:
                observation_dictionary['raw_exposure_data'] = json.dumps(observation_dictionary['raw_exposure_data'])
            if 'datacube_header' in observation_dictionary:
                observation_dictionary['datacube_header'] = json.dumps(observation_dictionary['datacube_header'])

        if new_expo:
            logging.info("Adding a new exposure...")
            Exposure(**observation_dictionary)
            db.commit()
            new_entries += 1
            logging.info("Done!")
        else:
            try:
                entry_id = database_entry.id
                logging.info("Modifying an entry with observation time " + observation_dictionary['observation_time'])
                logging.info("observations dictionary:")
                logging.info(observation_dictionary.keys())
                
                # We will remove some fhe entries from the dictionary, they should not be changed in any event.
                if 'target' in observation_dictionary.keys():
                    del observation_dictionary['target']
                if 'observation_time' in observation_dictionary.keys():
                    del observation_dictionary['observation_time']
                if 'obs_id' in observation_dictionary.keys():
                    del observation_dictionary['obs_id']
                if 'insMode' in observation_dictionary.keys():
                    del observation_dictionary['insMode']
                if 'datacube_header' in observation_dictionary.keys():
                    del observation_dictionary['datacube_header']

                # Create a pymysql connection to the database and start changing the values
                connection = pymysql.connections.Connection(host=host, user=user, password=passwd, database=database_name)
                cursor = connection.cursor()
                
                entry_modified = False
                for key, value in observation_dictionary.items():
                    # Sanity check, so we don't accidentally overwrite old good with an empty entry. This might happen if there is
                    # e.g. a pampelmuse file with no other associated files.
                    if value == None:
                        continue
                    if value == "None":
                        continue
                    if value == "null":
                        continue
                    if value == {}:
                        continue
                    if "fits.fz" in value:
                        continue
                    if value == "prm.fits" in value:
                        continue
                    if value == -1:
                        continue
                    try:
                        logging.info(f"Editing column {key} with phrase:")
                        logging.info(f"UPDATE exposure SET {key} = '{value}' WHERE id = {entry_id}")
                        cursor.execute(f"UPDATE exposure SET {key} = '{value}' WHERE id = {entry_id}")
                        entry_modified = True
                    except Exception as e:
                        logging.info("Encountered an error during update.")
                        logging.info(e)
                
                logging.info("Trying to commit...")
                connection.commit()
                if entry_modified:
                    modified_entries += 1
                logging.info("Done!")
            except Exception as e:
                logging.info("Commit caused an error")
                logging.info(e)

    print(f"Finished updating the database, {new_entries} new entries, {modified_entries} modified entries and {not_modified} unmodified entries.")
    logging.info(f"Finished updating the database, {new_entries} new entries, {modified_entries} modified entries and {not_modified} unmodified entries.")
    end = time.time()
    print(f"{warnings} warnings registered.")
    if(warnings > 0):
        print("Check log file for more detail")
    print('{:.3f}'.format(end-start) + " seconds to finish")




# Main part starts here
log_filename = datetime.now().strftime('%d-%m-%Y_%H:%M:%S.log')
logging.basicConfig(filename=log_filename,
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%d-%m-%Y %H:%M:%S.log',
                    level=logging.INFO)
print("Started a logfile called " + log_filename)
try:
    sql_debug(True)
    warnings.simplefilter('ignore', category=AstropyWarning)
    db.bind(provider=provider, host=host, user=user, passwd=passwd, db=database_name)
    db.generate_mapping(check_tables=False, create_tables=True)
    logging.info("Connected to SQL database successfully.")
except Exception as e:
    print("Error with SQL connection, check the log for more information.")
    logging.error("SQL connection problem")
    logging.error(e)
    quit(2)

muse_script()
