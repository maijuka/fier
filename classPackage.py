import tempfile
import os.path
import plistlib
import shutil
import zipfile
import re
import io
import sys

import pdb

try:
    from toolLog import Log
    from classCol import ColM, DataError
    from classRedescription import Redescription
    from classData import Data
    from classQuery import Query
    from classConstraints import Constraints, ActionsRegistry
    from classPreferencesManager import getPreferencesReader, TmpDirManager, EXT_SIREN, PATT_QSUFF
except ModuleNotFoundError:
    from .toolLog import Log
    from .classCol import ColM
    from .classRedescription import Redescription
    from .classData import Data
    from .classQuery import Query
    from .classConstraints import Constraints, ActionsRegistry
    from .classPreferencesManager import getPreferencesReader, TmpDirManager, EXT_SIREN, PATT_QSUFF


def make_filepath(filename,  basic_parts=None, default_suff="_alt", default_ext=""):
    # filename can be either
    # - -> just return
    # + -> filepath is built from the basic_parts if they are not None
    # a full path -> directly returned as filepath
    # a basename (no path, but an extension) -> combined to the path from basic_parts is it is not None
    # a suffix (no path, no extension) -> combined to basic_parts, instead of default_suff
    if filename is None or len(filename) == 0:
        return None
    elif filename == "-":
        return filename
    ext = None
    if basic_parts is not None:
        ext = basic_parts[-1]
    if (default_ext is not None and len(default_ext) > 0):
        ext = default_ext

    if basic_parts is not None and filename == "+":  # filename is + -> make filename from basic
        return os.path.join(basic_parts[0], basic_parts[1] + default_suff + ext)
    elif len(filename) > 0 or (default_ext is not None and len(default_ext) > 0):
        head_n, tail_n = os.path.split(filename)
        if len(head_n) > 0:  # if folder is provided, use it
            return filename
        elif basic_parts is not None:  # otherwise use from basic if possible
            root_n, ext_n = os.path.splitext(tail_n)
            if len(ext_n) > 0:  # if extension -> complete basename
                return os.path.join(basic_parts[0], root_n + ext_n)
            return os.path.join(basic_parts[0], basic_parts[1] + filename + ext)

def make_sid(count_rlists=0):
    return ("%s" % (count_rlists + 1))[::-1]
        
class Package(object):
    """Class to handle the zip packages that contain data, preferences, results, etc. for redescription mining.
    """

    # CONSTANTS
    # Names of the files in the package
    DATA_FILENAMES = ["data_LHS.csv",
                      "data_RHS.csv"]

    REDESCRIPTIONS_FILENAME = "redescriptions.csv"
    PREFERENCES_FILENAME = "preferences.xml"
    FDEFS_FILENAME = {"fields_vdefs": "fields_vdefs_custom.txt",
                      "fields_rdefs": "fields_rdefs_custom.txt",
                      "actions_rdefs": "actions_rdefs_custom.txt"}
    PLIST_FILE = "info.plist"
    PACKAGE_NAME = "siren_package"

    FILETYPE_VERSION = 6
    XML_FILETYPE_VERSION = 3

    FN_SEP = ";"

    CREATOR = "Clired/Siren Package"
    DEFAULT_EXT = EXT_SIREN
    DEFAULT_TMP = "siren"
    
    def __init__(self, filename, callback_mess=None, mode="r"):
        if filename is not None:
            filename = os.path.abspath(filename)
            if mode != "w" and not os.path.isfile(filename):
                raise IOError("File does not exist")
            if mode != "w" and not zipfile.is_zipfile(filename):
                raise IOError("File is of wrong type")
        self.filename = filename
        self.callback_mess = callback_mess
        self.plist = dict(creator=self.CREATOR,
                          filetype_version=self.FILETYPE_VERSION)
        self.tmp_dir = self.getTmpDir()

    def __str__(self):
        return "PACKAGE: %s" % self.filename

    def raiseMess(self):
        if self.callback_mess is not None:
            self.callback_mess()

    def getFilename(self):
        return self.filename

    def getPackagename(self):
        return self.plist.get("package_name")

    @classmethod
    def getTmpDir(tcl):
        return TmpDirManager.getTmpDir(tcl.DEFAULT_TMP)
    def getPackTmpDir(self):
        return self.tmp_dir

    
    def getFormatV(self):
        return self.plist.get("filetype_version", -1)

    def isOldXMLFormat(self):
        return self.getFormatV() <= self.XML_FILETYPE_VERSION

    def isLatestFormat(self):
        return self.getFormatV() == self.FILETYPE_VERSION

    def getSaveFilename(self):
        svfilename = self.filename
        if self.isOldXMLFormat():
            parts = self.filename.split(".")
            if len(parts) == 1:
                svfilename += "_new"
            elif len(parts) > 1:
                svfilename = ".".join(parts[:-1]) + "_new." + parts[-1]
        return svfilename

    def getNamelist(self):
        return self.package.namelist()

    def closePack(self):
        if self.package is not None:
            self.package.close()
            self.package = None

    def openPack(self):
        try:
            self.package = zipfile.ZipFile(self.filename, "r")
            # plist_fd = self.package.open(self.PLIST_FILE, "r")
            # self.plist = plistlib.readPlist(plist_fd)
            plist_ffd = self.package.read(self.PLIST_FILE)
            self.plist = plistlib.loads(plist_ffd)
        except Exception:
            self.package = None
            self.plist = {}
            self.raiseMess()
            raise

    # READING ELEMENTS
    ##########################

    def read(self, preferences_reader, options_args=None, pv=None, params_only=False):
        elements_read = {}
        self.openPack()
        try:
            preferences = self.readPreferences(preferences_reader, options_args, pv)
            
            if preferences is not None:
                elements_read["preferences"] = preferences

            if not params_only:
                if "actions_rdefs" in self.plist:
                    ar_fns = []
                    for f in self.plist["actions_rdefs"].split(Package.FN_SEP):
                        ff = f.strip()
                        if len(ff) > 0:
                            ar_fns.append(self.package.open(ff))
                    if len(ar_fns) > 0:
                        AR = ActionsRegistry(ar_fns)
                        if "preferences" in elements_read:
                            elements_read["preferences"]["AR"] = AR
                        else:
                            elements_read["preferences"] = {"AR": AR}
                    for fn in ar_fns:
                        fn.close()

                for fields_wich, fields_class in [("fields_vdefs", ColM), ("fields_rdefs", Redescription)]:
                    if fields_wich in self.plist:
                        fields_fns = []
                        for f in self.plist[fields_wich].split(Package.FN_SEP):
                            ff = f.strip()
                            if len(ff) > 0:
                                fields_fns.append(self.package.open(ff))
                        fields_class.extendRP(fields_fns)
                        for fn in fields_fns:
                            fn.close()

                add_info = IOTools.getDataAddInfo(preferences, plist=self.plist, version=self.getFormatV())
                data = self.readData(add_info)
                if data is not None:
                    if "ext_keys" in self.plist:
                        ext_keys = self.plist["ext_keys"].strip().split(Package.FN_SEP)
                        params_collect = data.loadExtensions(ext_keys=ext_keys, filenames=self.plist, params=preferences, details={"package": self.package})
                        preferences.update(params_collect)
                        data.recomputeCols()
                    elements_read["data"] = data
                    reds = self.readRedescriptions(data)
                    if reds is not None and len(reds) > 0:
                        elements_read["reds"] = reds
        finally:
            self.closePack()
        return elements_read

    def readPreferences(self, preferences_reader, options_args=None, pv=None):
        # Load preferences
        preferences = None
        if "preferences_filename" in self.plist:
            fd = None
            try:
                fd = self.package.open(self.plist["preferences_filename"], "r")
                preferences, _, _ = preferences_reader.getPreferences(options_args, conf_file=fd, pv=pv)
            except Exception:
                self.raiseMess()
                raise
            finally:
                if fd is not None:
                    fd.close()
        return preferences

    def readData(self, add_info):
        data = None
        if add_info is None:
            add_info = [{}, Data.NA_str]
        # Load data
        if "data_LHS_filename" in self.plist:
            try:
                fdLHS = io.TextIOWrapper(io.BytesIO(self.package.read(self.plist["data_LHS_filename"])))
                if self.plist.get("data_RHS_filename", self.plist["data_LHS_filename"]) != self.plist["data_LHS_filename"]:
                    fdRHS = io.TextIOWrapper(io.BytesIO(self.package.read(self.plist["data_RHS_filename"])))
                else:
                    fdRHS = None
                data = Data([fdLHS, fdRHS]+add_info, "csv")
            except Exception:
                data = None
                self.raiseMess()
                raise
            finally:
                fdLHS.close()
                if fdRHS is not None:
                    fdRHS.close()
        return data

    def readRedescriptions(self, data):
        reds = []
        # Load redescriptions
        rp = Redescription.getRP()
        if "redescriptions_filename" in self.plist:
            for file_red in self.plist["redescriptions_filename"].split(Package.FN_SEP):
                sid = make_sid(len(reds))
                try:
                    fd = io.TextIOWrapper(io.BytesIO(self.package.read(file_red)))
                    # fd = self.package.open(file_red, "r")
                    rs, _ = rp.parseRedList(fd, data, sid=sid)
                except Exception:
                    self.raiseMess()
                    raise
                finally:
                    fd.close()
                reds.append({"items": rs, "src": ("file", file_red, 1)})
        return reds

    # WRITING ELEMENTS
    ##########################
    # The saving function
    def writeToFile(self, filename, contents):
        # Store old package_filename
        old_package_filename = self.filename
        self.filename = os.path.abspath(filename)
        # Get a temp folder
        tmp_dir = self.getTmpDir()

        # Write plist
        plist, filens = self.makePlistDict(contents)
        try:
            plistlib.dump(plist, open(os.path.join(tmp_dir, self.PLIST_FILE), "bw"))
        except IOError:
            shutil.rmtree(tmp_dir)
            self.filename = old_package_filename
            self.raiseMess()
            raise

        # Write data files
        if "data" in contents:
            try:
                filenames = [os.path.join(tmp_dir, plist["data_LHS_filename"]), None]
                if plist.get("data_RHS_filename", plist["data_LHS_filename"]) != plist["data_LHS_filename"]:
                    filenames[1] = os.path.join(tmp_dir, plist["data_RHS_filename"])
                IOTools.writeData(contents["data"], filenames, toPackage=True)
                IOTools.writeDataExtensions(contents["data"], plist, tmp_dir)
            except IOError:
                shutil.rmtree(tmp_dir)
                self.filename = old_package_filename
                self.raiseMess()
                raise

        # Write redescriptions
        if "redescriptions" in contents:
            for rs in contents["redescriptions"]:
                try:
                    IOTools.writeRedescriptions(rs.get("items", []), os.path.join(tmp_dir, os.path.basename(rs["src"][1])),
                                                names=False, with_disabled=True, toPackage=True)
                except IOError:
                    shutil.rmtree(tmp_dir)
                    self.filename = old_package_filename
                    self.raiseMess()
                    raise

        # Write preferences
        if "preferences" in contents:
            try:
                IOTools.writePreferences(contents["preferences"], contents["preferences_reader"],
                                         os.path.join(tmp_dir, plist["preferences_filename"]), toPackage=True)
            except IOError:
                shutil.rmtree(tmp_dir)
                self.filename = old_package_filename
                self.raiseMess()
                raise

        for k in self.FDEFS_FILENAME.keys():
            if k in contents:
                fn = os.path.join(tmp_dir, plist[k])
                try:
                    with open(fn, "w") as f:
                        f.write(contents[k])
                except IOError:
                    shutil.rmtree(tmp_dir)
                    self.filename = old_package_filename
                    self.raiseMess()
                    raise

        # All"s there, so pack
        try:
            package = zipfile.ZipFile(self.filename, "w")
            package.write(os.path.join(tmp_dir, self.PLIST_FILE),
                          arcname=os.path.join(".", self.PLIST_FILE))
            for eln, element in filens.items():
                package.write(os.path.join(tmp_dir, element),
                              arcname=os.path.join(".", element),
                              compress_type=zipfile.ZIP_DEFLATED)
        except Exception:
            shutil.rmtree(tmp_dir)
            self.filename = old_package_filename
            self.raiseMess()
            raise
        finally:
            package.close()

        # All"s done, delete temp file
        shutil.rmtree(tmp_dir)

    def makePlistDict(self, contents):
        """Makes a dict to write to plist."""
        d = dict(creator=self.CREATOR,
                 filetype_version=self.FILETYPE_VERSION)

        if self.filename is None:
            d["package_name"] = self.PACKAGE_NAME
        else:
            (pn, suffix) = os.path.splitext(os.path.basename(self.filename))
            if len(pn) > 0:
                d["package_name"] = pn
            else:
                d["package_name"] = self.PACKAGE_NAME

        fns = {}
        if "data" in contents:
            d["NA_str"] = contents["data"].NA_str
            fns["data_LHS_filename"] = self.DATA_FILENAMES[0]
            if not contents["data"].isSingleD():
                fns["data_RHS_filename"] = self.DATA_FILENAMES[1]
            ext_keys = contents["data"].getActiveExtensionKeys()
            if len(ext_keys) > 0:
                d["ext_keys"] = Package.FN_SEP.join(ext_keys)
            fns.update(contents["data"].getExtensionsActiveFilesDict())

        if "preferences" in contents:
            fns["preferences_filename"] = self.PREFERENCES_FILENAME
        for k, fn in self.FDEFS_FILENAME.items():
            if k in contents:
                fns[k] = fn
        d.update(fns)

        if "redescriptions" in contents and len(contents["redescriptions"]) > 0:
            base_names = [os.path.basename(c["src"][1]) for c in contents["redescriptions"]]
            d["redescriptions_filename"] = Package.FN_SEP.join(base_names)
            for ci, c in enumerate(base_names):
                fns["redescriptions_filename_%d" % ci] = c
        return d, fns


class IOTools:
    NA_FILETYPE_VERSION = 4
    map_data_params = [{"trg": 0, "from": "delim_in", "to": "delimiter",
                        "vmap": {"(auto)": None, "TAB": "\t", "SPC": " "}},
                       {"trg": 1, "from": "NA_str", "to": "NA_str"},
                       {"trg": 1, "from": "time_dayfirst", "to": "time_dayfirst",
                        "vmap": {"(auto)": None, "yes": True, "no": False}},
                       {"trg": 1, "from": "time_yearfirst", "to": "time_yearfirst",
                        "vmap": {"(auto)": None, "yes": True, "no": False}}]

    @classmethod
    def getDataAddInfo(tcl, params={}, plist={}, version=None, add_info=None):
        if add_info is None:
            add_info = [{}, {"NA_str": Data.NA_str_def}]

        for p in tcl.map_data_params:
            for src in [params, plist]:
                if p["from"] in src:
                    val = src[p["from"]]
                    if "vmap" in p:
                        val = p["vmap"].get(val, val)
                        if val is not None:
                            add_info[p["trg"]][p["to"]] = val
                    else:
                        add_info[p["trg"]][p["to"]] = val
        if add_info[1]["NA_str"] is None and (version is not None and version <= tcl.NA_FILETYPE_VERSION):
            add_info[1]["NA_str"] = Data.NA_str_def
        # print("ADD_INFO", add_info)
        return add_info

    @classmethod
    def writeRedescriptions(tcl, reds, filename, names=[None, None], with_disabled=False, toPackage=False, style="", full_supp=False, nblines=1, supp_names=None, modifiers={}, fmts=[None, None, None]):
        if names is False:
            names = [None, None]
        red_list = [red for red in reds if red.isEnabled() or with_disabled]
        if toPackage:
            fields_supp = [-1, ":extra:status"]
        else:
            fields_supp = None
        rp = Redescription.getRP()
        if filename != "+" and filename != "-":
            f = open(filename, mode="w")
        else:
            f = sys.stdout
        if style == "tex":
            f.write(rp.printTexRedList(red_list, names, fields_supp, nblines=nblines, modifiers=modifiers, fmts=fmts))
        else:
            f.write(rp.printRedList(red_list, names, fields_supp, full_supp=full_supp, supp_names=supp_names, nblines=nblines, modifiers=modifiers, fmts=fmts))
        if filename != "+" and filename != "-":
            f.close()

    @classmethod
    def writeRedescriptionsFmt(tcl, reds, filename, data):
        rp = Redescription.getRP()
        params = tcl.getPrintParams(filename, data)
        params["modifiers"] = rp.getModifiersForData(data)
        tcl.writeRedescriptions(reds, filename, **params)

    @classmethod
    def writePreferences(tcl, preferences, preferences_reader, filename, toPackage=False, conf_filter=None, sections=False, helps=False, defaults=False, only_core=False):
        with open(filename, "w") as f:
            f.write(preferences_reader.dispParameters(preferences, conf_filter, sections, helps, defaults, only_core))

    @classmethod
    def writeData(tcl, data, filenames, toPackage=False):
        data.writeCSV(filenames)

    @classmethod
    def writeDataExtensions(tcl, data, plist=None, tmp_dir="./"):
        if plist is not None:
            data.saveExtensions(plist, {"tmp_dir": tmp_dir})

    @classmethod
    def saveAsPackage(tcl, filename, data, preferences=None, preferences_reader=None, reds=None, AR=None):
        package = Package(None, None, mode="w")

        (filename, suffix) = os.path.splitext(filename)
        contents = {}
        if preferences is not None:
            if preferences_reader is None:
                preferences_reader = getPreferencesReader()
            contents["preferences"] = dict(preferences)
            contents["preferences_reader"] = preferences_reader
        if data is not None:
            contents["data"] = data
        if reds is not None and len(reds) > 0:
            # contents["preferences"].pop("queries_file")
            # contents["preferences"].pop("queries_in")
            if isinstance(reds[0], Redescription):
                contents["redescriptions"] = [{"items": reds, "src": ('file', Package.REDESCRIPTIONS_FILENAME, 1)}]
            else:
                contents["redescriptions"] = reds

        # definitions
        vdefs = ColM.getRP().fieldsToStr()
        if len(vdefs) > 0:
            contents["preferences"].pop("fields_vdefs")
            contents["fields_vdefs"] = vdefs
        rdefs = Redescription.getRP().fieldsToStr()
        if len(rdefs) > 0:
            contents["preferences"].pop("fields_rdefs")
            contents["fields_rdefs"] = rdefs
        if AR is not None:
            adefs = AR.actionsToStr()
            if len(adefs) > 0:
                contents["preferences"].pop("actions_rdefs")
                contents["actions_rdefs"] = adefs

        package.writeToFile(filename+suffix, contents)

    @classmethod
    def getPrintParams(tcl, filename, data=None):
        basename = os.path.basename(filename)
        params = {"with_disabled": False, "style": "", "full_supp": False, "nblines": 1,
                  "names": [None, None], "supp_names": None}

        named = re.search("[^a-zA-Z0-9]named[^a-zA-Z0-9]", basename) is not None
        supp_names = (re.search("[^a-zA-Z0-9]suppnames[^a-zA-Z0-9]", basename) is not None) or \
                     (re.search("[^a-zA-Z0-9]suppids[^a-zA-Z0-9]", basename) is not None)

        params["with_disabled"] = re.search("[^a-zA-Z0-9]all[^a-zA-Z0-9]", basename) is not None
        params["full_supp"] = (re.search("[^a-zA-Z0-9]support[^a-zA-Z0-9]", basename) is not None) or supp_names

        if re.search(".tex$", basename):
            params["style"] = "tex"

        tmp = re.search("[^a-zA-Z0-9](?P<nbl>[1-3]).[a-z]*$", basename)
        if tmp is not None:
            params["nblines"] = int(tmp.group("nbl"))

        if named and data is not None:
            params["names"] = data.getNames()
            params["fmts"] = data.getFmts()
        if supp_names:
            params["supp_names"] = data.getRNames()
        return params

    @classmethod
    def isWritable(tcl, fname):
        if not os.path.isfile(fname):
            try:
                tfs = open(fname, "a")
                tfs.close()
                os.remove(fname)
            except IOError:
                return False
        return True

    @classmethod
    def checkQueriesFilename(tcl, fname, which="in", queries_basic=None, src_folder=None):
        if not tcl.isBaseFName(fname) or (queries_basic is not None and fname == queries_basic):
            queries_fname = fname
        else:
            queries_fname = tcl.makeUpAltFName(fname, queries_basic)

        # make absolute, so it's meaningful to check existence...
        if src_folder is not None and tcl.isAbsFName(src_folder) and tcl.isProperFName(queries_fname) and not tcl.isAbsFName(queries_fname):
            queries_fname = os.path.normpath(os.path.join(src_folder, queries_fname))

        if tcl.isProperFName(queries_fname):
            if os.path.isfile(queries_fname):
                if which == "in" or which == "?":
                    return (queries_fname, "in", None)
                else:
                    return (queries_fname, "out", "Queries file [%s] will be overwritten..." % queries_fname)
            else:
                if which == "in":
                    return (None, None, "Queries file [%s] not found..." % queries_fname)
                else:
                    if tcl.isWritable(queries_fname):
                        return (queries_fname, "out", None)
                    else:
                        return (None, None, "Queries file [%s] not writable..." % queries_fname)
        if queries_fname == "-":
            return (None, None, None)
        return (None, None, "Queries file [%s] has improper name..." % queries_fname)

    @classmethod
    def prepareFilenames(tcl, params, tmp_dir=None, src_folder=None, queries_basic_dest="out"):
        filenames = {"style_data": "csv",
                     "add_info": tcl.getDataAddInfo(params)
                     }

        # filenames from command-line
        params_filenames = params.get("filename", {})
        for p in ["result_rep", "data_rep", "extensions_rep"]:
            if p not in params:
                params[p] = ""
            params[p] = tcl.expandPath(params[p], tmp_dir, src_folder)
            
        # Make data file names
        filenames["LHS_data"] = ""
        if len(params_filenames.get("data_file", [])) > 0:
            filenames["LHS_data"] = params_filenames["data_file"][0]
        elif len(params.get("LHS_data", "")) != 0:
            filenames["LHS_data"] = params["LHS_data"]
        elif len(params.get("data_l", "")) != 0:
            filenames["LHS_data"] = os.path.join(params["data_rep"], params["data_l"]+params.get("ext_l", ""))

        filenames["RHS_data"] = ""
        if len(params_filenames.get("data_file", [])) > 0:
            if len(params_filenames.get("data_file")) == 1:
                filenames["RHS_data"] = params_filenames["data_file"][0]  # only one data file: use for both sides
            else:
                filenames["RHS_data"] = params_filenames["data_file"][1]  # more than one data file: first to RHS
                if len(params_filenames["data_file"]) > 2:  # if there is more than two data files provided on the command lines, consider the surplus as queries files
                    params_filenames["queries_file"] = params_filenames["data_file"][2:]+params_filenames.get("queries_file", [])  # the rest, queries
        elif len(params.get("RHS_data", "")) != 0:
            filenames["RHS_data"] = params["RHS_data"]
        elif len(params.get("data_r", "")) != 0:
            filenames["RHS_data"] = os.path.join(params["data_rep"], params["data_r"]+params.get("ext_r", ""))

        if len(params.get("trait_data", "")) != 0:
            filenames["traits_data"] = params["traits_data"]
        elif len(params.get("data_t", "")) != 0:
            filenames["traits_data"] = os.path.join(params["data_rep"], params["data_t"]+params.get("ext_t", ""))

        if os.path.splitext(filenames["LHS_data"])[1] != ".csv" or os.path.splitext(filenames["RHS_data"])[1] != ".csv":
            filenames["style_data"] = "multiple"
            filenames["add_info"] = []

        if len(params.get("extensions_names", "")) != 0:
            filenames["extensions"] = {}
            extkf = params.get("extensions_names", "")
            for blck in extkf.strip().split(Package.FN_SEP):
                parts = [p.strip() for p in blck.split("=")]
                if len(parts) == 2:
                    filenames["extensions"]["extf_"+parts[0]] = params["extensions_rep"] + parts[1]
                    
        # Prepare queries file names

        # First basic
        queries_basic = "-"
        queries_fnames = {"in": [], "out": []}
        if len(params.get("queries_file", "")) != 0:
            queries_basic = params["queries_file"]
        elif params.get("out_base", "-") != "-" and len(params["out_base"]) > 0 and tcl.isBaseFName(params.get("out_base", "./")):
            queries_basic = os.path.join(params["result_rep"], params["out_base"]+params.get("ext_queries", ".queries"))

        queries_fname, tow, mess = tcl.checkQueriesFilename(queries_basic, which=queries_basic_dest, queries_basic=queries_basic, src_folder=src_folder)
        if mess is not None:
            print(mess)
        if queries_fname is not None:
            queries_fnames[tow].append(queries_fname)

        cands_qnames = [("in", params.get("queries_in", "").split(Package.FN_SEP)),
                        ("out", params.get("queries_out", "").split(Package.FN_SEP)),
                        ("?", params_filenames.get("queries_file", []))]

        # Prepare additional query file names        
        for wio, cands in cands_qnames:
            qbasic = queries_basic
            for pi, org_part in enumerate(cands):
                if org_part != "":
                    part = tcl.expandPath(org_part, tmp_dir, params["result_rep"])
                    p_base, p_ext = os.path.splitext(part)                    
                    if len(p_ext) == 0:
                        part += params.get("ext_queries", ".queries")

                    if pi == 0:
                        queries_fname, tow, mess = tcl.checkQueriesFilename(part, which=wio, src_folder=src_folder)
                        if queries_fname is not None:
                            if len(queries_fnames[tow]) > 0 and queries_fnames[tow][-1] == qbasic:
                                queries_fnames[tow].pop()
                            qbasic = queries_fname
                        else:
                            queries_fname, tow, mess = tcl.checkQueriesFilename(part, which=wio, queries_basic=qbasic, src_folder=src_folder)
                    else:
                        queries_fname, tow, mess = tcl.checkQueriesFilename(part, which=wio, queries_basic=qbasic, src_folder=src_folder)
                    if mess is not None:
                        print(mess)
                    if queries_fname is not None:
                        queries_fnames[tow].append(queries_fname)

        filenames["queries_in"] = queries_fnames["in"]
        filenames["queries_out"] = queries_fnames["out"]

        if tcl.isProperFName(queries_basic):
            head_q, tail_q = os.path.split(queries_basic)
            root_q, ext_q = os.path.splitext(tail_q)
            basic_parts = (head_q, root_q, ext_q)
            filenames["basis"] = os.path.join(head_q, root_q)
        else:
            basic_parts = None
            filenames["basis"] = ""

        build_files = [("queries_named", "queries_named_file", "_named", None),
                       ("support", "support_file", "_support", params.get("ext_support")),
                       ("logfile", "logfile", "_"+params["task"] if params.get("task", "mine") != "mine" else "", params.get("ext_log")),
                       ("pairs_store", "pairs_store", "_pairs", params.get("ext_log")),
                       ("pairs_traces", "pairs_traces", "_pairs_traces", params.get("ext_log"))]
        # Make named queries file name
        filenames["logfile"] = "-"
        for (path_param, file_param, default_suff, default_ext) in build_files:
            fpath = make_filepath(params.get(file_param),  basic_parts, default_suff, default_ext)
            if fpath is not None:
                filenames[path_param] = fpath

        if len(params.get("series_id", "")) > 0:
            for k in filenames.keys():
                if type(filenames[k]) is str:
                    filenames[k] = filenames[k].replace("__SID__", params["series_id"])

        # make paths absolute
        if src_folder is not None and tcl.isAbsFName(src_folder):
            pk_filelist = ["queries_in", "queries_out"]
            pk_exclude = ["style_data", "add_info", "extensions"] + pk_filelist
            pk_file = [pk for pk in list(filenames.keys()) + list(filenames.get("extensions", {}).keys()) if pk not in pk_exclude]

            for k in pk_file:
                filenames[k] = tcl.expandPath(filenames[k], tmp_dir, src_folder)
                
        return filenames

    @classmethod
    def outputResults(tcl, filenames, results, data=None, with_headers=True, mode="w", data_recompute=None):
        default_stdout = tcl.outputResultsFList(filenames, results, data)
        tcl.outputResultsQNS(filenames, results, data, with_headers, mode, data_recompute, default_stdout)

    @classmethod
    def outputResultsFList(tcl, filenames, results, data=None):
        default_stdout = True
        for fname in filenames.get("queries_out", []):
            default_stdout = False
            IOTools.writeRedescriptionsFmt(results, fname, data)
        return default_stdout

    @classmethod
    def outputResultsQNS(tcl, filenames, results, data=None, with_headers=True, mode="w", data_recompute=None, default_stdout=True):
        filesfp = {"queries": None, "queries_named": None, "support": None}
        if "queries" in filenames or default_stdout:
            if filenames.get("queries", "-") == "-":
                filesfp["queries"] = sys.stdout
            else:
                filesfp["queries"] = open(filenames["queries"], mode)

        if data is not None and data.hasNames() and "queries_named" in filenames:
            filesfp["queries_named"] = open(filenames["queries_named"], mode)

        if "support" in filenames:
            filesfp["support"] = open(filenames["support"], mode)

        if all([v is None for v in filesfp.values()]):
            # nowhere to print
            return

        rp = Redescription.getRP()
        modifiers, modifiers_recompute = {}, {}
        fstyle = "basic"
        header_recompute = ""
        names = None
        if data is not None:
            names = data.getNames()
            modifiers = rp.getModifiersForData(data)

        all_fields = rp.getListFields(fstyle, modifiers)

        if data_recompute is not None:
            modifiers_recompute = rp.getModifiersForData(data_recompute)
            fields_recompute = rp.getListFields("stats", modifiers_recompute)
            header_recompute = "\t" + rp.dispHeaderFields(fields_recompute) + "\tacc_diff"

        if with_headers:
            if filesfp["queries"] is not None:
                filesfp["queries"].write(rp.dispHeaderFields(all_fields)+header_recompute+"\n")
            if filesfp["queries_named"] is not None:
                filesfp["queries_named"].write(rp.dispHeaderFields(all_fields)+header_recompute+"\n")

        # TO DEBUG: output all shown in siren, i.e. no filtering
        addto = ""
        for org in results:
            if data_recompute is not None:
                red = org.copy()
                red.recompute(data_recompute)
                acc_diff = (red.getAcc()-org.getAcc())/org.getAcc()
                addto = "\t"+red.disp(list_fields=fields_recompute)+"\t%f" % acc_diff
            if filesfp["queries"] is not None:
                filesfp["queries"].write(org.disp(list_fields=all_fields)+addto+"\n")
            if filesfp["queries_named"] is not None:
                filesfp["queries_named"].write(org.disp(names, list_fields=all_fields)+addto+"\n")
            if filesfp["support"] is not None:
                filesfp["support"].write(org.dispSupp()+"\n")

        for (ffi, ffp) in filesfp.items():
            if ffp is not None and filenames.get(ffi, "") != "-" and ffp != sys.stdout:
                ffp.close()

    @classmethod
    def loadUnpackaged(tcl, params, src_folder, pack_filename=None, tmp_dir=None, task_load={}, data=None, reds=None):
        filenames = tcl.prepareFilenames(params, tmp_dir, src_folder, task_load.get("queries_basic_dest", "out"))
        if pack_filename is not None:
            filenames["package"] = os.path.abspath(pack_filename)
        
        if task_load.get("with_log", True):
            logger = Log(verbosity=params["verbosity"], output=filenames["logfile"])
        else:
            logger = None

        if data is None and os.path.exists(filenames["LHS_data"]) and os.path.exists(filenames["RHS_data"]):
            data = Data([filenames["LHS_data"], filenames["RHS_data"]]+filenames["add_info"], filenames["style_data"])
            data.loadExtensions(ext_keys=params.get("activated_extensions", []), filenames=filenames.get("extensions"), params=params)

        ## ACTIONS AND FIELDS
        if len(params.get("actions_rdefs", "").strip()) != 0:
            ar_fns = []
            for f in params["actions_rdefs"].split(Package.FN_SEP):
                ff = tcl.expandPath(f.strip(), tmp_dir, src_folder)
                if len(ff) > 0:
                    ar_fns.append(ff)
            if len(ar_fns) > 0:
                if params.get("AR") is not None:
                    params["AR"].extend(ar_fns)
                else:
                    AR = ActionsRegistry(ar_fns)
                    params["AR"] = AR

        for fields_wich, fields_class in [("fields_vdefs", ColM), ("fields_rdefs", Redescription)]:
            if len(params.get(fields_wich, "").strip()) != 0:
                fields_fns = []
                for f in params[fields_wich].split(Package.FN_SEP):
                    ff = tcl.expandPath(f.strip(), tmp_dir, src_folder)
                    if len(ff) > 0:
                        fields_fns.append(ff)
                fields_class.extendRP(fields_fns)

        ## REDESCRIPTIONS
        if (data is not None) and len(filenames.get("queries_in", [])) > 0:
            off_rlists = 0
            if reds is None:
                reds = []
            if type(reds) is int: ##
                off_rlists = reds
                reds = []

            rp = Redescription.getRP()
            for file_red in filenames["queries_in"]:
                sid = make_sid(len(reds) + off_rlists)
                # attempt loading redescriptions from file
                try:
                    with open(file_red) as fd:
                        rs, _ = rp.parseRedList(fd, data, sid=sid)  # avoid rid collisions
                except (IOError, DataError):
                    # not able to load redescriptions from the file
                    rs = []
                if len(rs) > 0:
                    reds.append({"items": rs, "src": ("file", file_red, 0)})
        elif type(reds) is int: ##
            reds = None # offset for reading reds, but there aren't any to read
                    
        if logger is not None and task_load.get("log", True):
            logger.printL(2, data, "log")

        return {"params": params, 
                "data": data, "reds": reds,
                "logger": logger, "filenames": filenames}

    @classmethod
    def loadAll(tcl, arguments=[], conf_defs=None, tasks=[], tasks_default=None, tasks_load={}, conf_filter=None):
        preferences_reader = getPreferencesReader(conf_defs, tasks, tasks_default)
        params, leftover_args, preferences_mod = preferences_reader.getPreferences(arguments)

        if "task" in params and params["task"] in tasks_load:
            task_load = tasks_load[params["task"]]
        else:
            task_load = tasks_load
        params_only = task_load.get("params_only", False)

        if "task" in params and params["task"] in tasks_load:
            if conf_filter is None:
                conf_filter = task_load.get("conf_defs")
            else:
                conf_filter = conf_filter + task_load.get("conf_defs", [])
        else:
            conf_defs = None
        if params_only and len(preferences_mod) == 0:
            params["usage"] = True
        msg = preferences_reader.processHelp(params, conf_filter=conf_filter)
        if msg is not None:
            return False, msg

        config_filename = None
        package = None
        pack_filename = None
        elements_read = {}
        leftover_args = None
        
        tmp_dir = None
        exec_folder = os.path.dirname(os.path.abspath(__file__))
        if "src_folder" in params.get("filename", {}):
            src_folder = params["filename"]["src_folder"]
        else:
            src_folder = exec_folder

        if "pack_file" in params.get("filename", {}):
            pack_filename = params["filename"]["pack_file"][0]
            package = Package(pack_filename)
            elements_read = package.read(preferences_reader, params_only=params_only, pv={})
            tmp_dir = package.getPackTmpDir()
            for k,v in elements_read.get("preferences", {}).items():
                if k not in preferences_mod:
                    params[k] = v
                    
        loaded = {"params": params, "leftover_args": leftover_args,
                  "preferences_reader": preferences_reader, "preferences_mod": preferences_mod}
            
        if not params_only:
            loaded.update(tcl.loadUnpackaged(params, src_folder, pack_filename, tmp_dir, task_load, data=elements_read.get("data"), reds=elements_read.get("reds")))
            loaded["package"] = package
        return True, loaded

    @classmethod
    def collectLoadedReds(tcl, loaded):
        # pull all reds together
        all_reds = []
        srcs_reds = []
        all_queries_src = {}
        if loaded["reds"] is not None:
            count_red_lists = len(loaded["reds"])
            for i in range(len(loaded["reds"])):
                src = loaded["reds"][i]["src"][1]
                all_queries_src[src] = {"pos": i}
                all_queries_src[src]["nb_reds"] = len(loaded["reds"][i]["items"])
                all_reds.extend(loaded["reds"][i]["items"])
                srcs_reds.append(src)

        return all_reds, srcs_reds, all_queries_src

    @classmethod
    def isFSuff(tcl, qfname):
        return re.match(PATT_QSUFF, qfname) is not None

    @classmethod
    def isProperFName(tcl, qfname):
        return (qfname is not None and len(qfname.strip()) > 0 and qfname != "-" and not tcl.isFSuff(qfname))

    @classmethod
    def isAbsFName(tcl, qfname):
        return os.path.abspath(qfname) == qfname

    @classmethod
    def isBaseFName(tcl, qfname):
        return os.path.basename(qfname) == qfname

    @classmethod
    def isRelFName(tcl, qfname):  # start with ./ or such
        return not tcl.isBaseFName(qfname) and not tcl.isAbsFName(qfname)

    @classmethod
    def makeUpAltFName(tcl, qaltfname, qbasicfname):
        tmp = re.match(PATT_QSUFF, qaltfname)
        if tmp is not None:
            return tcl.makeUpFName(qbasicfname, tmp.group("suff"))
        return qaltfname

    @classmethod
    def makeUpFName(tcl, qfilename=None, alt_suff="_XXX"):
        if qfilename is None:
            trg = alt_suff
        elif qfilename == "-":
            trg = "-"
        else:
            qroot, qext = os.path.splitext(qfilename)
            if "." in alt_suff:
                trg = qroot + alt_suff
            else:
                trg = qroot + alt_suff + qext
        return trg

    @classmethod
    def expandPath(tcl, path, tmp_dir=None, src_folder=None):
        exp_path = path
        if tcl.isProperFName(path) and sys.platform != "win32":
            if sys.platform != "darwin":
                exp_path = re.sub("~", os.path.expanduser("~"), exp_path)

            if TmpDirManager.isTmpDirToken(exp_path):
                if tmp_dir is None:
                    tmp_dir = Package.getTmpDir()
                exp_path = os.path.join(tmp_dir, "")
            elif src_folder is not None and not tcl.isAbsFName(exp_path):
                exp_path = os.path.normpath(os.path.join(src_folder, exp_path))
        return exp_path
    
