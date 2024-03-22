import pdb
import os.path
import io
import numpy

pck_imports = True

try:
    from classCol import DataError, ColM, BoolColM, CatColM, NumColM
    from classRedescription import Redescription
    from classSParts import SSetts
    import csv_reader
    # import prepare_polygons
except ModuleNotFoundError:
    from .classCol import DataError, ColM, BoolColM, CatColM, NumColM
    from .classRedescription import Redescription
    from .classSParts import SSetts
    from . import csv_reader
    # from . import prepare_polygons
    pck_imports = False


class DataExtension(object):

    extension_key = "-"
    extras_map = {}
    filing = []
    params_keys = []

    def __init__(self, data, filenames=None, params=None, details={}):
        self.parent_data = data
        self.own_data = {}
        self.params_changed = set()
        pp = {}
        if details is not None and "params" in details:
            pp.update(details["params"])
        if params is not None:
            pp.update(params)
        if len(pp) > 0:
            self.setParams(pp)
        if filenames is not None:
            self.doWithFiles("load", filenames, details)

    @classmethod
    def getKey(tcl):
        return tcl.extension_key

    @classmethod
    def getParamKeys(tcl):
        return tcl.params_keys

    def getParams(self):
        params_kv = {}
        for k in self.getParamKeys():
            if self.hasElement(k):
                params_kv[k] = self.getElement(k)
        return params_kv

    @classmethod
    def getExtensionDetails(tcl):
        return {}

    @classmethod
    def getExtrasKeys(tcl):
        return list(tcl.extras_map.keys())

    @classmethod
    def getFilenamesKeys(tcl):
        fks = []
        for f in tcl.filing:
            fks.extend(f.get("filenames", {}).keys())
        return fks

    @classmethod
    def getFilesDict(tcl):
        fks = {}
        for f in tcl.filing:
            fks.update(f.get("filenames", {}))
        return fks

    def getActiveFilesDict(self):
        fks = {}
        for f in self.filing:
            if ("active" not in f) or eval("self.%s()" % f["active"]):
                fks.update(f.get("filenames", {}))
        return fks

    def getParentData(self):
        return self.parent_data

    def hasElement(self, key):
        return key in self.own_data

    def getElement(self, key, default=None):
        return self.own_data.get(key, default)

    def setElement(self, key, elem):
        self.own_data[key] = elem

    def delElement(self, key):
        if self.hasElement(key):
            del self.own_data[key]

    def updateElements(self, vs={}):
        self.own_data.update(vs)

    def resetElements(self, replacement={}):
        self.own_data = replacement

    def retainElements(self, keys=[]):
        replacement = {}
        for k in keys:
            if self.hasElement(k):
                replacement[k] = self.getElement(k)
        self.resetElements(replacement)

    def hasChangedParams(self, pks):
        for pk in pks:
            if pk in self.params_changed:
                return True
        return False

    def checkedParams(self, pks):
        for pk in pks:
            self.params_changed.discard(pk)

    def setParams(self, params):
        if "params" in params:
            params = params.get("params")
        for k in self.params_keys:
            if k in params and self.getElement(k) != params[k]:
                self.setElement(k, params[k])
                self.params_changed.add(k)

    def closeFp(self, fp, details={}):
        if details is None or details.get("package") is None:
            fp.close()

    def getFp(self, fk, filenames=None, details={}, mode="r"):
        fn, folder, pck = (None, "", None)
        if filenames is None or fk not in filenames:
            fn = self.getFilesDict().get(fk)
        else:
            fn = filenames[fk]

        if fn is not None:
            if details is not None and len(details) > 0:
                if "tmp_dir" in details:
                    folder = details["tmp_dir"]
                if "dir" in details:
                    folder = details["dir"]
                if "package" in details:
                    pck = details["package"]
            if len(folder) > 0:
                fn = os.path.join(folder, fn)

            if pck is not None:
                if mode == "r":
                    return io.TextIOWrapper(io.BytesIO(pck.read(fn)))
                else:
                    return pck.open(fn, mode)
            try:
                return open(fn, mode)
            except FileNotFoundError:
                return None

    def computeExtras(self, item, extra_keys=None, details={}):
        if extra_keys is None:
            extra_keys = list(self.extras_map.keys())
        out = dict([(e, None) for e in extra_keys])
        map_meths = {}
        for extra_key in extra_keys:
            if extra_key in self.extras_map:
                meth = self.extras_map[extra_key]
                if meth not in map_meths:
                    map_meths[meth] = []
                map_meths[meth].append(extra_key)
        for meth, ekeys in map_meths.items():
            tmp = self.call_meth(meth, item, details)
            for ek in ekeys:
                out[ek] = tmp[ek]
        return out

    def computeExtra(self, item, extra_key=None, details={}):
        out = {}
        if extra_key in self.extras_map:
            out = self.call_meth(self.extras_map[extra_key], item, details)
        return out.get(extra_key)

    def call_meth(self, method_name, item, details={}):
        return getattr(self, method_name)(item, details)

    def doWithFiles(self, action="save", filenames={}, details={}):
        set_fks = None
        if filenames is not None:
            set_fks = set(filenames.keys())
        for f in self.filing:
            ffiles = f.get("filenames", {}).keys()
            if set_fks is None or len(set_fks.intersection(ffiles)) == len(ffiles):
                if f.get(action) is not None:
                    self.call_meth(f.get(action), filenames, details)

try:
    if pck_imports:
        import prepare_polygons
    else:
        from . import prepare_polygons

    class GeoPlusExtension(DataExtension):

        extension_key = "geoplus"
        extras_map = {}
        filing = []
        params_keys = ["hgrid_percentile", "wgrid_percentile", "after_cut", "dst_type"]

        def reset(self):
            self.retainElements(self.params_keys)

        def setCoordsBckg(self, coords_bckg):
            self.reset()
            self.setElement("coords_bckg", coords_bckg)

        def needRecomputeEdges(self):
            loc_params = ["wgrid_percentile", "hgrid_percentile", "dst_type"]
            return not self.hasElement("list_edges") or (self.getElement("after_cut", True) and self.hasChangedParams(loc_params))

        # ######################################
        # #### COMPUTATIONS
        # ######################################
        def gatherPoints(self):
            PointsMap = {}
            PointsIds = {}
            PointsPolys = None
            coords = self.getParentData().getCoordPoints()
            rnames = self.getParentData().getRNames()

            map_merge = None
            key = None
            if self.hasElement("coords_bckg"):
                map_merge = {}
                if self.hasElement("rnames_bckg"):
                    key = "rnames"
                else:
                    key = "coords"

            for i, coord in enumerate(rnames):
                PointsIds[i] = rnames[i]
                PointsMap[i] = (coords[i, 0], coords[i, 1])
                if map_merge is not None:
                    if key == "rnames":
                        map_merge[PointsIds[i]] = i
                    else:
                        map_merge[PointsMap[i]] = i
            org_nb = len(PointsMap)
            added_bckg = False
            if self.hasElement("coords_bckg"):
                i = len(PointsIds)
                rnames = self.getElement("rnames_bckg")
                for ci, coord in enumerate(self.getElement("coords_bckg")):
                    kv = coord
                    rname = "bckg_%d" % i
                    if key == "rnames":
                        kv = rnames[ci]
                        rname = rnames[ci]

                    if kv in map_merge:
                        next_id = map_merge[kv]
                    else:
                        added_bckg = True
                        next_id = i
                        map_merge[kv] = next_id
                        i += 1

                    PointsIds[next_id] = rname
                    PointsMap[next_id] = coord
            ## print("NBS", org_nb, len(PointsMap))

            # if self.getParentData().hasGeoPoly() and not added_bckg:
            #     PointsPolys = self.getParentData().getCoordPoints()

            return PointsMap, PointsIds, PointsPolys

        def computeEdges(self, details={}, force=False):
            if details is None:
                details = {}
            self.setParams(details)
            if force or self.needRecomputeEdges():
                loc_params = ["wgrid_percentile", "hgrid_percentile", "dst_type"]
                self.checkedParams(["after_cut"]+loc_params)
                PointsMap, PointsIds, PointsPolys = self.gatherPoints()
                # TODO: Use data polygons if available, needs finding intersection between sides of polygons
                if PointsPolys is not None:
                    map_edges, list_edges, polys, polys_cut, bbox = prepare_polygons.prepare_edges_polys(PointsMap, PointsPolys)
                else:
                    map_edges, list_edges, polys, polys_cut, bbox = prepare_polygons.prepare_edges_dst(PointsMap, self.getElement("hgrid_percentile", -1), self.getElement("wgrid_percentile", -1), self.getElement("after_cut", True), self.getElement("dst_type", "globe"))
                dets = {"p_ids": PointsIds, "p_map": PointsMap,
                        "map_edges": map_edges, "list_edges": list_edges,
                        "polys": polys, "polys_cut": polys_cut}

                self.updateElements(dets)
            return {"list_edges": self.getElement("list_edges")}

        def prepNodePairs(self, details={}, force=False):
            if details is None:
                details = {}
            self.setParams(details)
            recomputed = False
            if force or self.needRecomputeEdges():
                self.computeEdges(details=details, force=force)
                recomputed = True
            if recomputed or not self.hasElement("node_pairs"):
                node_pairs = prepare_polygons.prepare_node_pairs(self.getElement("list_edges"), self.getParentData().nbRows(), self.getElement("after_cut", True))
                self.setElement("node_pairs", node_pairs)
            return {"node_pairs": self.getElement("node_pairs")}

        def prepEdgesCoordsFlatten(self, seids=None, details={}, force=False):
            if details is None:
                details = {}
            self.setParams(details)
            if force or self.needRecomputeEdges():
                self.computeEdges(details=details, force=force)
            return prepare_polygons.prepare_edges_coords_flatten(self.getElement("list_edges"), seids, after_cut=self.getElement("after_cut", True))

        def prepAreasData(self, cells_colors, details={}, force=False):
            if details is None:
                details = {}
            self.setParams(details)
            recomputed = False
            if force or self.needRecomputeEdges():
                self.computeEdges(details=details, force=force)
                recomputed = True

            if recomputed or not self.hasElement("nodes_graph"):
                nodes_graph = prepare_polygons.prepare_nodes_graph(self.getElement("list_edges"), after_cut=self.getElement("after_cut", True))
                self.updateElements({"nodes_graph": nodes_graph})

            pp = prepare_polygons.prepare_areas_polys(self.getElement("polys"), self.getElement("polys_cut"), self.getElement("after_cut", True))
            if recomputed or not self.hasElement("border_info"):
                border_info = prepare_polygons.prepare_border_info(self.getElement("list_edges"), self.getElement("map_edges"), pp, after_cut=self.getElement("after_cut", True))
                self.updateElements({"border_info": border_info})

            ccs_data, cks, adjacent = prepare_polygons.prepare_areas_data(cells_colors, self.getElement("list_edges"), pp, self.getElement("border_info"), self.getElement("nodes_graph"))
            return {"ccs_data": ccs_data, "cks": cks, "adjacent": adjacent}

        def computeCohesion(self, item, details={}, force=False):
            if details is None:
                details = {}
            cohesion = -1
            back_v = details.get("back_v")
            if isinstance(item, Redescription):
                vec = numpy.array(item.supports().getVectorABCD())
                back_v = SSetts.Eoo
                if details.get("spids") is not None:
                    back_v = 0
                    vec_org = vec
                    vec = numpy.zeros(vec_org.shape, dtype=int)
                    for spid in details.get("spids"):
                        vec[vec_org == spid] = 1
            elif isinstance(item, ColM):
                vec = item.getVector()
                if item.simpleBool():
                    back_v = 0
            elif type(item) is tuple and len(items) >= 2:
                var = self.col(item[0], item[1])
                vec = var.getVector()
                if var.simpleBool():
                    back_v = 0
            else:
                vec = item
            cohesion = self.computeVectorCohesion(vec, back_v, details, force)
            return {"cohesion": cohesion}
        extras_map["cohesion"] = "computeCohesion"

        def computeVectorCohesion(self, vec_org, back_v=None, details={}, force=False):
            if details is None:
                details = {}
            self.setParams(details)
            self.prepNodePairs(details=details, force=force)
            cohesion = -1
            bv = back_v if back_v is not None else float("Inf")
            node_pairs = self.getElement("node_pairs")
            vec = numpy.concatenate([vec_org, [bv, bv]])
            nb_inner_edges = numpy.sum((node_pairs[:, 2] >= 0) & (vec[node_pairs[:, 1]] == vec[node_pairs[:, 2]]))
            nb_outer_edges = numpy.sum(vec[node_pairs[:, 1]] != vec[node_pairs[:, 2]])
            if nb_outer_edges > 0:
                cohesion = nb_inner_edges/float(nb_inner_edges+nb_outer_edges)
            return cohesion

        # ######################################
        # #### SAVING / LOADING
        # ######################################
        F_COORD_BCKG = "coords_bckg"
        filing.append({"filenames": {"extf_"+F_COORD_BCKG: F_COORD_BCKG+".csv"},
                       "save": "writeBckgCoords", "load": "readBckgCoords", "active": "containsBckgCoords"})

        def containsBckgCoords(self):
            return self.hasElement("coords_bckg")

        def writeBckgCoords(self, filenames, details={}):
            if self.hasElement("coords_bckg"):
                fp = self.getFp("extf_"+self.F_COORD_BCKG, filenames, details, "w")
                csv_reader.write_coords_csv(fp, self.getElement("coords_bckg"), self.getElement("rnames_bckg"))
                self.closeFp(fp, details)

        def readBckgCoords(self, filenames, details={}):
            fp = self.getFp("extf_"+self.F_COORD_BCKG, filenames, details, "r")
            coords_bckg, rnames_bckg = csv_reader.read_coords_csv(fp)
            self.closeFp(fp, details)
            if coords_bckg is not None:
                self.setCoordsBckg(coords_bckg)
                if rnames_bckg is not None:
                    self.setElement("rnames_bckg", rnames_bckg)
            return coords_bckg, rnames_bckg

        F_EDGES = "list_edges"
        filing.append({"filenames": {"extf_"+F_EDGES: F_EDGES+".csv"},
                       "save": "writeEdges", "load": "readEdges", "active": "containsEdges"})

        def containsEdges(self):
            return self.hasElement("list_edges")

        def writeEdges(self, filenames, details={}):
            if self.hasElement("list_edges"):
                fp = self.getFp("extf_"+self.F_EDGES, filenames, details, "w")
                prepare_polygons.write_edges(fp, self.getElement("list_edges"))
                self.closeFp(fp, details)

        def readEdges(self, filenames, details={}):
            fp = self.getFp("extf_"+self.F_EDGES, filenames, details, "r")
            if fp is not None:
                ## id_int = details.get("id_int", False)
                map_edges, list_edges, polys, polys_cut = prepare_polygons.read_edges_and_co(fp)
                params_up = {}
                if "grid_percentile" in list_edges[0]:
                    params_up["hgrid_percentile"], params_up["wgrid_percentile"] = list_edges[0]["grid_percentile"]
                if "dst_type" in list_edges[0]:
                    params_up["dst_type"] = list_edges[0]["dst_type"]
                # params_up["after_cut"] = list_edges[0].get("last_org", 0) < list_edges[0].get("last_cut", 0)
                if len(params_up) > 0:
                    self.setParams(params_up)
                if list_edges[0].get("source") == "polys":
                    self.checkedParams(["hgrid_percentile", "wgrid_percentile", "dst_type"])
                else:
                    self.checkedParams(params_up.keys())
                self.closeFp(fp, details)
                dets = {"map_edges": map_edges, "list_edges": list_edges,
                        "polys": polys, "polys_cut": polys_cut}
                self.updateElements(dets)
                return {"list_edges": self.getElement("list_edges")}
            return {}

except ImportError:
    print("GeoPlusExtension is not available (requires shapely package)...")
