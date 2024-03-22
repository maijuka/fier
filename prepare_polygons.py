import codecs
import re
import io
import numpy
import scipy.spatial
from shapely.geometry import polygon, Point, Polygon
from shapely.geometry.linestring import LineString
from shapely.geometry.multilinestring import MultiLineString

try:
    from csv_reader import getFp
except ModuleNotFoundError:
    from .csv_reader import getFp

import pdb

IN_PACK = True
if __name__ == "__main__":
    IN_PACK = False

BUFFER = 10**-6
TIK = "types_index"
# TYPES legend
# -----------
# outer: cut from an edge from the frame
# inner: cut from edges not from the frame
# org: from the voronoi tessalation
# frame: frame edge, external bounding box
# cut: cut edge
# isolated: part of a cell unconnected to surrounding (all edges cut)
# projected: indirectly joining existing end points
# fill: directly joining two existing end points
# residue: side of a cut edge
TYPES_MAP = {"outer": -10, "inner": -11, "org": 0, "frame": -2, "cut": 1, "isolated": 10, "projected": 20, "fill": 21, "residue": 22}
OUTER_MAP = {True: "outer", False: "inner"}

LOG_ON = False

# PLOTTING
#####################
if IN_PACK:
    # OFF
    def plot_nodes(coords, nids=None, marker="+", color="k"):
        pass

    def plot_edges(eids, list_edges, lbl=None, linestyle="-", marker="", color="k", linewidth=1, kedge="edge"):
        pass

    def plot_filled(eids, list_edges, lbl=None, linestyle="-", marker="", color="k", linewidth=1, kedge="edge"):
        pass

    def plot_edges_colordered(eids, list_edges, lbl=None, linestyle="-", marker="o", linewidth=1, ccolor=None, simple=False, kedge="edge"):
        pass

    def plot_esets_colordered(esets, list_edges, lbls=None, linestyle="-", marker="o", linewidth=1, simple=False, kedge="edge"):
        pass

    def plot_show(block=False):
        pass

else:
    # ON
    import matplotlib.pyplot as plt
    import matplotlib.cm as mcm

    def plot_nodes(coords, nids=None, marker="+", color="k"):
        x, y = zip(*coords)
        plt.plot(x, y, marker=marker, color=color, linestyle="")
        if nids is not None:
            for ii, nid in enumerate(nids):
                plt.text(x[ii], y[ii], nid, color=color)

    def plot_edges(eids, list_edges, lbl=None, linestyle="-", marker="", color="k", linewidth=1, kedge="edge"):
        for pi, pp in enumerate(eids):
            edge = get_ordered_edge(pp, list_edges, kedge=kedge)
            x, y = zip(*edge)
            plt.plot(x, y, marker=marker, color=color, linewidth=linewidth, linestyle=linestyle)
            if lbl == "IDS":
                plt.text((x[0]+x[1])/2., (y[0]+y[1])/2., "%d:%d" % (pi, pp), color=color)
            elif lbl is not None and pi == 0:
                plt.text(x[0], y[0], lbl)

    def plot_filled(eids, list_edges, lbl=None, linestyle="-", marker="", color="k", linewidth=1, kedge="edge"):
        fe = get_ordered_edge(eids[0], list_edges, kedge=kedge)
        points = [fe[0], fe[1]]+[get_ordered_edge(pp, list_edges, kedge=kedge)[1] for pp in eids[1:]]
        xs, ys = zip(*points)
        plt.fill(xs, ys, color=color, linewidth=linewidth, linestyle=linestyle)
        if lbl is not None and pi == 0:
            plt.text(xs[0], ys[0], lbl)

    def plot_edges_colordered(eids, list_edges, lbl=None, linestyle="-", marker="o", linewidth=1, ccolor=None, simple=False, kedge="edge"):
        if ccolor is None:
            color = None
            cmap = mcm.get_cmap('rainbow')
        else:
            color = ccolor
        lp = len(eids)-1.
        for pi, pp in enumerate(eids):
            if ccolor is None:
                color = cmap(pi/lp)
            edge = get_ordered_edge(pp, list_edges, kedge=kedge)
            x, y = zip(*edge)
            plt.plot(x, y, marker="", color=color, linewidth=linewidth, linestyle=linestyle)
            if not simple:
                plt.plot(x[0], y[0], marker=marker, color=color)
                plt.text((x[0]+2*x[1])/3., (y[0]+2*y[1])/3., "%d:%d" % (pi, pp))
            if lbl is not None and pi == 0:
                plt.text(x[0], y[0], lbl)

    def plot_esets_colordered(esets, list_edges, lbls=None, linestyle="-", marker="o", linewidth=1, simple=False, kedge="edge"):
        cmap = mcm.get_cmap('rainbow')
        lp = numpy.maximum(1., len(esets)-1.)
        for pi, pp in enumerate(esets):
            color = cmap(pi/lp)
            lbl = None
            if lbls is not None:
                lbl = lbls[pi]
            plot_edges_colordered(pp, list_edges, lbl, linestyle, marker, linewidth, color, simple, kedge=kedge)

    def plot_show(block=False):
        plt.show(block)

#######################################################

# READING
#####################
# parameters = {
#     "NAMES_CID": 0,
#     "LAT_CID": 2,
#     "LNG_CID": 1,
#     "SEP": ",",
#     "ID_INT": False
#     }

if True:
    # if IN_PACK:
    #     def read_coords_csv(filename, csv_params={}, unknown_string=None):
    #         return [], []
    # else:
    import csv

    LATITUDE = ('lat', 'latitude', 'Lat', 'Latitude', 'lats', 'latitudes', 'Lats', 'Latitudes')
    LONGITUDE = ('long', 'longitude', 'Long', 'Longitude', 'longs', 'longitudes', 'Longs', 'Longitudes')
    IDENTIFIERS = ('id', 'identifier', 'Id', 'Identifier', 'ids', 'identifiers', 'Ids', 'Identifiers', 'ID', 'IDS')
    COND_COL = ('cond_var', 'cond_col', 'cond_time')

    ENABLED_ROWS = ('enabled_row', 'enabled_rows')
    ENABLED_COLS = ('enabled_col', 'enabled_cols')
    GROUPS_COLS = ('groups_col', 'groups_cols')

    COLVAR = ['cid', 'CID', 'cids', 'CIDS', 'variable', 'Variable', 'variables', 'Variables']
    COLVAL = ['value', 'Value', 'values', 'Values']

    def read_coords_csv(filename, csv_params={}, unknown_string=None):
        f, fcl = getFp(filename)
        if f is not None:
            try:
                dialect = csv.Sniffer().sniff(f.read(2048))
            except Exception:
                dialect = "excel"
            f.seek(0)
            #header = csv.Sniffer().has_header(f.read(2048))
            # f.seek(0)
            csvreader = csv.reader(f, dialect=dialect, **csv_params)
            # Try to read headers
            head = next(csvreader)

            cpos = {}
            for i, h in enumerate(head):
                for clbls in [LATITUDE, LONGITUDE, IDENTIFIERS]:
                    if h in clbls:
                        cpos[clbls[0]] = i

            if not (LATITUDE[0] in cpos and LONGITUDE[0] in cpos):
                return None, None
            cmax = max(cpos.values())
            coords = []
            if IDENTIFIERS[0] in cpos:
                rnames = []
            else:
                rnames = None

            for row in csvreader:
                if re.match("\s*#", row[0]) or row[0] in ENABLED_ROWS+ENABLED_COLS+GROUPS_COLS:
                    continue
                if len(row) < cmax+1:
                    raise ValueError('number of columns does not match (is ' +
                                     str(len(row))+', should be at least' +
                                     str(cmax+1)+')')

                cc = (float(row[cpos[LONGITUDE[0]]].strip()), float(row[cpos[LATITUDE[0]]].strip()))
                # WARNING! filter a rectangle by coordinates
                # if cc[0] > -20 or cc[1] > -40:
                #     continue
                coords.append(cc)
                # coords.append((float(row[cpos[LATITUDE[0]]].strip()), float(row[cpos[LONGITUDE[0]]].strip())))
                if rnames is not None:
                    tmp = row[cpos[IDENTIFIERS[0]]].strip()
                    if tmp != type(tmp)(unknown_string):
                        rnames.append(tmp)
                    else:
                        rnames.append(None)
        if fcl:
            f.close()
        return coords, rnames
#######################################################

######################################
# PREPARING POLYGON FROM COORDS
######################################


def getInterYV(x1, y1, x2, y2, xv):
    if x2 == x1:
        return numpy.nan
    # get intersection of vertical line x=xv, with line going through (x1, y1) and (x2, y2)
    return y1+(xv-x1)*(y2-y1)/(x2-x1)


def getInterXH(x1, y1, x2, y2, yh):
    if y2 == y1:
        return numpy.nan
    # get intersection of horizontal line y=yh, with line going through (x1, y1) and (x2, y2)
    return x1+(yh-y1)*(x2-x1)/(y2-y1)


def add_points(polygon, points):
    new_coords = list(polygon.boundary.coords)
    changed = False
    for point in points:
        if point not in new_coords:
            i = 1
            while i < len(new_coords):
                if Point(point).within(LineString((new_coords[i-1], new_coords[i])).buffer(BUFFER)):
                    changed = True
                    new_coords.insert(i, point)
                    i = len(new_coords)+1
                else:
                    i += 1
            # if i == len(new_coords):
            #     raise Exception("No insertion point found!")
    if changed:
        # print("Changed", new_coords)
        return Polygon(new_coords)
    return polygon


FACT_CUT = 0.8


def getCutPoint(PointsMap, edge_data, end_cut, gridHW):
    # cut edge on end_cut
    f = FACT_CUT
    edge = edge_data["edge"]
    nA, nB = edge_data["nodes"]
    midX = (PointsMap[nA][0] + PointsMap[nB][0])/2.
    midY = (PointsMap[nA][1] + PointsMap[nB][1])/2.

    cancel = {"V": False, "H": False}
    (xv, yv) = (edge[end_cut][0], edge[end_cut][1])
    (xh, yh) = (edge[end_cut][0], edge[end_cut][1])
    if edge[end_cut][0] < edge[1-end_cut][0]:  # end to cut is to the left
        xv = midX - f*getRectW(midX, midY, gridHW)
        yv = getInterYV(edge[0][0], edge[0][1], edge[1][0], edge[1][1], xv)
    elif edge[end_cut][0] > edge[1-end_cut][0]:  # end to cut is to the right
        xv = midX + f*getRectW(midX, midY, gridHW)
        yv = getInterYV(edge[0][0], edge[0][1], edge[1][0], edge[1][1], xv)
    else:  # vertical edge
        cancel["V"] = True
        # print("No horizontal intersection")
    if edge[end_cut][1] < edge[1-end_cut][1]:  # end to cut is at the bottom
        yh = midY - f*getRectH(midX, midY, gridHW)
        xh = getInterXH(edge[0][0], edge[0][1], edge[1][0], edge[1][1], yh)
    elif edge[end_cut][1] > edge[1-end_cut][1]:  # end to cut is at the top
        yh = midY + f*getRectH(midX, midY, gridHW)
        xh = getInterXH(edge[0][0], edge[0][1], edge[1][0], edge[1][1], yh)
    else:  # horizontal edge
        cancel["H"] = True
        # print("No vertical intersection")

    if xv < numpy.minimum(edge[0][0], edge[1][0]) or xv > numpy.maximum(edge[0][0], edge[1][0]) or yv < numpy.minimum(edge[0][1], edge[1][1]) or yv > numpy.maximum(edge[0][1], edge[1][1]):
        (xv, yv) = (edge[end_cut][0], edge[end_cut][1])
        cancel["V"] = True
    if xh < numpy.minimum(edge[0][0], edge[1][0]) or xh > numpy.maximum(edge[0][0], edge[1][0]) or yh < numpy.minimum(edge[0][1], edge[1][1]) or yh > numpy.maximum(edge[0][1], edge[1][1]):
        (xh, yh) = (edge[end_cut][0], edge[end_cut][1])
        cancel["H"] = True
    if cancel["H"] and cancel["V"]:
        ### print("Cancel cut", nA, nB)
        if edge_data["n_closer"] == end_cut:
            (xh, yh) = (edge[edge_data["n_closer"]][0], edge[edge_data["n_closer"]][1])
        else:
            (xh, yh) = ((2*edge[end_cut][0]+1*edge[1-end_cut][0])/3., (2*edge[end_cut][1]+1*edge[1-end_cut][1])/3.)

    dh = (xh-edge[1-end_cut][0])**2 + (yh-edge[1-end_cut][1])**2
    dv = (xv-edge[1-end_cut][0])**2 + (yv-edge[1-end_cut][1])**2
    dorg = (edge[end_cut][0]-edge[1-end_cut][0])**2 + (edge[end_cut][1]-edge[1-end_cut][1])**2
    select_cut = None
    if dh < dorg:
        if dv < dh:
            select_cut = (xv, yv)
        else:
            select_cut = (xh, yh)
    elif dv < dorg:
        select_cut = (xv, yv)
    return select_cut


def getRectW(x, y, gridHW):
    return gridHW["W"]


def getRectH(x, y, gridHW):
    # get height of standard edge for coordinates (x, y)
    return gridHW["H"]


def getPolyRect(x, y, gridHW, f):
    # get standard rectangle for coordinates (x, y)
    w = getRectW(x, y, gridHW)
    h = getRectH(x, y, gridHW)
    return [(x+fx*f*w, y+fy*f*h) for (fx, fy) in [(1, 1), (1, -1), (-1, -1), (-1, 1), (1, 1)]]

# tool function
####################


def ortho_prj(x1, y1, x2, y2, xn, yn):
    xA, yA = (x2-x1, y2-y1)
    xB, yB = (xn-x1, yn-y1)
    lA = numpy.sqrt(xA**2+yA**2)
    lB = numpy.sqrt(xB**2+yB**2)
    cosA = (xA*xB+yA*yB)/(lA*lB)
    return (x1 + cosA*lB*xA/lA, y1 + cosA*lB*yA/lA)


def flatten_poly(poly):
    # turn bunch of edges into polygon as ordered sequence of edges
    gpoly = {}
    map_pis = {}
    for pi, (p1, p2) in enumerate(poly):
        map_pis[(p1, p2)] = pi+1
        map_pis[(p2, p1)] = -(pi+1)
        if p1 not in gpoly:
            gpoly[p1] = [p2]
        else:
            gpoly[p1].append(p2)
        if p2 not in gpoly:
            gpoly[p2] = [p1]
        else:
            gpoly[p2].append(p1)

    ks = gpoly.keys()
    for k in ks:
        gpoly[k].sort()

    fps = []
    fpis = []
    prev_prev = None
    prev_reuse = False
    while len(gpoly) > 0:
        if prev_prev is not None and len(gpoly.get(prev_prev, [])) > 0:
            prev_reuse = True
            prev_node = prev_prev
            fps[-1].reverse()
            fpis[-1].reverse()
        else:
            prev_node = sorted(gpoly.keys(), key=lambda x: (-len(gpoly[x]), x))[0]
            prev_reuse = False
            prev_prev = prev_node
            fps.append([prev_node])
            fpis.append([])

        while len(gpoly[prev_node]) > 0:
            next_node = gpoly[prev_node].pop(0)
            gpoly[next_node].remove(prev_node)
            fps[-1].append(next_node)
            fpis[-1].append(map_pis[(prev_node, next_node)])
            prev_node = next_node

        ks = list(gpoly.keys())
        for k in ks:
            if len(gpoly[k]) == 0:
                del gpoly[k]
        # if len(gpoly) > 0:
        #     print("More than one polygon")
        #     pdb.set_trace()
        #     print(gpoly)
        if (len(fps[-1]) == 1) or (len(fps[-1]) == 2 and fps[-1][0] == fps[-1][1]):
            fps.pop()
            fpis.pop()
    return fps, fpis


def round_poly(poly, nb_dec=None):
    # drop successive equal dots in polygons
    if nb_dec is None:
        return poly
    return [(numpy.around(p[0], nb_dec), numpy.around(p[1], nb_dec)) for p in poly]


def dedup_poly(poly):
    # drop successive equal dots in polygons
    i = 1
    while i < len(poly):
        if poly[i] == poly[i-1]:
            poly.pop(i)
        else:
            i += 1
    return poly


def decomplex_poly(poly):
    # break up loops in polygon to separate polygons
    copy_poly = list(poly)
    i = 1
    subpolys = []
    while i < len(poly):
        if poly[i] in poly[:i]:
            j = poly[:i].index(poly[i])
            tmp_poly = []
            for k in range(i, j, -1):
                tmp_poly.append(poly.pop(k))
            subpolys.append([poly[j]]+tmp_poly[::-1])
            i = j
        else:
            i += 1
    if len(poly) > 1:
        subpolys.append(poly)
    return subpolys


def clean_poly(poly):
    return decomplex_poly(dedup_poly(poly))


def smoothen_polys(polys, fact=1):
    if fact <= 1:
        return polys
    return [smoothen_poly(poly, fact) for poly in polys]


def smoothen_poly(poly, fact=1):
    if fact <= 1 or len(poly) / fact < 10:
        return poly
    return poly[:-(fact-1):fact]+[poly[-1]]


###########################################################################
def get_ordered_edge_flatten(seid, list_edges, sign=None, kedge="edge"):
    if seid == 0:
        return [0, 1, 1, 0]
    eA, eB = get_ordered_edge(seid, list_edges, sign, kedge)
    return [eA[0], eA[1], eB[0], eB[1]]
    # if (sign is None and seid > 0) or (sign is not None and sign > 0):
    #     return [list_edges[abs(seid)][kedge][0][0], list_edges[abs(seid)][kedge][0][1], list_edges[abs(seid)][kedge][1][0], list_edges[abs(seid)][kedge][1][1]]
    # else:
    #     return [list_edges[abs(seid)][kedge][1][0], list_edges[abs(seid)][kedge][1][1], list_edges[abs(seid)][kedge][0][0], list_edges[abs(seid)][kedge][0][1]]


def get_ordered_edge(seid, list_edges, sign=None, kedge="edge"):
    if (sign is None and seid > 0) or (sign is not None and sign > 0):
        return (list_edges[abs(seid)][kedge][0], list_edges[abs(seid)][kedge][1])
    else:
        return (list_edges[abs(seid)][kedge][1], list_edges[abs(seid)][kedge][0])


def get_cut_edge(list_edges, seid, node=None, kedge="edge"):
    if "cut_eid" in list_edges[abs(seid)]:
        ceid = list_edges[abs(seid)]["cut_eid"]
        # if ceid < 0 or seid < 0: pdb.set_trace()
        edge = get_ordered_edge(ceid, list_edges, kedge=kedge)
        if seid < 0:
            ceid *= -1
            edge = [edge[1], edge[0]]
    else:
        ceid = seid
        edge = get_ordered_edge(seid, list_edges, kedge=kedge)

    if node is not None and node not in list_edges[abs(ceid)]["nodes"]:
        list_edges[abs(ceid)]["nodes"].append(node)
    return ceid, edge


def order_edge(prev_point, next_point=None):
    if next_point is None:
        prev_point, next_point = prev_point
    if prev_point < next_point:
        sign = 1
        new_edge = (prev_point, next_point)
    else:
        sign = -1
        new_edge = (next_point, prev_point)
    return new_edge, sign


def update_types(list_edges, neid, type_edge=None, from_edge=False):
    if type_edge is not None:
        tte = type_edge
        if not type(type_edge) in [list, set]:
            tte = [type_edge]
        if not from_edge:
            list_edges[neid]["types"].update(tte)
        if list_edges[0] is not None and TIK in list_edges[0]:
            for te in tte:
                assert(te in TYPES_MAP)
                if te not in list_edges[0][TIK]:
                    list_edges[0][TIK][te] = set()
                list_edges[0][TIK][te].add(neid)


def create_new_edge(prev_point,  next_point, node, map_edges, list_edges, type_edge=None, add_data={}):
    new_edge, sign = order_edge(prev_point, next_point)
    if new_edge in map_edges:
        # if not add_force: pdb.set_trace() ### something strange
        neid = map_edges[new_edge]
    else:
        neid = len(list_edges)
        map_edges[new_edge] = neid
        list_edges.append({"edge": new_edge, "nodes": [], "pos": [], "types": set()})
    if node is not None and node not in list_edges[neid]["nodes"]:
        list_edges[neid]["nodes"].append(node)
    list_edges[neid].update(add_data)
    update_types(list_edges, neid, type_edge)
    return sign*neid


def query_edge(prev_point,  next_point, map_edges):
    new_edge, sign = order_edge(prev_point, next_point)
    if new_edge in map_edges:
        # if not add_force: pdb.set_trace() ### something strange
        neid = map_edges[new_edge]
        return sign*neid
    return None

##############################
# ACTUAL VORONOI COMPUTATION
##############################


def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """
    # https://stackoverflow.com/questions/20515554/colorize-voronoi-diagram/20678647#20678647

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()*2

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all([v >= 0 for v in vertices]):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1]  # tangent
            t /= numpy.linalg.norm(t)
            n = numpy.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = numpy.sign(numpy.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = numpy.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = numpy.arctan2(vs[:, 1] - c[1], vs[:, 0] - c[0])
        new_region = numpy.array(new_region)[numpy.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, numpy.asarray(new_vertices)


def VoronoiPolygons(points, radius=None):
    vor = scipy.spatial.Voronoi(points)
    last_finite_vertice = len(vor.vertices)
    new_regions, new_vertices = voronoi_finite_polygons_2d(vor, radius)
    return new_regions, new_vertices, last_finite_vertice


def clip_bbox(points, regions, vertices, bbox):

    bbox_poly = Polygon([(bbox[0], bbox[1]), (bbox[2], bbox[1]), (bbox[2], bbox[3]), (bbox[0], bbox[3])])

    points_in_bbox = numpy.array([bbox_poly.contains(Point(p)) for p in vertices], dtype=bool)
    regions_in_box = numpy.array([all(points_in_bbox[r]) for r in regions], dtype=bool)

    new_vertices = vertices.tolist()
    last_existing_vertice = len(new_vertices)
    points_map = dict([(tuple(p), pi) for pi, p in enumerate(new_vertices)])

    ris_cut = numpy.where(~regions_in_box)[0]
    for ri in ris_cut:
        tmp = [tuple(vertices[r, :]) for r in regions[ri]]
        P = Polygon([tuple(vertices[r, :]) for r in regions[ri]])
        X = bbox_poly.intersection(P)
        new_region = []
        for x in X.exterior.coords:
            if x not in points_map:
                points_map[x] = len(new_vertices)
                new_vertices.append(x)
            new_region.append(points_map[x])
        regions[ri] = new_region
    vertices = numpy.array(new_vertices)
    return regions, vertices, last_existing_vertice


def clip_squares(points, regions, vertices, sw, sh=None):
    if sh is None:
        sh = sw
    new_vertices = vertices.tolist()
    last_existing_vertice = len(new_vertices)
    points_map = dict([(tuple(p), pi) for pi, p in enumerate(new_vertices)])

    nb_regions = len(regions)
    for ri in range(nb_regions):
        point = points[ri, :]
        region = regions[ri]
        bbox_poly = Polygon([(point[0]-sw/2, point[1]-sh/2), (point[0]+sw/2, point[1]-sh/2), (point[0]+sw/2, point[1]+sh/2), (point[0]-sw/2, point[1]+sh/2)])

        if not all([bbox_poly.contains(Point(vertices[p, :])) for p in region]):
            P = Polygon([tuple(vertices[r, :]) for r in region])
            X = bbox_poly.intersection(P)
            new_region = []
            for x in X.exterior.coords:
                if x not in points_map:
                    points_map[x] = len(new_vertices)
                    new_vertices.append(x)
                new_region.append(points_map[x])
        regions[ri] = new_region
    vertices = numpy.array(new_vertices)
    return regions, vertices, last_existing_vertice


def get_ebbox(bbox, bbox_ext=None, marg_ratio=100):
    # bbox = (numpy.min(vor.vertices[:,0]), numpy.min(vor.vertices[:,1]),
    #         numpy.max(vor.vertices[:,0]), numpy.max(vor.vertices[:,1]))
    dbx, dby = (bbox[2] - bbox[0])/marg_ratio, (bbox[3] - bbox[1])/marg_ratio
    ebbox = [bbox[0]-dbx, bbox[1]-dby, bbox[2]+dbx, bbox[3]+dby]
    if bbox_ext is not None:
        assert(bbox[0] >= bbox_ext[0])
        assert(bbox[1] >= bbox_ext[1])
        assert(bbox[2] <= bbox_ext[2])
        assert(bbox[3] <= bbox_ext[3])
        ebbox = [numpy.maximum(ebbox[0], bbox_ext[0]), numpy.maximum(ebbox[1], bbox_ext[1]),
                 numpy.minimum(ebbox[2], bbox_ext[2]), numpy.minimum(ebbox[3], bbox_ext[3])]
    return ebbox


def borders_exterior(poly_points, map_edges, list_edges, polys, key_nodes="nodes"):
    O = list(polygon.orient(Polygon(poly_points)).exterior.coords)
    oeids = [query_edge(O[i-1], O[i], map_edges) for i in range(1, len(O))]
    if len(list_edges[abs(oeids[0])].get(key_nodes, [])) == 1:
        peids = polys[list_edges[abs(oeids[0])][key_nodes][0]]
        E = list(polygon.orient(Polygon(get_poly_coords(list_edges, peids))).exterior.coords)
        i0, i1 = E.index(O[0]), E.index(O[1])
        if abs(i1-i0) != 1:
            if i0 == (len(E)-2):  # i0 indexes the one to last point in E (last point is repeat from first)
                i0 = -1
            elif i1 == (len(E)-2):  # i1 indexes the one to last point in E (last point is repeat from first)
                i1 = -1
    else:
        raise Exception("Can't find the polygon for exterior edge %s!" % oeids[0])
    return i1 > i0


def borders_exterior_check(poly_points, map_edges, list_edges, polys, key_nodes="nodes"):
    O = list(polygon.orient(Polygon(poly_points)).exterior.coords)
    oeids = [query_edge(O[i-1], O[i], map_edges) for i in range(1, len(O))]
    info_edges = []
    for opos, oeid in enumerate(oeids):
        if len(list_edges[abs(oeid)].get(key_nodes, [])) == 1:
            peids = polys[list_edges[abs(oeid)][key_nodes][0]]
            E = list(polygon.orient(Polygon(get_poly_coords(list_edges, peids))).exterior.coords)
            i0, i1 = E.index(O[opos]), E.index(O[opos+1])
            if abs(i1-i0) != 1:
                if i0 == (len(E)-2):  # i0 indexes the one to last point in E (last point is repeat from first)
                    i0 = -1
                elif i1 == (len(E)-2):  # i1 indexes the one to last point in E (last point is repeat from first)
                    i1 = -1
            info_edges.append((i1-i0, i0, i1, list_edges[abs(oeid)][key_nodes]))

        else:
            raise Exception("Can't find the polygon for exterior edge %s!" % oeid)
    if not all([e[0] == info_edges[0][0] for e in info_edges[1:]]):
        raise Exception("Orientation of all edges don't agree! %s" % info_edges)
    return info_edges[0][0] > 0


def get_poly_coords(list_edges, peids, kedge="edge"):
    return [get_ordered_edge(peid, list_edges, kedge=kedge)[0] for peid in peids]


def compute_voronoi(PointsMap, bbox_ext=None):
    # MAIN FUNCTION for computing polygons from coordinates
    # Compute the voronoi diagram
    #################################
    # coordinates of the points
    lbls, coords = zip(*PointsMap.items())

    points = numpy.array(coords)
    vor = scipy.spatial.Voronoi(points)

    bbox = (numpy.min(points[:, 0]), numpy.min(points[:, 1]), numpy.max(points[:, 0]), numpy.max(points[:, 1]))
    ebbox = get_ebbox(bbox, bbox_ext)

    radius = numpy.sqrt((ebbox[2]-ebbox[0])**2+(ebbox[3]-ebbox[1])**2)*2
    regions, vertices, lfinite = VoronoiPolygons(coords, radius)
    regions, vertices, lexist = clip_bbox(points, regions, vertices, ebbox)

    # cmap = mcm.get_cmap('rainbow')
    # # nb_regions = len(regions)
    # nb_regions = 10
    # fig, ax = plt.subplots()
    # plt.plot([ebbox[0], ebbox[2], ebbox[2], ebbox[0], ebbox[0]], [ebbox[1], ebbox[1], ebbox[3], ebbox[3], ebbox[1]], "-k", linewidth=.5)
    # # ax.set_aspect(1.0)
    # for ri, region in enumerate(regions[:nb_regions]):
    #     color = cmap(ri/nb_regions)
    #     plt.plot(vertices[region+region[:1],0], vertices[region+region[:1],1], marker="", color=color, linestyle=":") #, linewidth=linewidth)
    #     # plt.plot(vertices[lfinite:,0], vertices[lfinite:,1], marker="o", color='k', linestyle='none') #, linewidth=linewidth)
    #     plt.plot(coords[ri][0], coords[ri][1], marker="s", color=color)
    #     plt.text(coords[ri][0], coords[ri][1], "%d:%s" % (ri, lbls[ri]))
    # plt.show()

    map_edges = {}
    list_edges = [{"edge": None, TIK: {}, "nb_nodes": len(points)}]  # so no real edge is indexed with 0, which has no sign and does not allow marking direction
    polys = {}
    for ri, region in enumerate(regions):
        pids = []
        for i in range(len(region)):
            neid = create_new_edge((vertices[region[i-1], 0], vertices[region[i-1], 1]),
                                   (vertices[region[i], 0], vertices[region[i], 1]),
                                   ri, map_edges, list_edges, type_edge="org")
            pids.append(neid)
        polys[ri] = pids
        for pi, eid in enumerate(pids):
            list_edges[abs(eid)]["pos"].append((pi+1)*numpy.sign(eid))

    # exactly two adjacent nodes, except for border edges that have only one
    ## assert all([(len(x["nodes"])==2 or len(x["nodes"])==1) for (k,x) in map_edges.items()])
    return map_edges, list_edges, polys, ebbox


ANGLE_TOL = 0.05


def compute_edges_distances(PointsMap, list_edges):
    horz_edges = []
    vert_edges = []
    for eid in range(1, len(list_edges)):
        dt = list_edges[eid]
        if len(dt["nodes"]) == 1:  # border edge
            update_types(list_edges, eid, "frame")
            dt["far"] = -1
        else:  # other edge
            edge = dt["edge"]
            midx = numpy.abs(PointsMap[dt["nodes"][0]][0] + PointsMap[dt["nodes"][1]][0])/2.
            midy = numpy.abs(PointsMap[dt["nodes"][0]][1] + PointsMap[dt["nodes"][1]][1])/2.
            dx = numpy.abs(PointsMap[dt["nodes"][0]][0] - PointsMap[dt["nodes"][1]][0])
            dy = numpy.abs(PointsMap[dt["nodes"][0]][1] - PointsMap[dt["nodes"][1]][1])
            closer = 0
            if (midx-edge[1][0])**2+(midy-edge[1][1])**2 < (midx-edge[0][0])**2+(midy-edge[0][1])**2:
                closer = 1
            dt.update({"dx": dx, "dy": dy, "far": 0, "n_closer": closer})
            horz_edges.append(dx)
            vert_edges.append(dy)

    vert_edges = numpy.array(vert_edges)
    horz_edges = numpy.array(horz_edges)
    return vert_edges, horz_edges, list_edges


def compute_grid_size(vert_edges, horz_edges, hgrid_percentile=-1, wgrid_percentile=-1):
    # return {'H': 0.007839320868257893, 'W': 0.0078403303939800773}
    if hgrid_percentile <= 0:
        hgrid_percentile = 50
    if wgrid_percentile <= 0:
        wgrid_percentile = hgrid_percentile

    # gridH = numpy.mean(vert_edges)
    # gridW = numpy.mean(horz_edges)
    gridH = numpy.percentile(vert_edges, hgrid_percentile)
    gridW = numpy.percentile(horz_edges, wgrid_percentile)

    gridHW = {"H": gridH, "W": gridW}
    # print("GRID", gridHW)
    return gridHW


def update_far(PointsMap, list_edges, gridHW):
    for eid in range(1, len(list_edges)):
        if list_edges[eid]["far"] == 0:
            dt = list_edges[eid]
            if dt["dx"] > 1.5*gridHW["W"] or dt["dy"] > 1.2*gridHW["H"]:
                dt["far"] = 1
    return list_edges


def gather_require_cut(PointsMap, map_edges, list_edges, polys):
    require_cut = {}
    polys_cut_info = {}
    isolated_nodes = []
    # go over each node in the voronoi diagram, prepare polygon
    for node, poly in polys.items():
        cut_info = [-int(list_edges[abs(seid)]["far"] != 0) for seid in poly]
        if numpy.sum(cut_info) == 0:  # contains at least one edge far==0
            continue
        elif numpy.prod(cut_info) == 0:  # contains at least one edge far==0
            # collect edges that are not far but adjacent to far, which will need to be cut
            for pi, seid in enumerate(poly):
                if list_edges[abs(seid)]["far"] != 0:
                    if list_edges[abs(poly[pi-1])]["far"] == 0:
                        oei = abs(poly[pi-1])
                        wend = int(numpy.sign(poly[pi-1]) > 0)
                        if oei not in require_cut:
                            require_cut[oei] = []
                        require_cut[oei].append((node, pi, wend))
                        cut_info[pi-1] = 1

                    aft = pi+1
                    if aft == len(poly):
                        aft = 0
                    if list_edges[abs(poly[aft])]["far"] == 0:
                        oei = abs(poly[aft])
                        wend = int(numpy.sign(poly[aft]) < 0)
                        if oei not in require_cut:
                            require_cut[oei] = []
                        require_cut[oei].append((node, pi, wend))
                        cut_info[aft] = 1
            polys_cut_info[node] = cut_info
        else:
            isolated_nodes.append(node)
    return require_cut, polys_cut_info, isolated_nodes


def cut_edges(PointsMap, map_edges, list_edges, require_cut, gridHW):
    # cutting edges
    recut_nodes = {}
    for ei, dt in require_cut.items():
        edge = [list_edges[ei]["edge"][0], list_edges[ei]["edge"][1]]
        end_points = set([d[-1] for d in dt])
        end_cut = []
        for end_point in end_points:
            cut_point = getCutPoint(PointsMap, list_edges[ei], end_point, gridHW)
            if cut_point is not None:
                end_cut.append(end_point)
                edge[end_point] = cut_point
        if len(end_cut) > 0:
            neid = create_new_edge(edge[0], edge[1], None, map_edges, list_edges, type_edge=["cut", OUTER_MAP[False]])
            nodes = set([d[0] for d in dt])
            recut_n = set(list_edges[ei]["nodes"]).difference(nodes)
            for rn in recut_n:
                if rn not in recut_nodes:
                    recut_nodes[rn] = {}
                recut_nodes[rn][ei] = end_cut
        else:
            neid = ei
        list_edges[ei]["cut_eid"] = neid
        list_edges[abs(neid)]["uncut_eid"] = ei
    return recut_nodes


FACT_ALONE = 0.5


def assemble_isolated(PointsMap, map_edges, list_edges, polys, isolated_nodes, gridHW, assembled={}):
    for node in isolated_nodes:
        epoly = getPolyRect(PointsMap[node][0], PointsMap[node][1], gridHW, FACT_ALONE)
        map_to = []
        # does it connects to the outer border?
        outer = any([list_edges[abs(k)]["far"] == -1 for k in polys[node]])
        for pi in range(len(epoly)-1):
            neid = create_new_edge(epoly[pi], epoly[pi+1], node, map_edges, list_edges, type_edge=["isolated", OUTER_MAP[outer]])
            map_to.append(neid)
        assembled[node] = map_to
    return assembled


def clip_edge(edge, bbox, c=False):
    marg = 10**-6
    new_edge = [[edge[0][0], edge[0][1]], [edge[1][0], edge[1][1]]]
    for end_point in [0, 1]:
        # if new_edge[end_point][0] < bbox[0] or new_edge[end_point][1] < bbox[1] or new_edge[end_point][0] > bbox[2] or new_edge[end_point][1] > bbox[3]:
        if new_edge[end_point][0] < bbox[0]:
            new_edge[end_point][1] = getInterYV(new_edge[0][0], new_edge[0][1], new_edge[1][0], new_edge[1][1], bbox[0]+marg)
            new_edge[end_point][0] = bbox[0]+marg
            if numpy.isnan(new_edge[end_point][1]):  # vertical line left of lower bound
                new_edge[1-end_point][0] = bbox[0]+marg
                new_edge[end_point][1] = edge[end_point][1]

        if new_edge[end_point][1] < bbox[1]:
            new_edge[end_point][0] = getInterXH(new_edge[0][0], new_edge[0][1], new_edge[1][0], new_edge[1][1], bbox[1]+marg)
            new_edge[end_point][1] = bbox[1]+marg
            if numpy.isnan(new_edge[end_point][0]):  # horizontal line below lower bound
                new_edge[1-end_point][1] = bbox[1]+marg
                new_edge[end_point][0] = edge[end_point][0]

        if new_edge[end_point][0] > bbox[2]:
            new_edge[end_point][1] = getInterYV(new_edge[0][0], new_edge[0][1], new_edge[1][0], new_edge[1][1], bbox[2]-marg)
            new_edge[end_point][0] = bbox[2]-marg
            if numpy.isnan(new_edge[end_point][1]):  # vertical line right of upper bound
                new_edge[1-end_point][0] = bbox[2]-marg
                new_edge[end_point][1] = edge[end_point][1]

        if new_edge[end_point][1] > bbox[3]:
            new_edge[end_point][0] = getInterXH(new_edge[0][0], new_edge[0][1], new_edge[1][0], new_edge[1][1], bbox[3]-marg)
            new_edge[end_point][1] = bbox[3]-marg
            if numpy.isnan(new_edge[end_point][0]):  # horizontal line below lower bound
                new_edge[1-end_point][1] = bbox[3]-marg
                new_edge[end_point][0] = edge[end_point][0]

    return ((new_edge[0][0], new_edge[0][1]), (new_edge[1][0], new_edge[1][1]))


def assemble_cut(PointsMap, map_edges, list_edges, polys, polys_cut_info, recut_nodes, assembled={}):
    # LOG_ON = False
    # ND = 279 # 3274
    # if ND in polys_cut_info:
    #     LOG_ON = True
    #     polys_cut_info = {ND: polys_cut_info[ND]}

    for node, info in polys_cut_info.items():
        poly = polys[node]
        modified = False
        dropped = 0
        if node in recut_nodes:
            recuts = recut_nodes.pop(node)
            for ii, p in enumerate(poly):
                if abs(p) in recuts:
                    info[ii] = 2

        outer = False
        pids = []
        process_order = []

        if numpy.prod(info) != 0:
            modified = True
            ones_pos = [zid for zid, i in enumerate(info) if i > 0]
            first_zid = ones_pos.pop(0)
            process_order = list(range(first_zid, len(info)))+list(range(first_zid))
            if len(ones_pos) == 0:
                outer = any([list_edges[abs(k)]["far"] == -1 for k in poly])
                seid = poly[first_zid]
                ceid, edge = get_cut_edge(list_edges, seid, node)

                xp, yp = ortho_prj(edge[0][0], edge[0][1], edge[1][0], edge[1][1], PointsMap[node][0], PointsMap[node][1])
                tx, ty = (2.*(PointsMap[node][0]-xp), 2.*(PointsMap[node][1]-yp))
                new_edge = (edge[1][0]+tx, edge[1][1]+ty), (edge[0][0]+tx, edge[0][1]+ty)
                new_edge = clip_edge(new_edge, list_edges[0]["ebbox"])
                neid = create_new_edge(new_edge[0], new_edge[1], node, map_edges, list_edges, type_edge=["projected", OUTER_MAP[outer]])
                pids = [neid]
            elif len(ones_pos) == 1:
                first_seid = poly[first_zid]
                corner = None
                if ones_pos[0] == first_zid + 1:
                    corner = get_ordered_edge(first_seid, list_edges)[1]
                    if "cut_eid" in list_edges[abs(first_seid)]:
                        ccc = get_cut_edge(list_edges, first_seid)
                        if ccc[1][1] != corner:
                            corner = None
                elif first_zid == 0 and ones_pos[0] == len(info)-1:
                    corner = get_ordered_edge(first_seid, list_edges)[0]
                    if "cut_eid" in list_edges[abs(first_seid)]:
                        ccc = get_cut_edge(list_edges, first_seid)
                        if ccc[1][0] != corner:
                            corner = None
                    process_order = [len(info)-1] + list(range(len(info)-1))

                if corner is not None:
                    outer = any([list_edges[abs(k)]["far"] == -1 for k in poly])
                    new_corner = (corner[0] + 2*(PointsMap[node][0]-corner[0]), corner[1] + 2*(PointsMap[node][1]-corner[1]))
                    new_edge = clip_edge((new_corner, corner), list_edges[0]["ebbox"], c=True)
                    pids = [new_edge[0], None]

        else:
            first_zid = 0
            while first_zid < len(info) and info[first_zid] != 0:
                first_zid += 1
            process_order = list(range(first_zid, len(info)))+list(range(first_zid))

        for i in process_order:
            if info[i] > 0:
                seid = poly[i]
                ceid, ord_edge = get_cut_edge(list_edges, seid, node)
                if seid != ceid:
                    modified = True
                elif dropped == 1:
                    dropped = 0
                prev_point = ord_edge[0]
                if len(pids) > 0:
                    if pids[-1] is None:
                        prev_point = pids[0]
                        pids = []
                    else:
                        prev_point = get_ordered_edge(pids[-1], list_edges)[1]
                        # prev_point = list_edges[abs(pids[-1])]["edge"][int(pids[-1] > 0)]
                if prev_point != ord_edge[0]:
                    neid = create_new_edge(prev_point, ord_edge[0], node, map_edges, list_edges, type_edge=["fill", OUTER_MAP[outer]])
                    pids.append(neid)
                    outer = False

                if len(pids) > 0 and ceid == -pids[-1]:
                    pids.pop()
                else:
                    pids.append(ceid)

            elif info[i] == 0:
                if len(pids) > 0:
                    prev_point = get_ordered_edge(pids[-1], list_edges)[1]
                    curr_point = get_ordered_edge(poly[i], list_edges)[0]
                    if prev_point != curr_point:
                        neid = create_new_edge(prev_point, curr_point, node, map_edges, list_edges, type_edge=["fill", OUTER_MAP[False]])
                        pids.append(neid)
                        modified = True
                pids.append(poly[i])
            else:  # info[i] == -1:
                dropped += 1
                if list_edges[abs(poly[i])]["far"] == -1:
                    outer = True

        # make sure it's closed
        prev_point = get_ordered_edge(pids[-1], list_edges)[1]
        first_point = get_ordered_edge(pids[0], list_edges)[0]
        if prev_point != first_point:
            # print("Doesn't close")
            neid = create_new_edge(prev_point, first_point, node, map_edges, list_edges, type_edge=["fill", OUTER_MAP[outer]])
            if neid == -pids[-1]:
                pids.pop()
            else:
                pids.append(neid)
            outer = False

        if dropped > 1:
            modified = True
        if not modified and sorted(pids) != sorted(poly):
            print("MODIFIED?", node, modified, "actually:", sorted(pids) != sorted(poly))
            print(pids)
            print(poly)
            raise ValueError("Cutting modification error")
        if modified:
            assembled[node] = pids

        if LOG_ON:
            plot_edges(poly, list_edges, lbl="IDS", color="b", linestyle="-", marker="x")
            plot_edges(assembled.get(node, []), list_edges, lbl="IDS", color="r", linestyle=":", marker="+")
            plot_show(True)
    return assembled


def make_list_cutdets(poly, eids):
    list_recut = []
    for p, v in enumerate(poly):
        ext = [(p, 0, 0)]
        if abs(v) in eids:
            if 0 in eids[abs(v)]:
                if v > 0:
                    ext.insert(0, (p, -1))
                else:
                    ext.append((p, 1))
            if 1 in eids[abs(v)]:
                if v < 0:
                    ext.insert(0, (p, -1))
                else:
                    ext.append((p, 1))
        list_recut.extend(ext)
    list_recut.insert(0, list_recut.pop())
    return list_recut


def assemble_recut(PointsMap, map_edges, list_edges, polys, recut_nodes, assembled={}):
    global LOG_ON
    # DEBUG
    # ND = 169 #1735 # 3274
    # if ND in recut_nodes:
    #     LOG_ON = True
    #     recut_nodes = {ND: recut_nodes[ND]}
    #     pdb.set_trace()

    for node, eids in recut_nodes.items():
        poly = polys[node]
        map_poly = [[pi] for pi in poly]
        modified = False
        list_recut = make_list_cutdets(poly, eids)

        i = 0
        if LOG_ON:
            print("NODE", node, list_recut)
            print("INIT MAP_POLY", map_poly)
        while i < len(list_recut):
            if list_recut[i][-1] != 0:
                if LOG_ON:
                    print("Check", list_recut[i])
                pos, where_insert = list_recut[i]
                cut_end = (where_insert+1)//2
                seid = poly[pos]
                ceid = list_edges[abs(seid)]["cut_eid"]
                if seid < 0:
                    ceid *= -1
                edge = get_ordered_edge(seid, list_edges)
                cedge = get_ordered_edge(ceid, list_edges)

                if cedge[cut_end] != edge[cut_end]:
                    jump = False

                    if i+1 < len(list_recut) and list_recut[i][-1]*list_recut[i+1][-1] != 0:
                        jump = True
                        if LOG_ON:
                            print("Check other", list_recut[i+1])

                        pos_other, where_insert_other = list_recut[i+1]
                        cut_end_other = (where_insert_other+1)//2
                        seid_other = poly[pos_other]
                        ceid_other = list_edges[abs(seid_other)]["cut_eid"]
                        if seid_other < 0:
                            ceid_other *= -1
                        edge_other = get_ordered_edge(seid_other, list_edges)
                        cedge_other = get_ordered_edge(ceid_other, list_edges)
                        if cedge_other[cut_end_other] != edge_other[cut_end_other]:
                            # corner recut
                            if LOG_ON:
                                print("corner cut", edge, edge_other, "->", cedge, cedge_other)

                            if cut_end == 0 or cut_end_other == 1:
                                raise Exception("Something unexpected during recut!")

                            neid = create_new_edge(cedge[cut_end], cedge_other[cut_end_other], node, map_edges, list_edges, type_edge=["residue", OUTER_MAP[False]])

                            map_poly[pos][-1] = ceid
                            map_poly[pos_other][-1] = ceid_other
                            map_poly[pos_other].insert(0, neid)
                            if LOG_ON:
                                print(pos, pos_other, ceid, neid, ceid_other)
                                print("MAP_POLY", map_poly)
                            modified = True
                            i += 1
                        else:
                            jump = False

                    if not jump:
                        # single recut
                        if LOG_ON:
                            print("single cut", edge, "->", cedge)

                        neid = create_new_edge(edge[cut_end], cedge[cut_end], node, map_edges, list_edges, type_edge=["residue", OUTER_MAP[False]])
                        if cut_end == 0:
                            map_poly[pos][0] = ceid
                            map_poly[pos].insert(0, neid)
                        else:
                            map_poly[pos][-1] = ceid
                            map_poly[pos].append(-1*neid)
                        modified = True
                else:
                    print("cut redacted")
            i += 1
            ####################

        if modified:
            assembled[node] = []
            for pp in map_poly:
                assembled[node].extend(pp)
        if LOG_ON:
            print("MAP_POLY", map_poly)
            plot_edges(poly, list_edges, lbl="IDS", color="b", linestyle="-", marker="x")
            plot_edges(assembled.get(node, []), list_edges, lbl="IDS", color="r", linestyle=":", marker="+")
            plot_show(True)

    LOG_ON = False
    return assembled


def compute_edges_data(PointsMap, hgrid_percentile=-1, wgrid_percentile=-1, after_cut=True, bbox=None):
    map_edges, list_edges, polys, ebbox = compute_voronoi(PointsMap, bbox)

    # plot_esets_colordered(list(polys.values())[:30], list_edges)
    # plt.plot([ebbox[0], ebbox[2], ebbox[2], ebbox[0], ebbox[0]], [ebbox[1], ebbox[1], ebbox[3], ebbox[3], ebbox[1]], ":b")
    # for ii, (k,v) in enumerate(list(polys.items())[:30]):
    #     # plt.subplot(4, 4, ii+1)
    #     # plot_esets_colordered([list(v)], list_edges)
    #     # plt.plot([ebbox[0], ebbox[2], ebbox[2], ebbox[0], ebbox[0]], [ebbox[1], ebbox[1], ebbox[3], ebbox[3], ebbox[1]], ":b")
    #     plt.plot(PointsMap[k][0], PointsMap[k][1], "+k")
    #     plt.text(PointsMap[k][0], PointsMap[k][1], "%d" % k)
    # plt.show()

    list_edges[0]["last_org"] = len(list_edges)
    list_edges[0]["last_cut"] = len(list_edges)
    list_edges[0]["last_isolated"] = len(list_edges)
    list_edges[0]["dst_type"] = "flat"
    list_edges[0]["ebbox"] = ebbox
    if bbox is not None:
        list_edges[0]["bbox"] = bbox

    polys_cut = {}
    if hgrid_percentile > 0 and after_cut:
        list_edges[0]["grid_percentile"] = (hgrid_percentile, wgrid_percentile)
        vert_edges, horz_edges, _ = compute_edges_distances(PointsMap, list_edges)
        # gridHW = {'H': 0.0078393208682577542, 'W': 0.0078393208682577542} #
        gridHW = compute_grid_size(vert_edges, horz_edges, hgrid_percentile, wgrid_percentile)
        update_far(PointsMap, list_edges, gridHW)
        require_cut, polys_cut_info, isolated_nodes = gather_require_cut(PointsMap, map_edges, list_edges, polys)

        recut_nodes = cut_edges(PointsMap, map_edges, list_edges, require_cut, gridHW)
        list_edges[0]["last_cut"] = len(list_edges)

        assemble_isolated(PointsMap, map_edges, list_edges, polys, isolated_nodes, gridHW, polys_cut)
        list_edges[0]["last_isolated"] = len(list_edges)

        assemble_cut(PointsMap, map_edges, list_edges, polys, polys_cut_info, recut_nodes, polys_cut)
        assemble_recut(PointsMap, map_edges, list_edges, polys, recut_nodes, polys_cut)

    for node, poly in polys.items():
        if node in polys_cut:
            pp = polys_cut[node]
        else:
            pp = poly
        for pi, eid in enumerate(pp):
            if "nodes_cut" not in list_edges[abs(eid)]:
                list_edges[abs(eid)]["nodes_cut"] = []
                list_edges[abs(eid)]["pos_cut"] = []
            list_edges[abs(eid)]["nodes_cut"].append(node)
            list_edges[abs(eid)]["pos_cut"].append((pi+1)*numpy.sign(eid))

    return map_edges, list_edges, polys, polys_cut, ebbox

# Compute the edges data when points coordinates are provided, using distances on globe or flat


def prepare_edges_dst(PointsMap, hgrid_percentile=-1, wgrid_percentile=-1, after_cut=True, dst_type="globe"):
    if dst_type == "globe":
        PointsRad = coordsPointsToGlobe(PointsMap)
        deg_bbox = [-180, -90, 180, 90]
        bbox = numpy.radians(deg_bbox)
    else:
        PointsRad = PointsMap
        deg_bbox = None
        bbox = None

    map_edges, list_edges, polys, polys_cut, ebbox = compute_edges_data(PointsRad, hgrid_percentile=hgrid_percentile, wgrid_percentile=wgrid_percentile, after_cut=after_cut, bbox=bbox)
    if dst_type == "globe":
        coordsEdgesFromGlobe(list_edges, clip=deg_bbox)
        list_edges[0]["dst_type"] = "globe"
    list_edges[0]["source"] = "voronoi"
    return map_edges, list_edges, polys, polys_cut, ebbox

# Compute the edges data when polygons coordinates are provided


def prepare_edges_polys(PointsMap, poly_coords):
    lbls, coords = zip(*PointsMap.items())
    points = numpy.array(coords)

    polys = {}
    poly_xcoords, poly_ycoords = poly_coords[0], poly_coords[1]
    assert(len(poly_xcoords) == len(poly_ycoords))
    nb_regions = len(poly_xcoords)

    polygons = [Polygon(zip(poly_xcoords[ri], poly_ycoords[ri])) for ri in range(nb_regions)]

    vor = scipy.spatial.Voronoi(points)
    for nA, nB in vor.ridge_points:
        X = polygons[nA].intersection(polygons[nB])
        if isinstance(X, Polygon):
            take_from = nA if (polygons[nA].area > polygons[nB].area) else nB
            polytmp = polygons[take_from].difference(X)
            if isinstance(polytmp, Polygon):
                polygons[take_from] = polytmp
            else:
                raise Exception("Diff (%d \ (%d inter %d)) not a polygon!" % (take_from, nA, nB))
            X = polygons[nA].intersection(polygons[nB])

        points_inter = set()
        if isinstance(X, MultiLineString):
            for x in X:
                points_inter.update(x.coords)
        elif isinstance(X, Point) or isinstance(X, LineString):
            points_inter.update(X.coords)
        if len(points_inter) > 0:
            for ni in [nA, nB]:
                polygons[ni] = add_points(polygons[ni], points_inter)

    map_edges = {}
    list_edges = [{"edge": None, TIK: {}, "nb_nodes": len(points)}]  # so no real edge is indexed with 0, which has no sign and does not allow marking direction
    for ri in range(nb_regions):
        pids = []
        polygon_xy = polygons[ri].boundary.coords
        if polygon_xy[0] == polygon_xy[-1]:
            bottom = 1
        else:
            bottom = 0
        for i in range(bottom, len(polygon_xy)):
            neid = create_new_edge(polygon_xy[i-1], polygon_xy[i],
                                   ri, map_edges, list_edges, type_edge="org")
            pids.append(neid)
        polys[ri] = pids
        for pi, eid in enumerate(pids):
            list_edges[abs(eid)]["pos"].append((pi+1)*numpy.sign(eid))

    tmp_edges = list(zip(*map_edges.keys()))
    all_x, all_y = zip(*(tmp_edges[0]+tmp_edges[1]))
    bbox = (numpy.min(all_x), numpy.min(all_y), numpy.max(all_x), numpy.max(all_y))
    ebbox = get_ebbox(bbox)
    list_edges[0]["last_org"] = len(list_edges)
    list_edges[0]["last_cut"] = len(list_edges)
    list_edges[0]["last_isolated"] = len(list_edges)
    list_edges[0]["dst_type"] = "globe"
    list_edges[0]["source"] = "polys"
    list_edges[0]["ebbox"] = ebbox
    if bbox is not None:
        list_edges[0]["bbox"] = bbox

    # ## FULL MAP BY TYPES
    # bounds = [1, list_edges[0]["last_org"], list_edges[0]["last_cut"], list_edges[0]["last_isolated"], len(list_edges)]
    # elements = [(bounds[0], bounds[1], "-", "k"), (bounds[1], bounds[2], "-", "g"), (bounds[2], bounds[3], "-", "m"), (bounds[3], bounds[4], "-", "r")]
    # for (bdw, bup, linestyle, color) in elements:
    #     plot_edges(range(bdw, bup), list_edges, linestyle=linestyle, color=color) #, kedge="flat_edge")
    # xs, ys = zip(*PointsMap.values())
    # plt.plot(xs, ys, "bo")
    # for pi, cc in PointsMap.items():
    #     plt.text(cc[0], cc[1], pi)
    # plot_show()

    return map_edges, list_edges, polys, polys, ebbox

# Compute the edges coordinates in a flattened format


def prepare_edges_coords_flatten(list_edges, seids=None, after_cut=True):
    kedge = "edge"
    if list_edges[0].get("dst_type") == "globe" and "flat_edge" in list_edges[-1]:
        kedge = "flat_edge"
    if seids is None:
        if after_cut:
            up_to = len(list_edges)
        else:
            up_to = list_edges[0]["last_org"]
        return [get_ordered_edge_flatten(eid, list_edges, kedge=kedge) for eid in range(0, up_to)]
    else:
        return [get_ordered_edge_flatten(seid, list_edges, kedge=kedge) for seid in seids]


def get_all_grps(list_edges):
    return [None] + list(range(-2, 2))


def get_edge_grp(list_edges, eid=None):
    if eid is None:
        return None
    grp = 0
    if "nodes_cut" in list_edges[eid]:
        grp += 1
    if "far" in list_edges[eid]:
        grp -= 1
        if not "cut_eid" in list_edges[eid]:
            grp -= 1
    return grp


def make_edges_graph_grps(map_edges, list_edges):
    edges_graph = {}
    all_grps = get_all_grps(list_edges)
    grp_to_edges = dict([(s, []) for s in all_grps])
    for edge, eid in map_edges.items():
        grp = get_edge_grp(list_edges, eid)
        grp_to_edges[grp].append(eid)
        for ni, n in enumerate(edge):
            if n not in edges_graph:
                edges_graph[n] = dict([(s, []) for s in all_grps])
            edges_graph[n][grp].append(edge[1-ni])
        if "cut_eid" in list_edges[eid]:
            ceid = list_edges[eid]["cut_eid"]
            cedge = get_ordered_edge(ceid, list_edges)
            grp_sp = get_edge_grp(list_edges, None)
            added = False
            for ni, n in enumerate(edge):
                if cedge[ni] != n:
                    added = True
                    edges_graph[n][grp_sp].append(cedge[ni])
                    if cedge[ni] not in edges_graph:
                        edges_graph[cedge[ni]] = dict([(s, []) for s in all_grps])
                    edges_graph[cedge[ni]][grp_sp].append(n)
            if added:
                grp_to_edges[grp_sp].extend([eid, abs(ceid)])
    return edges_graph, grp_to_edges


def get_ordered_pairs(v, vs):
    return [order_edge(v, vv)[0] for vv in vs]


# Compute the exterior and interior border polygons
def prepare_border_info(list_edges, map_edges, polys, after_cut=True):
    if after_cut and list_edges[0].get("last_org", 0) < list_edges[0].get("last_cut", 0):
        key_nodes = "nodes_cut"
        up_to = len(list_edges)
    else:
        key_nodes = "nodes"
        up_to = list_edges[0].get("last_org", len(list_edges))

    beids = [eid for eid in range(1, up_to) if len(list_edges[eid].get(key_nodes, [])) == 1]
    out = [list_edges[eids]["edge"] for eids in beids]
    fps, fpis = flatten_poly(out)
    border_info = {}
    for ppi, pp in enumerate(fps):
        tt = clean_poly(pp)
        for ttt in tt:
            peids = [query_edge(ttt[i-1], ttt[i], map_edges) for i in range(1, len(ttt))]
            if None in peids:
                raise Exception("Edge(s) not found! %s -> %s" % (ttt, peids))
            ci = -(len(border_info)+1)
            bext = not borders_exterior(ttt, map_edges, list_edges, polys, key_nodes)
            border_info[ci] = {"ci": ci, "cells": [], "polys": [peids], "color": -1, "level": -int(bext)}
    return border_info

# Prepare adjacency graph


def prepare_nodes_graph(list_edges, after_cut=True):
    nodes_graph = {}
    if after_cut and list_edges[0].get("last_org", 0) < list_edges[0].get("last_cut", 0):
        key_nodes = "nodes_cut"
        up_to = len(list_edges)
    else:
        key_nodes = "nodes"
        up_to = list_edges[0].get("last_org", len(list_edges))

    for eid in range(1, up_to):
        if key_nodes not in list_edges[eid]:
            continue
        nodes = list_edges[eid][key_nodes]
        for node in nodes:
            if node not in nodes_graph:
                nodes_graph[node] = {}
            for nodex in nodes:
                if node != nodex:
                    nodes_graph[node][nodex] = eid
    return nodes_graph

# Prepare list of adjacent nodes


def prepare_node_pairs(list_edges, max_nid=None, after_cut=True):
    node_pairs = []
    if after_cut and (list_edges[0].get("last_org", 0) < list_edges[0].get("last_cut", 0)):
        key_nodes = "nodes_cut"
        up_to = len(list_edges)
    else:
        key_nodes = "nodes"
        up_to = list_edges[0].get("last_org", len(list_edges))

    for eid in range(1, up_to):
        if key_nodes not in list_edges[eid]:
            continue
        nodes = list_edges[eid][key_nodes]
        if len(nodes) == 1:
            if max_nid is None or nodes[0] < max_nid:
                node_pairs.append([eid, nodes[0], -1])
        elif len(nodes) == 2:
            if max_nid is None or (nodes[0] < max_nid and nodes[1] < max_nid):
                node_pairs.append([eid, nodes[0], nodes[1]])
            elif nodes[0] < max_nid:
                node_pairs.append([eid, nodes[0], -2])
            elif nodes[1] < max_nid:
                node_pairs.append([eid, nodes[1], -2])
    return numpy.array(node_pairs, dtype=int)


def compute_cc_polys(list_edges, nodes_sel, polys):
    edge_counts = {}
    for node in nodes_sel:
        for eid in polys.get(node, []):
            edge_counts[abs(eid)] = edge_counts.get(abs(eid), 0) + 1
    borders_eids = [eid for eid, c in edge_counts.items() if c == 1]
    borders_edges = [list_edges[eid]["edge"] for eid in borders_eids]
    fps, fpis = flatten_poly(borders_edges)
    pps = [[numpy.sign(ee)*borders_eids[abs(ee)-1] for ee in pp] for pp in fpis]
    return pps


def ss_splength(source, nodes_graph, nodes_colors=None):
    seen = {}                  # level (number of hops) when seen in BFS
    level = 0                  # the current level
    nextlevel = set([source])  # dict of nodes to check at next level
    while nextlevel:
        thislevel = nextlevel  # advance to next level
        nextlevel = set()         # and start a new list (fringe)
        for v in thislevel:
            if v not in seen:
                seen[v] = level  # set the level of vertex v
                nextlevel.update([c for c in nodes_graph[v].keys() if nodes_colors is None or (c >= 0 and nodes_colors[v] == nodes_colors[c])])  # add neighbors of v
        level = level+1
    return seen  # return all path lengths as dictionary


def connected_components(nodes_graph, nodes_colors, color):
    seen = {}
    pool = numpy.where(nodes_colors == color)[0]
    for v in pool:
        if v not in seen:
            c = ss_splength(v, nodes_graph, nodes_colors)
            yield list(c)
            seen.update(c)

# Collects the polygons (eids)


def prepare_areas_polys(polys, polys_cut, after_cut=True):
    if after_cut:
        pp = dict([(p, polys_cut.get(p) or polys[p]) for p in polys])
    else:
        pp = polys
    return pp

# Prepare the polygons for plotting map with given nodes highlighted


def prepare_areas_data(nodes_colors, list_edges, polys, border_info, nodes_graph):
    # nodes_colors_bckg: contains ids of parts all nodes belong to, i.e. support part, variable value, -2 means outside support, bckg cell
    nodes_colors_bckg = -2*numpy.ones(len(polys), dtype=int)
    if type(nodes_colors) is list:
        nodes_colors_bckg[nodes_colors] = 0
    elif type(nodes_colors) is dict:
        for i, v in nodes_colors.items():
            nodes_colors_bckg[i] = v
    else:
        nodes_colors_bckg[:nodes_colors.shape[0]] = nodes_colors

    # compute color connected components
    ccs_data = {}
    for si in set(nodes_colors_bckg):
        for cc in connected_components(nodes_graph, nodes_colors_bckg, si):
            ci = len(ccs_data)
            cc_polys = compute_cc_polys(list_edges, cc, polys)
            ccs_data[ci] = {"ci": ci, "nodes": cc, "polys": cc_polys, "color": si, "level": -1}
    ccs_data.update(border_info)

    cis = ccs_data.keys()
    eids_to_ccs = {}
    adjacent = dict([(ci, set()) for ci in cis])
    for ci in cis:
        ccs_data[ci]["poly_sets"] = []
        for p in ccs_data[ci]["polys"]:
            ccs_data[ci]["poly_sets"].append(set([abs(eid) for eid in p]))
            for eid in ccs_data[ci]["poly_sets"][-1]:
                if eid not in eids_to_ccs:
                    eids_to_ccs[eid] = []
                else:
                    adjacent[eids_to_ccs[eid][0]].add(ci)
                    adjacent[ci].add(eids_to_ccs[eid][0])
                eids_to_ccs[eid].append(ci)

    # compute reachability level starting from outside level=0
    queue = [k for (k, vs) in ccs_data.items() if vs["level"] == 0]
    nextlevel = set()
    level = 1
    while len(queue) > 0:
        for ci in queue:
            for cj in adjacent[ci]:
                if ccs_data[cj]["level"] < 0:
                    ccs_data[cj]["level"] = level
                    nextlevel.add(cj)
        queue = nextlevel
        nextlevel = set()
        level += 1

    # From outside to inner most components find exterior border of component
    cks = sorted(ccs_data.keys(), key=lambda x: ccs_data[x]["level"])
    for ck, data in ccs_data.items():
        data["exterior_polys"] = [i for i in range(len(data["polys"]))]
        if data["level"] == 0:
            data["exterior_polys"] = []
            continue
        if len(data["polys"]) < 2:
            continue
        interior_dots = set()
        for k in adjacent[ck]:
            if data["level"] < ccs_data[k]["level"]:
                interior_dots.update(*ccs_data[k]["poly_sets"])

        if len(interior_dots) > 0:
            data["exterior_polys"] = [pi for pi, poly in enumerate(data["poly_sets"]) if len(set(poly).difference(interior_dots)) > 0]
    return ccs_data, cks, adjacent


# EDGES_FIELDS_LIST = []
EDGES_FIELDS_LIST = [('eid', "%d"), ('edge', "(%s,%s)"), ('flat_edge', "(%s,%s)"), ("types", "%s"),
                     ('nodes', "(%s,%s)"), ('pos', "(%s,%s)"),
                     ('nodes_cut', "(%s,%s)"), ('pos_cut', "(%s,%s)"),
                     ('uncut_eid', "%d"), ('cut_eid', "%d"), ('far', "%d"), ('n_closer', "%d")]
EDGES_FIELDS_MAP = dict(EDGES_FIELDS_LIST)


def preamble_for_edges(list_edges):
    entries = [(k, v) for (k, v) in list_edges[0].items() if re.match("last_", k) or (k == "nb_nodes")]
    entries_more = [(k, list_edges[0][k]) for k in ["dst_type", "source"]]
    if "grid_percentile" in list_edges[0]:
        entries_more.append(("grid_percentile", "(%s,%s)" % tuple(list_edges[0]["grid_percentile"])))
    for k in ["bbox", "ebbox"]:
        if k in list_edges[0]:
            entries_more.append((k, "(%s,%s,%s,%s)" % tuple(list_edges[0][k])))
    entries.sort(key=lambda x: x[1])
    return "# LIST OF EDGES\t%s" % " ".join(["%s=%d" % (k, v) for (k, v) in entries]+["%s=%s" % (k, v) for (k, v) in entries_more])


def header_for_edges(list_edges):
    return "\t".join([k for k, f in EDGES_FIELDS_LIST])


def parse_preamble(line):
    list_edges = [{"edge": None, TIK: {}, "nb_nodes": 0}]
    if re.match("# LIST OF EDGES\t", line):
        parts = line.strip().split("\t")[1].split(" ")
        for part in parts:
            pp = part.strip().split("=")
            if pp[0] in ["bbox", "ebbox"]:
                list_edges[0][pp[0]] = numpy.array(eval(pp[1]))
            elif pp[0] in ["grid_percentile"]:
                list_edges[0][pp[0]] = tuple(eval(pp[1]))
            elif pp[0] in ["dst_type", "source"]:
                list_edges[0][pp[0]] = pp[1]
            elif pp[0] in ["globe"]:
                list_edges[0][pp[0]] = (pp[1] == "True")
            else:
                list_edges[0][pp[0]] = int(pp[1])
    return list_edges


def parse_header(line):
    parts = line.strip().split("\t")
    return [(p, i) for i, p in enumerate(parts)]


def str_for_edge(list_edges, eid):
    edge = {"eid": eid}
    edge.update(list_edges[eid])
    for field in ['types']:
        if field in edge:
            edge[field] = "(%s)" % ",".join(edge[field])
    for field in ['nodes', 'pos', 'nodes_cut', 'pos_cut']:
        if field in edge:
            if len(edge[field]) == 1:
                edge[field] = (edge[field][0], "-")
            elif len(edge[field]) == 2:
                edge[field] = (edge[field][0], edge[field][1])
            elif len(edge[field]) > 2:
                print("EDGE", eid, edge)
                raise ValueError("More than two nodes to an edge!")
            else:
                edge[field] = None
    return "\t".join([f % edge[k] if edge.get(k) is not None else "" for k, f in EDGES_FIELDS_LIST])


def parse_edge(line, head):
    parts = line.strip().split("\t")
    eid = None
    edge = {"edge": None, "nodes": [], "types": set()}
    for p, i in head:
        if i < len(parts) and len(parts[i]) > 0:
            if p in ["edge", "flat_edge"]:
                v = eval(parts[i])
            elif p in ["types"]:
                v = set(parts[i].strip("()").split(","))
            elif EDGES_FIELDS_MAP[p] == "(%s,%s)":
                v = [int(vx) for vx in parts[i].strip("()").split(",") if vx != "-"]
            elif re.match("%[0-9.]*f$", EDGES_FIELDS_MAP[p]):
                v = float(parts[i])
            elif re.match("%[0-9.]*d$", EDGES_FIELDS_MAP[p]):
                v = int(parts[i])
            else:
                v = parts[i]

            if p == "eid":
                eid = v
            else:
                edge[p] = v
    return eid, edge


def write_edges(fp, list_edges):
    close = False
    if not isinstance(fp, io.TextIOBase):
        fp = open(fp, "w")
        close = True
    fp.write(preamble_for_edges(list_edges)+"\n")
    fp.write(header_for_edges(list_edges)+"\n")
    for eid in range(1, len(list_edges)):
        fp.write(str_for_edge(list_edges, eid)+"\n")
    if close:
        fp.close()


def read_edges(filename):
    list_edges = None
    close = False
    f, fcl = getFp(filename)
    if f is not None:
        head = None
        for line in f:
            if list_edges is None:
                list_edges = parse_preamble(line)
            elif head is None:
                head = parse_header(line)
            else:
                eid, edge = parse_edge(line, head)
                if eid is not None and eid == len(list_edges):
                    list_edges.append(edge)
    if fcl:
        f.close()
    return list_edges


def build_from_edges(list_edges):
    map_edges, polys, polys_cut = ({}, {}, {})
    tmp_polys, tmp_polys_cut = ({}, {})
    for eid in range(1, len(list_edges)):
        edge = list_edges[eid]
        map_edges[edge["edge"]] = eid
        if len(edge["types"]) > 0:
            update_types(list_edges, eid, edge["types"], from_edge=True)
        for nkey, pkey, trg in [("nodes_cut", "pos_cut", tmp_polys_cut), ("nodes", "pos", tmp_polys)]:
            if edge.get(nkey) is not None and edge.get(pkey) is not None and len(edge.get(nkey)) > 0 and len(edge.get(nkey)) == len(edge.get(pkey)):
                for ni, n in enumerate(edge[nkey]):
                    p = edge[pkey][ni]
                    if n not in trg:
                        trg[n] = []
                    trg[n].append((abs(p), numpy.sign(p)*eid))
    nodes_polys = set().union(tmp_polys.keys(), tmp_polys_cut.keys())
    for node in nodes_polys:
        if node in tmp_polys:
            tmp_polys[node].sort()
        if node in tmp_polys_cut:
            tmp_polys_cut[node].sort()
        if tmp_polys.get(node) == tmp_polys_cut.get(node):
            polys[node] = [p[1] for p in tmp_polys[node]]
        else:
            if node in tmp_polys:
                polys[node] = [p[1] for p in tmp_polys[node]]
            if node in tmp_polys_cut:
                polys_cut[node] = [p[1] for p in tmp_polys_cut[node]]
    return map_edges, list_edges, polys, polys_cut


def read_edges_and_co(filename):
    list_edges = read_edges(filename)
    return build_from_edges(list_edges)


def flatToGlobe(coord):
    rx, ry = numpy.radians(coord)
    return (numpy.cos(ry) * rx, ry)


def globeToFlat(coord, clip=False):
    x, y = (numpy.degrees(coord[0]/numpy.cos(coord[1])), numpy.degrees(coord[1]))
    if clip is not None:
        if x < clip[0]:
            x = clip[0]
        elif x > clip[2]:
            x = clip[2]
        if y < clip[1]:
            y = clip[1]
        elif y > clip[3]:
            y = clip[3]
    return x, y


def coordsPointsToGlobe(PointsMap):
    return dict([(i, flatToGlobe(coord)) for i, coord in PointsMap.items()])


def coordsEdgesFromGlobe(list_edges, clip=None):
    for i in range(1, len(list_edges)):
        globe_edge = list_edges[i]["edge"]
        flat_edge = (globeToFlat(globe_edge[0], clip=clip), globeToFlat(globe_edge[1], clip=clip))
        list_edges[i]["flat_edge"] = flat_edge
    return list_edges


# Example of preparing and plotting polygons
if __name__ == "__main__":

    COORDS_FILE = "coords_bckg-PETIT.csv"
    # COORDS_FILE = "coords_bckg.csv"
    coords_bckg, rnames_bckg = read_coords_csv(COORDS_FILE)

    PointsMap = {}
    PointsIds = {}
    for i, coord in enumerate(rnames_bckg):
        PointsIds[i] = rnames_bckg[i]
        PointsMap[i] = coords_bckg[i]

    dst_type = "globe"
    # dst_type = "flat"
    hgrid_percentile = 95
    map_edges, list_edges, polys, polys_cut, ebbox = prepare_edges_dst(PointsMap, hgrid_percentile=hgrid_percentile, dst_type=dst_type)

    # ## FULL MAP BY TYPES
    # bounds = [1, list_edges[0]["last_org"], list_edges[0]["last_cut"], list_edges[0]["last_isolated"], len(list_edges)]
    # elements = [(bounds[0], bounds[1], "-", "k"), (bounds[1], bounds[2], "-", "g"), (bounds[2], bounds[3], "-", "m"), (bounds[3], bounds[4], "-", "r")]
    # for (bdw, bup, linestyle, color) in elements:
    #     plot_edges(range(bdw, bup), list_edges, linestyle=linestyle, color=color, kedge="flat_edge")
    # # xs, ys = zip(*PointsMap.values())
    # # plt.plot(xs, ys, "bo")
    # # for pi, cc in PointsMap.items():
    # #     plt.text(cc[0], cc[1], pi)
    # plot_show(True)

    # fname = "/home/egalbrun/short/test_edges.txt"
    # write_edges(fname, list_edges)
    # rmap_edges, rlist_edges, rpolys, rpolys_cut = read_edges_and_co(fname)

    # npairs = prepare_node_pairs(list_edges, max_nid=200)

    after_cut = True
    pp = prepare_areas_polys(polys, polys_cut, after_cut)
    # for p in [4, 12]:
    pds, pps = zip(*pp.items())
    # PointsRad = coordsPointsToGlobe(PointsMap)
    for p in pds:
        # plt.figure()
        # plt.text(PointsRad[p][0], PointsRad[p][1], p)
        # plt.plot(PointsRad[p][0], PointsRad[p][1], "ro")
        plt.text(PointsMap[p][0], PointsMap[p][1], p)
        plt.plot(PointsMap[p][0], PointsMap[p][1], "ro")
        plot_esets_colordered([pp[p]], list_edges, simple=False, kedge="flat_edge")
        # plot_show(False)
    #plot_esets_colordered(pps, list_edges, simple=False, kedge="flat_edge")
    plot_show(True)

    border_info, nodes_graph = prepare_areas_helpers(map_edges, list_edges, after_cut)
    # Xs = [o["polys"][0] for o in border_info.values()]
    # plot_esets_colordered(Xs, list_edges, simple=True, kedge="flat_edge")
    # plot_show()

    ccs_data, cks, adjacent = prepare_areas_data([10, 8], list_edges, pp, border_info, nodes_graph)

    colors = {-2: "#AAAAAA", -1: "white", 0: "r"}
    pp_data = {"ccs_data": ccs_data, "cks": cks, "adjacent": adjacent}
    if pp_data is not None and "cks" in pp_data:
        for cii, ck in enumerate(pp_data["cks"]):
            pp_polys = pp_data["ccs_data"][ck]["polys"]

            for pi in pp_data["ccs_data"][ck]["exterior_polys"]:
                plot_filled(pp_polys[pi], list_edges, color=colors.get(pp_data["ccs_data"][ck]["color"], "k"))
    plot_show(True)

    if dst_type == "globe":
        colors = {-2: "#AAAAAA", -1: "white", 0: "r"}
        pp_data = {"ccs_data": ccs_data, "cks": cks, "adjacent": adjacent}
        if pp_data is not None and "cks" in pp_data:
            for cii, ck in enumerate(pp_data["cks"]):
                pp_polys = pp_data["ccs_data"][ck]["polys"]

                for pi in pp_data["ccs_data"][ck]["exterior_polys"]:
                    plot_filled(pp_polys[pi], list_edges, color=colors.get(pp_data["ccs_data"][ck]["color"], "k"), kedge="flat_edge")

    # for ci, dt in border_info.items():
    #     if dt["level"] == 0:
    #         linewidth = 2
    #     else:
    #         linewidth = 1
    #     plot_edges(dt["polys"][0], list_edges, linewidth=linewidth)

    # # eids_sel =  get_border_eids(list_edges)
    # eids_sel = list_edges[0][TIK][OUTER_MAP[True]]
    # for eid in range(1, len(list_edges)):
    #     linewidth = 1
    #     if eid in eids_sel:
    #         linewidth = 2
    #     nbN = len(list_edges[eid].get("nodes_cut", []))
    #     if nbN > 0:
    #         if nbN == 2:
    #             color = "b"
    #         elif nbN == 1:
    #             color = "g"
    #         else:
    #             color = "r"
    #         plot_edges([eid], list_edges, color=color, linewidth=linewidth)

    # plot_edges_colordered(list_edges[0][TIK][OUTER_MAP[True]], list_edges, "outer")
    # plot_edges(get_border_eids(list_edges), list_edges, linestyle=":", color="k")

    plot_show(True)
