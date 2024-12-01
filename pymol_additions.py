from pymol import cmd
from pymol.cgo import *
from glob import glob
try: from supercell import *
except: pass

cylinder_template = """
obj = [CYLINDER]
obj.extend(%s)
obj.extend(%s)
obj.extend([%f])
obj.extend(%s)
obj.extend(%s)
"""

def reverse_states(objstr):
    numstates = cmd.count_states(objstr)
    cmd.set_state_order(objstr, range(numstates,0,-1))

def color_each_spectrum():
    for i in cmd.get_names('objects'): cmd.spectrum('count',selection=i)

def scalebar(start, end, radius=3, color=[0.0, 0.0, 0.0], name='scalebar'):
    ## Usage:
    ## scalebar((220,200,36), (320,200,36))
    outstr = ""
    loA = '[%s]' % ','.join(map(str, start))
    hiB = '[%s]' % ','.join(map(str, end))
    outstr += cylinder_template % (loA, hiB, radius, color, color)
    outstr += "cmd.load_cgo(obj, '%s')" % name
    open('scalebar.pml','w').write(outstr)
    print('Wrote scalebar.pml')

def splitchains(selection,chainlabels):
    for chain in chainlabels:
        print(chain)
        cmd.save('%s.%s.pdb' % (selection,chain), '%s and chain %s' % (selection,chain))
        cmd.load('%s.%s.pdb' % (selection,chain))

def globload(pattern='*.pdb'):
    files = glob(pattern)
    for file in files: cmd.load(file)

def hide_bb_sticks():
  cmd.hide("sticks","name C or name O or (name N and not resn PRO)")

def labeldna():
  cmd.label("name C5'", 'chain+resi')

def hide_bb_lines():
  cmd.hide("lines","name C or name O or (name N and not resn PRO)")

def hide_bb():
  cmd.hide("sticks","name C or name O or (name N and not resn PRO)")
  cmd.hide("lines","name C or name O or (name N and not resn PRO)")

def max_cartoon():
  cmd.hide("everything")
  cmd.set("cartoon_cylindrical_helices","on")
  cmd.set("cartoon_helix_radius",1.5)
  cmd.set("cartoon_smooth_loops","on")
  cmd.set("cartoon_flat_sheets","on")
  cmd.set("cartoon_discrete_colors","on")
  cmd.color("wheat","ss l+")
  cmd.color("raspberry","ss h")
  cmd.color("orange","ss s")
  cmd.show("cartoon")

def modern_cartoon():
  util.performance(0)
  cmd.bg_color('white')
  cmd.hide("everything")
  cmd.set("cartoon_helix_radius",1.5)
  cmd.set("cartoon_flat_sheets","on")
  cmd.set("cartoon_loop_radius",0.25)
  util.chainbow()
  cmd.show("cartoon")
  cmd.set("ambient",1)
  cmd.set("specular",0)
  cmd.set("ray_shadow",0)
  cmd.set("ray_trace_gain",0.0001)

def blob(obj, resolution=8):
    cmd.set('gaussian_b_floor', 50)
    cmd.set('gaussian_resolution', resolution)
    mapname = 'map.%s'%str(obj)
    surfname = 'surf.%s'%str(obj)
    #cmd.map_new(mapname,'gaussian', grid=3, selection=obj, buffer=10)
    cmd.map_new(mapname,'gaussian', grid=1, selection=obj, buffer=10)
    cmd.isosurface(surfname,mapname,1.0)

def cgo_arrow(atom1='pk1', atom2='pk2', radius=0.5, gap=0.0, hlength=-1, hradius=-1,
              color='blue red', name=''):
    '''
DESCRIPTION
    Create a CGO arrow between two picked atoms.

ARGUMENTS
    atom1 = string: single atom selection or list of 3 floats {default: pk1}
    atom2 = string: single atom selection or list of 3 floats {default: pk2}
    radius = float: arrow radius {default: 0.5}
    gap = float: gap between arrow tips and the two atoms {default: 0.0}
    hlength = float: length of head
    hradius = float: radius of head
    color = string: one or two color names {default: blue red}
    name = string: name of CGO object
    '''
    from chempy import cpv

    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)

    try:
        color1, color2 = color.split()
    except:
        color1 = color2 = color
    color1 = list(cmd.get_color_tuple(color1))
    color2 = list(cmd.get_color_tuple(color2))

    def get_coord(v):
        if not isinstance(v, str):
            return v
        if v.startswith('['):
            return cmd.safe_list_eval(v)
        return cmd.get_atom_coords(v)

    xyz1 = get_coord(atom1)
    xyz2 = get_coord(atom2)
    normal = cpv.normalize(cpv.sub(xyz1, xyz2))

    if hlength < 0:
        hlength = radius * 3.0
    if hradius < 0:
        hradius = hlength * 0.6

    if gap:
        diff = cpv.scale(normal, gap)
        xyz1 = cpv.sub(xyz1, diff)
        xyz2 = cpv.add(xyz2, diff)

    xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

    obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
          [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
          [1.0, 0.0]

    if not name:
        name = cmd.get_unused_name('arrow')

    cmd.load_cgo(obj, name)

## Expose to the pymol command line
cmd.extend("blob",blob)
cmd.extend("max_cartoon",max_cartoon)
cmd.extend("modern_cartoon",modern_cartoon)
cmd.extend("hide_bb_sticks",hide_bb_sticks)
cmd.extend("hide_bb_lines",hide_bb_lines)
cmd.extend("hide_bb",hide_bb)
cmd.extend("loadglob",globload)
cmd.extend("globload",globload)
cmd.extend("splitchains",splitchains)
cmd.extend('cgo_arrow', cgo_arrow)
cmd.extend('scalebar',scalebar)
cmd.extend('color_each_spectrum',color_each_spectrum)
cmd.extend('reverse_states',reverse_states)
