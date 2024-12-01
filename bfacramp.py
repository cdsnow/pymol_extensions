from pymol import cmd, stored, math

def global_bounds():
    objects = cmd.get_object_list()
    print('Inspecting',len(objects),'objects')

    stored.bfacts=[]
    for obj in objects:
        cmd.iterate(obj, 'stored.bfacts.append(b)')
    print(min(stored.bfacts),'to',max(stored.bfacts))
    round_min_to_10 = 10*int(min(stored.bfacts)/10)
    #round_max_to_10 = 10*round(max(stored.bfacts)/10)
    round_max_to_10 = 10*math.ceil(max(stored.bfacts)/10)
    return round_min_to_10, round_max_to_10

def Bramp(mol=None, order='rainbow'):
    """
    mol = any object selection (within one single object though)
    
    example: Bramp 1LVM and chain A
    """
    objects = cmd.get_object_list()
    if len(objects) == 1 and mol is None:
        obj=objects[0]
    elif mol in objects:
        obj=cmd.get_object_list(mol)[0]
    elif len(objects)>1 and mol is None:
        print(f'More than one object. Please specify one')
        return
    else:
        print(f'Could not find {mol} in the list of objects')
        return

    stored.bfacts=[]
    cmd.iterate(obj, 'stored.bfacts.append(b)')
    print(min(stored.bfacts),'to',max(stored.bfacts))
    round_min_to_10 = 10*int(min(stored.bfacts)/10)
    #round_max_to_10 = 10*round(max(stored.bfacts)/10)
    round_max_to_10 = 10*math.ceil(max(stored.bfacts)/10)

    if order == 'rainbow': 
        cmd.ramp_new(f'scalebar_{obj}', obj, [round_min_to_10,round_max_to_10], 'rainbow')
        print(f'coloring {obj} by b-factor from {round_min_to_10} to {round_max_to_10}')
        cmd.spectrum('b',selection=obj, minimum=round_min_to_10,maximum=round_max_to_10, palette='rainbow')
    elif order == 'AF2':
        colorlist = ['red','yellow','green','cyan','blue']
        rgblist = [list(cmd.get_color_tuple(cc)) for cc in colorlist]
        print('rgblist',rgblist)
        cmd.ramp_new(f'scalebar_{obj}', obj, [round_min_to_10,round_max_to_10], rgblist)
        mypalette = '_'.join(colorlist)
        print('mypalette',mypalette)
        print(f'coloring {obj} by b-factor from {round_min_to_10} to {round_max_to_10}')
        cmd.spectrum('b',selection=obj, minimum=round_min_to_10,maximum=round_max_to_10, palette=mypalette)
    else:
        print('color order is not recognized')
        return
    cmd.recolor()

def BrampAll(order='rainbow'):
    """
    example: BrampAll
    """
    objects = cmd.get_object_list()
    round_min_to_10, round_max_to_10 = global_bounds()

    if order == 'rainbow': 
        cmd.ramp_new('Bscalebar', objects[0], [round_min_to_10,round_max_to_10], 'rainbow')
        for obj in objects:
            print(f'coloring {obj} by b-factor from {round_min_to_10} to {round_max_to_10}')
            cmd.spectrum('b',selection=obj, minimum=round_min_to_10,maximum=round_max_to_10, palette='rainbow')
    elif order == 'AF2':
        colorlist = ['red','yellow','green','cyan','blue']
        rgblist = [list(cmd.get_color_tuple(cc)) for cc in colorlist]
        print('rgblist',rgblist)
        mypalette = '_'.join(colorlist)
        print('mypalette',mypalette)
        cmd.ramp_new('Bscalebar', objects[0], [round_min_to_10,round_max_to_10], rgblist)
        for obj in objects:
            print(f'coloring {obj} by b-factor from {round_min_to_10} to {round_max_to_10}')
            cmd.spectrum('b',selection=obj, minimum=round_min_to_10,maximum=round_max_to_10, palette=mypalette)
    else:
        print('color order is not recognized')
        return
    cmd.recolor()

def AFscalebar(mol=None):
    Bramp(mol=mol, order='AF2')

def AFscalebarAll():
    BrampAll(order='AF2')

cmd.extend("Bramp", Bramp);
cmd.extend("BrampAll", BrampAll);
cmd.extend("AFscalebar", AFscalebar);
cmd.extend("AFscalebarAll", AFscalebarAll);

