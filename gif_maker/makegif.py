import os
#import imageio as io
import re
from PIL import Image, ImageDraw, ImageFont
import argparse #argument parser, handles the arguments passed by command line
sys.path.append('/home/jkhadka/transferdata/scripts/simulation_functions/')
import simulation_functions as sf
#############################################################
#setting up the arguments to be passed 
parser = argparse.ArgumentParser()#parser
#parser.add_argument('-l','--location', help="The location of targetformmatrix and coordinate file",type = string)
parser.add_argument('-s',"--starttime", help="Start of simulation time",default = 1, type = int)
parser.add_argument('-e',"--endtime", help="End of simulation area",default = None,  type = int)
parser.add_argument('-d',"--timestep", help="time step for simulation", type = int,default = 1)
parser.add_argument('-u',"--duration",help="duration for each frame to be displayed in gif", default  = 10, type = int)
parser.add_argument('-c',"--surfacearea", help = "make gif for simulation till this surfacearea", type = float, default = None)


## Getting the arguments 
args = parser.parse_args()

duration = args.duration
startstep = args.starttime
endstep= args.endtime
timestep = args.timestep
surfacearea = args.surfacearea
#################################################
# if surface area is used, calculate timstep
#################################################
if surfacearea:
    endStep = 2000
    startStep = 1000
    endarea,tissueSurfaceArea = sf.getTimeStep(surfacearea, endStep, startStep, stepsize = 10)
    areastep = timestep#use area step instead of time for plotting
#####################################################################################################
#key to sort the file_names in order
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts
############################################################


if surfacearea:
    file_names = []
    area_list = []
    startarea = sf.getSurfaceArea(1)
    for steparea in range(startarea, endarea, int(areastep)):
            step,tissueSurfaceArea = sf.getTimeStep(steparea, endStep, laststep, stepsize = 5)
            ########################################################################
            if not os.path.isfile("qdObject_step=%03d.obj"%step):#check if file exists
                break
            file_names.append('surface_time=%03.d.png'%(step))
            area_list.append(tissueSurfaceArea)
else:
    file_names = sorted((fn for fn in os.listdir('.') if fn.startswith('surface')), key = numericalSort)
"""#gif writer
with io.get_writer('python_growth.gif', mode='I', duration=0.1) as writer:
    for filename in file_names:
        image = io.imread(filename)
        writer.append_data(image)
#writer.close()
"""

######################################################################
def gen_frame(path,surfacearea=None, counter=None):
    im = Image.open(path)
    alpha = im.getchannel('A')

    # Convert the image into P mode but only use 255 colors in the palette out of 256
    im = im.convert('RGB').convert('P', palette=Image.ADAPTIVE, colors=255)

    # Set all pixel values below 128 to 255 , and the rest to 0
    mask = Image.eval(alpha, lambda a: 255 if a <=20 else 0)

    # Paste the color of index 255 and use alpha as a mask
    im.paste(255, mask)

    # The transparency index is 255
    im.info['transparency'] = 255
    
    # adding area detail if needed
    if surfacearea:
        draw = ImageDraw.Draw(im)
        font = ImageFont.truetype("sans-serif.ttf", 16)
        draw.text((50,10),r"$A_T = %.1f$"%surfacearea,font = font)
    return im
######################################################################
frames = []
file_names = file_names[startstep:endstep:timestep]#[::5]
counter  = 0
for filename in file_names:
    frames.append(gen_frame(filename,surfacearea=surfacearea, counter=counter))
    counter += 1
######################################################################
name = os.path.basename(os.getcwd())
if surfacearea:
    frames[0].save(name+"_surfacearea"+'.gif', save_all=True, append_images=frames[1::],loop = 5, duration = duration)
else:
    frames[0].save(name+"_timestep"+'.gif', save_all=True, append_images=frames[1::],loop = 5, duration = duration)
