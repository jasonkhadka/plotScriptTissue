import os
#import imageio as io
import re
from PIL import Image
#############################################################
#key to sort the file_names in order
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts
############################################################
file_names = sorted((fn for fn in os.listdir('.') if fn.startswith('surface')), key = numericalSort)

"""#gif writer
with io.get_writer('python_growth.gif', mode='I', duration=0.1) as writer:
    for filename in file_names:
        image = io.imread(filename)
        writer.append_data(image)
#writer.close()
"""

######################################################################
def gen_frame(path):
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

    return im
######################################################################
frames = []
file_names = file_names[::5]
for filename in file_names:
	frames.append(gen_frame(filename))
######################################################################
frames[0].save('tissue_growth.gif', save_all=True, append_images=frames[1::],loop = 5, duration = 10)
