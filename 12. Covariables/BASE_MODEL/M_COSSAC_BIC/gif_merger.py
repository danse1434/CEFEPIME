import imageio
import numpy as np    

#Create reader object for the gif
gif1 = imageio.get_reader('G1_anim.gif')
gif2 = imageio.get_reader('G2_anim.gif')

#If they don't have the same number of frame take the shorter
number_of_frames = min(gif1.get_length(), gif2.get_length()) 

#Create writer object
new_gif = imageio.get_writter('output.gif')

for frame_number in range(number_of_frame):
    img1 = gif1.get_next_data()
    img2 = gif2.get_next_data()
    #here is the magic
    new_image = np.hstack((img1, img2))
    new_gif.append_data(new_image)

gif1.close()
gif2.close()    
new_gif.close()