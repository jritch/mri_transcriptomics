import os,sys,itertools,shutil

neuro_expresso_folder = "other gene lists/NeuroExpresso"
region = "Cortex"

pyramidal_deep = os.path.join(neuro_expresso_folder,"PyramidalDeep",region)

other_cell_lists = filter(lambda x: "Pyra" not in x, os.listdir(pyramidal_deep))
other_cell_lists = [os.path.join(pyramidal_deep,x) for x in other_cell_lists]

pyramidal_cell_list = os.path.join(neuro_expresso_folder,"CellTypes",region,"Pyramidal")

pyramidal_cell_list

cell_lists = other_cell_lists + [pyramidal_cell_list]

new_names = ["other gene lists/" + ".".join(itertools.compress(x.split("/"),[0,1,0,1,1])) for x in cell_lists];

for i in range(len(new_names)):
    name = new_names[i]

    words = name.split(".")

    if "Pyra" in words[-1] or "Gaba" in words[-1]:
        words.append(words[-1])
        words[-2] = "Neuron"
    new_names[i] = ".".join(words)


for i in range(len(cell_lists)):
    shutil.copy(cell_lists[i],new_names[i])
