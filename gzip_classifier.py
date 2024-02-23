import gzip
import numpy as np
import os
from random import shuffle
import pandas as pd
import os.path
#to-dos
#the major change this code needs is in how it iterates over test and trainign set
#and what the test and training set are
#it needs to take in folders containing .fa files
#iterate over them to get ncd between the folders, then output
#the ncd to a csv file
#the knn classifier is not necessary, but we can include if we want, focus
#on getting ncd calcs first tho

splatter_dirname = "../data/splatter_500_fastas"
astrocyte_dirname = "../data/sc_astrocyte_500_fastas"
epithileal_dirname= "../data/epithileal_fastas"
ext = ".fa"
k = 2
training_splatter = {}
training_astrocyte = {}
test_splatter = {}
test_astrocyte = {}
training_set = []
test_set = []
ordered_list = []

for fa_files in os.listdir(splatter_dirname):
    if fa_files.endswith(ext):
        ordered_list.append(splatter_dirname + "/" + fa_files)
    else:
        continue

for fa_files in os.listdir(astrocyte_dirname):
    if fa_files.endswith(ext):
        ordered_list.append(astrocyte_dirname + "/" + fa_files)
    else:
        continue

shuffle(ordered_list)
counter = 0
for test_fa in ordered_list[0:50]: 
    with open(test_fa, "r") as file:
    #training_fa.split('_')[0]
        temp_tuple = (file.read(), os.path.basename(test_fa))
        if counter == 0:
            counter += 1
    file.close()
    test_set.append(temp_tuple)

training_list = []

for fa_files in os.listdir(epithileal_dirname):
    if fa_files.endswith(ext):
        training_list.append(epithileal_dirname + "/" + fa_files)
    else:
        continue

shuffle(training_list)
for training_fa in training_list[0:100]:
    with open(training_fa, "r") as file:
        temp_tuple = (file.read(), "test_epithileal_cell") 
    file.close()
    training_set.append(temp_tuple)

for (x1, name) in test_set:
    print(name)
    Cx1 = len(gzip.compress(x1.encode()))
    distance_from_x1=[]
    for (x2, name) in training_set:
        Cx2 = len(gzip.compress(x2.encode()))
        x1x2 = " ".join([x1, x2])
        Cx1x2 = len(gzip.compress(x1x2.encode()))
        ncd = (Cx1x2 - min(Cx1, Cx2)) / max(Cx1, Cx2)
        distance_from_x1.append(ncd)
    sorted_idx = np.argsort(np.array(distance_from_x1))
    df = pd.DataFrame(distance_from_x1)
    df.to_csv("../data/ncd_calcs/" + name + "_ncd_calc.csv")
    #print(distance_from_x1)
    #top_k_class = training_set[sorted_idx[:k]]
    #predict_class = max(set(top_k_class), key=top_k_class.count)

#np.savetxt("../data/sorted_idx.txt", sorted_idx)
#np.savetxt("../data/top_k_class.txt", top_k_class)
#np.savetxt("../data/predict_class", predict_class)
