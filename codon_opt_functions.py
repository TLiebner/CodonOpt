from numpy.random import randint
from numpy.random import rand
import time
from scipy.optimize import minimize
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
import pandas as pd
import numpy as np
import re
from collections import Counter
import math
import itertools
import pprint
import random
import copy 
import pygad
import statistics

def parse(file):
    file = re.sub(r'\( ', '', file)
    file = re.sub(r'\) ', '', file)
    file = re.sub(r'\n', '', file)
    file = re.sub(r'[(|)]', ' ', file)
    file_split = file.split(" ")

    rows = []
    row = []
    counter = 1
    for i in file_split:
        #print(i, counter)
        if i == "":
            continue
        if counter == 5:
            #print(row)
            rows.append(row)
            row = []
            counter = 1
            continue
        if counter < 5:
            row.append(i)
            counter += 1

    profile = pd.DataFrame(columns = ('Codon', 'AA', 'Frequency', 'Usage'))
    profile = pd.DataFrame(rows).rename({0: "Codon", 1: "AA", 2: "Frequency", 3: "Usage"}, axis='columns')
    return(profile)

def read_and_parse(file_path):
    file = open(file_path, 'r').read()
    file = re.sub(r'\( ', '', file)
    file = re.sub(r'\) ', '', file)
    file = re.sub(r'\n', '', file)
    file = re.sub(r'[(|)]', ' ', file)
    file_split = file.split(" ")

    rows = []
    row = []
    counter = 1
    for i in file_split:
        #print(i, counter)
        if i == "":
            continue
        if counter == 5:
            #print(row)
            rows.append(row)
            row = []
            counter = 1
            continue
        if counter < 5:
            row.append(i)
            counter += 1

    profile = pd.DataFrame(columns = ('Codon', 'AA', 'Frequency', 'Usage'))
    profile = pd.DataFrame(rows).rename({0: "Codon", 1: "AA", 2: "Frequency", 3: "Usage"}, axis='columns')
    return(profile)

def AA_look_up(AA_seq, profile):
    AA_seq = [x for x in AA_seq]
    AA_seq = pd.Series(AA_seq)
    
    counts = pd.DataFrame(AA_seq.value_counts())
    counts.index.name = 'AA'
    counts.reset_index(inplace = True)
    counts = counts.rename(columns = {0 : "Total"})
    look_up_df = pd.merge(counts, profile, left_on = "AA", right_on = "AA")
    look_up_df_type_dict = {
        'AA' : str,
        'Total': float,
        'Codon' : str,
        'Frequency': float,
        'Usage' : str
                          }
    look_up_df = look_up_df.astype(look_up_df_type_dict)
    look_up_df['#'] = np.round(look_up_df['Total'] * look_up_df['Frequency'])
    look_up_df['Used'] = 0

    look_up_arr = np.asarray(look_up_df)
    #print(np.unique(AA_seq))
    look_up_arr_list = []
    for AA in np.unique(AA_seq):
        sub_look_up_arr = look_up_arr[look_up_arr[:,0] == AA]
        #print(sub_look_up_arr)
        total = sub_look_up_arr[0,1]
        col_sum = sub_look_up_arr.sum(axis = 0)[5]
        #print(sub_look_up_arr, total, col_sum)
        #col_sum = col_sum - 5
        counter = total - col_sum
        if counter != 0:
            modifier = counter/abs(counter)
        #counter = 10
        #print(counter, modifier)
        #print(counter)
        #y = 
        freq = np.sort(sub_look_up_arr[:,3])
        freq[-1]
        if counter > 0:
            used_highest = sub_look_up_arr[sub_look_up_arr[:,5] == max(sub_look_up_arr[:,5]), 6] 
            sub_look_up_arr[sub_look_up_arr[:,5] == max(sub_look_up_arr[:,5]), 6] = used_highest + 1
        if counter < 0:
            used_highest = sub_look_up_arr[sub_look_up_arr[:,5] == min(sub_look_up_arr[:,5]), 6] 
            sub_look_up_arr[sub_look_up_arr[:,5] == min(sub_look_up_arr[:,5]), 6] = used_highest - 1    
        x = 0*len(sub_look_up_arr[:,6])
        #print(sub_look_up_arr[:,6] + 1, x)

        while counter != 0:
            for row in sub_look_up_arr:
                #print(row)
                if counter == 0:
                    #print(sub_look_up_arr)
                    break
                if row[6] != 0:
                    continue


                sub_look_up_arr[sub_look_up_arr[:,2] == row[2], 5] += 1 * modifier
                sub_look_up_arr[sub_look_up_arr[:,2] == row[2], 6] += 1 * modifier

                #print(sub_look_up_arr, counter)

                if counter > 0:
                    counter = counter - 1
                else:
                    counter = counter + 1

            sub_look_up_arr[:,6] = 0
        sub_look_up_arr[:,6] = 100
        look_up_arr_list.append(sub_look_up_arr) 
    look_up_array = np.concatenate(look_up_arr_list)
    
    return(look_up_array)

def subsequence_array(look_up_array, AA_seq):
    AA_uniq = np.unique([x for x in AA_seq])
    df_list = []

    for AA in AA_uniq:
        sub_look_up = look_up_array[look_up_array[:,0] == AA]
        index = np.where(sub_look_up[:,3] == max(sub_look_up[:,3]))
        index = np.concatenate(index)[0]
        col_sum = sub_look_up.sum(axis = 0)[4]
        #print(sub_look_up)


        #print(len(df[:,0]))
        #print(col_sum)
        subset = []
        #print(sub_look_up)
        #print(range(len(sub_look_up))
        for i in range(len(sub_look_up)):
            subset_AA = sub_look_up[i]
            #print(subset_AA)
            if subset_AA[5] != 0:
                sequence = np.arange(start = 0, stop = subset_AA[5])
                for x in sequence:
                    subset.append(subset_AA[2])
        matches = [i for i,x in enumerate(AA_seq) if x == AA]
        data = {'Index':  matches,
                'Codon': subset}        
        df = pd.DataFrame(data)
        df = np.asarray(df)
        df_list.append(df)
    return(df_list)

def mean_dist_fitness_func_arr(array):
    distance_list = []
    for i in array:
        a = i
        array = array[1:]
        for j in array:
            if a[1] == j[1]:
                distance = j[0] - a[0]
                distance_list.append(distance)
    if len(distance_list) == 0:
        mean_distance = None
    else:
        mean_distance = statistics.mean(distance_list)
    fitness = mean_distance
    return(fitness)

def mutation_func(offspring, r_mut):
    if rand() < r_mut:
        index_1 = random.randint(0, offspring.shape[0] - 1)
        index_2 = random.randint(0, offspring.shape[0] - 1)
        codon_1 = offspring[index_1,1]
        codon_2 = offspring[index_2,1]
        offspring[index_1,1] = codon_2
        offspring[index_2,1] = codon_1

    return offspring

def create_initial_pop(array, n_mut, n_pop = 10):
    pop = []
    for j in range(n_pop):
        for i in range(n_mut):
            array = mutation_func(array, r_mut = 1)
        pop.append(array)
    return(pop)

def selection(pop, scores, k=3):
    # first random selection
    selection_ix = randint(len(pop))
    for ix in randint(0, len(pop), k-1):
        # check if better (e.g. perform a tournament)
        if scores[ix] or scores[selection_ix] is None:
            break
        if scores[ix] < scores[selection_ix]:
            selection_ix = ix
    return pop[selection_ix]

def genetic_algorithm(input_array, fitness_function, n_iter, n_pop, r_mut, max_s = 20):
    pop = create_initial_pop(array = input_array, n_mut = len(input_array), n_pop = n_pop)
    t_0 = time.perf_counter()
    best = pop[0]
    best_eval = fitness_function(pop[0])
    for gen in range(n_iter):
        # evaluate all candidates in the population
        scores = [fitness_function(c) for c in pop]
        # check for new best solution
        for i in range(n_pop):
            if best_eval is not None:
                if scores[i] > best_eval:
                    best, best_eval = pop[i], scores[i]
                    print(">New best, Generation: %s with mean: %.3f" % (gen, scores[i]))
        selected = [selection(pop, scores) for _ in range(n_pop)]
        # create the next generation
        children = list()
        for i in range(0, n_pop, 2):
            # get selected parents in pairs
            p1, p2 = selected[i], selected[i+1]
            # crossover and mutation
            for c in [p1, p2]:
                # mutation
                mutation_func(c, r_mut)
                # store for next generation
                children.append(c)
        # replace population
        pop = children
        t_n = time.perf_counter()
        #if (t_0 - t_n) > max_s:
            #break
    #print((t_n - t_0))
    return [best, best_eval]

def codon_opt(AA_seq, profile, n_pop = 100, n_iter = 100, fitness_function = mean_dist_fitness_func_arr):
    look_up_array = AA_look_up(AA_seq = AA_seq, profile = profile)
    sub_arrays = subsequence_array(AA_seq = AA_seq, look_up_array = look_up_array)
    complete_seq = []
    for sub_arr in sub_arrays:
        sub_seq = genetic_algorithm(input_array = sub_arr,
                                    fitness_function = fitness_function, 
                                    n_iter = n_iter, 
                                    n_pop = n_pop, 
                                    r_mut = (1.0 / float(len(sub_arr[:,0]))) * 0.5)
        complete_seq.append(sub_seq[0])
    final = np.concatenate(complete_seq, axis = 0)
    final = final[final[:, 0].argsort()]
    final = pd.DataFrame(final).rename(columns={0: "Index", 1: "Codon"})
    final = final.merge(profile.iloc[:,0:2], how = 'left', on = 'Codon')
    return(final)