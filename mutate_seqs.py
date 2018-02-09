import sys
import numpy
import shutil
import os
from datetime import datetime
import random


fas_file = "all"
stat_file = "all"
use_poisson = False

for i in range(1,len(sys.argv)):
	arg = sys.argv[i]
	if arg == "-f":
		fas_file = sys.argv[i+1]
	elif arg == "-s":
		stat_file = sys.argv[i+1]
	elif arg == "-p":
		use_poisson = True


print()
if use_poisson:
	print("USING POISSON DISTRIBUTION")
else:
	print("USING MARKOV CHAINS")
print()

substitution_matrix = {
	"A":{"G":0,"C":0,"T":0},
	"G":{"A":0,"C":0,"T":0},
	"C":{"A":0,"G":0,"T":0},
	"T":{"A":0,"G":0,"C":0}
	}
sub_probs = {"A":0.001,"G":0.001,"C":0.001,"T":0.001}

average_mutation_rate_per_year = 0

date_loc = -1
date_format = '%m/%d/%y'


def mutate_seq(original_seq,days):
	""" Takes a sequence and mutates it according to the amount of time that has passed and the probabilities in the substitution matrix """
	""" Months are measured in 30 day periods """
	new_seq = ""
	mutation_count = 0
	for letter in original_seq:
		prob = sub_probs[letter] * average_mutation_rate_per_year/365.0
		for day in range(days):
			r = random.random()
			if r < prob:
				sub_sum = 0
				r = random.random()
				for sub in substitution_matrix[letter]:
					row_prob = substitution_matrix[letter][sub]
					if r < row_prob+sub_sum:
						letter = sub
						mutation_count += 1
						break
					else:
						sub_sum += row_prob
					
		new_seq += letter
	print("Mutated sequence for " + str(days) + " days. Made " + str(mutation_count) + " mutations. Expected " + str((average_mutation_rate_per_year/365)*len(original_seq)*days) + " mutations.")
	return new_seq
			
		
def mutate_seq_poisson(original_seq,days):
	""" Takes a sequence and mutates it according to the amount of time that has passed and the probabilities in the substitution matrix using a poisson estimation"""
	""" Months are measured in 30 day periods """
	new_seq = original_seq
	mutation_count = numpy.random.poisson((average_mutation_rate_per_year/365)*len(original_seq)*days,1)[0]

	seq_list = list(range(len(original_seq)))
	seq_weights = []
	for letter in original_seq:
		seq_weights.append(sub_probs[letter])
	
	seq_weights_sum = sum(seq_weights)
	seq_weights = [x/seq_weights_sum for x in seq_weights]

	replace_indicies = list(numpy.random.choice(seq_list,mutation_count, p=seq_weights,replace=False))

	new_seq = original_seq
	for entry in replace_indicies:
		letter = new_seq[entry]
		letter_choices = [x for x in substitution_matrix[letter]]
		letter_weights = [substitution_matrix[letter][x] for x in letter_choices]
		new_letter = numpy.random.choice(letter_choices,p=letter_weights)
		new_seq = new_seq[:entry] + new_letter + new_seq[entry+1:]

	print("Mutated sequence for " + str(days) + " days. Made " + str(mutation_count) + " mutations. Expected " + str((average_mutation_rate_per_year/365)*len(original_seq)*days) + " mutations.")
	return new_seq

	



	

	for letter in original_seq:
		prob = sub_probs[letter]
		for day in range(days):
			r = random.random()
			if r < prob:
				sub_sum = 0
				r = random.random()
				for sub in substitution_matrix[letter]:
					row_prob = substitution_matrix[letter][sub]
					if r < row_prob+sub_sum:
						letter = sub
						mutation_count += 1
						break
					else:
						sub_sum += row_prob
					
		new_seq += letter
	print("Mutated sequence for " + str(days) + " days. Made " + str(mutation_count) + " mutations. Expected " + str((average_mutation_rate_per_year/365)*len(original_seq)*days) + " mutations.")
	return new_seq
	
def make_substitution_matrix(file_name):
	""" Reads the input substitution matrix and applies it """
	global average_mutation_rate_per_year
	global date_loc
	global date_format
	given_sub_matrix = {
		"A":{"G":1,"C":1,"T":1},
		"G":{"A":1,"C":1,"T":1},
		"C":{"A":1,"G":1,"T":1},
		"T":{"A":1,"G":1,"C":1}
		}
	given_sub_probs = {"A":1,"G":1,"C":1,"T":1}
	with open(file_name,"r") as f:
		lines = f.readlines()
		reversible = lines[0].strip() == "reversible"
		start = 0
		if reversible:
			start = 1
		for line in lines[start:]:
			if line.startswith("mutation rate"):
				average_mutation_rate_per_year = float(line.strip().split(":")[1])
				continue

			if line.startswith("date location"):
				date_loc = int(line.strip().split(":")[1])
				continue	

			if line.startswith("date format"):
				date_format = line.strip().split(":")[1]
				continue
				
			line_parts = line.strip().split(":")
			first = line_parts[0][0]
			second = line_parts[0][1]
			number = float(line_parts[1])
			given_sub_matrix[first][second] = number
			if reversible:
				given_sub_matrix[second][first] = number

		matrix_sum = sum([sum([given_sub_matrix[x][y] for y in given_sub_matrix[x]]) for x in given_sub_matrix])
		for x in given_sub_matrix:
			row_sum = sum([given_sub_matrix[x][y] for y in given_sub_matrix[x]])
			given_sub_probs[x] = (row_sum/matrix_sum)*4
			for y in given_sub_matrix[x]:
				given_sub_matrix[x][y] = given_sub_matrix[x][y]/row_sum



		return given_sub_probs,given_sub_matrix


fas_files = []
stat_files = []
if fas_file == "all":
	for f in os.listdir():
		if f.endswith(".fas") and not f.endswith("_mutated.fas"):
			fas_files.append(f)
	if stat_file == "all":
		for f in os.listdir():
			if f.endswith(".sta"):
				stat_files.append(f)
	else:
		stat_files = [stat_file]
else:
	fas_files = [fas_file]
		

if os.path.exists('mutated/'):
        shutil.rmtree('mutated/')
os.mkdir('mutated/')

for fas_f in fas_files:
	if stat_files == "all":
		cur_stat = [x for x in stat_files if x[:-4] == f[:-4]][0]
	else:
		cur_stat = stat_file

	print("Mutating " + fas_f + " under the stats found in " + stat_file + " ...")
	print("Average probality of mutation per base pair per year: " + str(average_mutation_rate_per_year))


	sub_probs, substitution_matrix = make_substitution_matrix(cur_stat)
			

	seq_names = []
	seq_letters = []
	with open(fas_f,"r") as f:
		for line in f:
			line = line.strip()
			if line.startswith(">"):
				seq_names.append(line[1:])
			else:
				seq_letters.append(line.replace("-",""))

	seq_dates = []
	for name in seq_names:
		raw_date = name.split("_")[date_loc]
		seq_dates.append(datetime.strptime(raw_date, date_format))


	index_list = [i[0] for i in sorted(enumerate(seq_dates), key=lambda x:x[1])]


	prev_date = seq_dates[index_list[0]]
	prev_letters = seq_letters[index_list[0]] 
	prev_name = seq_names[index_list[0]]

	new_file_name = "mutated/" + fas_f[:-4] + "_mutated.fas"
	with open(new_file_name,"w") as f:
		f.write(">" + prev_name + "\n")
		f.write(prev_letters + "\n")

	print("Sequence length: " + str(len(prev_letters)))

	for i in index_list[1:]:
		cur_date = seq_dates[i]
		cur_name = seq_names[i]
		
		day_number = (cur_date - prev_date).days

		if use_poisson:
			cur_letters = mutate_seq_poisson(prev_letters,day_number)
		else:
			cur_letters = mutate_seq(prev_letters,day_number)

		with open(new_file_name,"a") as f:
			f.write(">" + cur_name + "\n")
			f.write(cur_letters + "\n")
		
		prev_date = cur_date
		prev_letters = cur_letters
	print()



