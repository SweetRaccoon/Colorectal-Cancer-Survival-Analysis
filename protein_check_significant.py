import json
import os
import csv
from pprint import pprint
from xml.etree import ElementTree
from operator import itemgetter
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from matplotlib import pyplot as plt

# class definition 

class Patient_data(object):
	pat_count = 0;

	def __init__(self,patient_id,max_followup,min_new_tumor,max_death,last_vital_status):
		self.patient_id = patient_id
		self.max_followup = max_followup
		self.min_new_tumor = min_new_tumor
		self.max_death = max_death
		self.last_vital_status = last_vital_status
		self.protein_data = None
		Patient_data.pat_count += 1

	def display_count(self):
		print(Patient_data.pat_count)

	def display_data(self):
		pprint(self.__dict__)

	def output_protein(self,index):
		if self.protein_data != None:
			return float(self.protein_data[protein_names[index]])



# data definition

dict_patient = {}

# initialize

with open('data.csv', 'r') as f:
	f_csv = csv.DictReader(f)
	for row in f_csv:
		dict_patient[row['patient_id']] = Patient_data(row['patient_id'],row['max_followup'],row['min_new_tumor'],row['max_death'],row['last_vital_status'])

with open('L3.csv', 'r') as f:
	f_csv = csv.reader(f)
	protein_names = next(f_csv)
	for x in range(0,3):
		protein_names.pop(0)

with open('L3.csv', 'r') as f:
	f_csv = csv.DictReader(f)
	i = 0
	for row in f_csv:
		pat_id = row['Sample_ID'][8:12]
		if pat_id in dict_patient:
			i += 1
			dict_patient[pat_id].protein_data = row
	print('Total:',i,'found match')




	

def km_test(protein_index,mode):
	partition_type = int(mode%2)
	curve_type = int(mode/2)
	expressions = []
	for key in dict_patient:
		a = (key,dict_patient[key].output_protein(protein_index))
		if a[1] != None:
			expressions.append(a)

	expressions.sort(key=itemgetter(1))

	T_high = []
	E_high = []
	T_low = []
	E_low = []

	# partition


	if partition_type == 0:
		range_low = range(int(len(expressions)/2))
		range_high = range(int(len(expressions)/2),len(expressions)) 
	elif partition_type == 1:
		range_low = range(int(len(expressions)/4))
		range_high = range(3*int(len(expressions)/4),len(expressions))


	# T and E array

	if curve_type == 0:
		for index in range_low:
			if dict_patient[expressions[index][0]].max_death != 'unknown':
				T_low.append(int(dict_patient[expressions[index][0]].max_death))
				if dict_patient[expressions[index][0]].last_vital_status == 'Dead':
					E_low.append(1)
				else:
					E_low.append(0)
		for index in range_high:
			if dict_patient[expressions[index][0]].max_death != 'unknown':
				T_high.append(int(dict_patient[expressions[index][0]].max_death))
				if dict_patient[expressions[index][0]].last_vital_status == 'Dead':
					E_high.append(1)
				else:
					E_high.append(0)
	elif curve_type == 1:
		for index in range_low:
			if dict_patient[expressions[index][0]].max_followup != 'unknown':
				T_low.append(int(dict_patient[expressions[index][0]].max_followup))
				if dict_patient[expressions[index][0]].min_new_tumor == 'unknown':
					E_low.append(0)
				else:
					E_low.append(1)
		for index in range_high:
			if dict_patient[expressions[index][0]].max_followup != 'unknown':
				T_high.append(int(dict_patient[expressions[index][0]].max_followup))
				if dict_patient[expressions[index][0]].min_new_tumor == 'unknown':
					E_high.append(0)
				else:
					E_high.append(1)


	result = logrank_test(T_high, T_low, E_high, E_low, alpha=0.95)
	return result



rows = []
headings = ['Protein','Partition_type','Curve_type','p_value','test_statistic','test_result','is_significant']
curve_name = ['Overall','Disease_Free']
partition_name = ['half','quarter']

i = 0
for index in range(len(protein_names)):
	for mode in range(0,4):
		result = km_test(index,mode)
		row = []
		row.append(protein_names[index])
		row.append(partition_name[int(mode%2)])
		row.append(curve_name[int(mode/2)])
		row.append(result.p_value)
		row.append(result.test_statistic)
		row.append(result.test_result)
		row.append(result.is_significant)
		rows.append(row)

with open('protein.csv', 'w') as f:
	f_csv = csv.writer(f)
	f_csv.writerow(headings)
	f_csv.writerows(rows)
