import json
import os
import csv
from pprint import pprint
from xml.etree import ElementTree
from operator import itemgetter
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from matplotlib import pyplot as plt

class Patient_data(object):
	pat_count = 0;

	def __init__(self,patient_id,last_vital_status,overall_survival,disease_free_survival,min_new_tumor):
		self.patient_id = patient_id
		self.last_vital_status = last_vital_status
		self.overall_survival = overall_survival
		self.disease_free_survival = disease_free_survival
		self.min_new_tumor = min_new_tumor
		self.rna_expressions = {}
		self.expressions = {}
		self.filename = None
		Patient_data.pat_count += 1

	def display_count(self):
		print(Patient_data.pat_count)

	def display_data(self):
		pprint(self.__dict__)

	def set_filename(self,fn):
		self.filename = fn

	def retrieve_rna_expression(self):
		if self.filename != None:
			with open('RNA-Seq/Exp_Data/'+self.filename,'r') as f:
				expressions = csv.reader(f, delimiter='\t')
				for row in expressions:
					self.rna_expressions[rna_names[row[0][:15]]] = float(row[1])

				return True
		else:
			return False

	def retrive_expressions(self):
		if self.filename != None:
			with open('RNA-Seq/Exp_Data/'+self.filename,'r') as f:
				expressions = csv.reader(f, delimiter='\t')
				for row in expressions:
					self.rna_expressions[rna_names[row[18][:30]]] = float(row[1])

				return True
		else:
			return False



# data definition

dict_patient = {}

# initialize

with open('data.csv', 'r') as f:
	f_csv = csv.DictReader(f)
	for row in f_csv:
		dict_patient[row['patient_id']] = Patient_data(row['patient_id'],row['last_vital_status'],row['overall_survival'],row['disease_free_survival'],row['min_new_tumor'])

with open('RNA-Seq/metadata.json','r') as f_data:
	data = json.load(f_data)

	for index in range(len(data)):
		cdata = data[index]
		fn = cdata['file_name']
		if fn.endswith('FPKM-UQ.txt.gz'):
			fn = fn[0:48]
			pat_id = cdata['cases'][0]['submitter_id'][8:12]
			if pat_id in dict_patient:
				dict_patient[pat_id].set_filename(fn)

rna_names = {}
i = 0

with open('gene_table.csv','r') as f:
	rows = csv.reader(f)
	for row in rows:
		rna_names[row[0]] = row[1]

i = 0
for key in dict_patient:
	result = dict_patient[key].retrieve_rna_expression()
	if result == True:
		i += 1
		print('Patient ID:',key,'retrieve successful')
	else:
		print('Patient ID:',key,'No data')

print(i,'found match')


expressions = []
for key in dict_patient:
	a = (key,dict_patient[key].retrieve_expressions())
	if a[1] != None:
		expressions.append(a)

expressions.sort(key=itemgetter(1))



def km_test(rna_index,mode):
	partition_type = int(mode%2)
	curve_type = int(mode/2)
	expressions = []
	for key in dict_patient:
		a = (key,dict_patient[key].retrieve_expressions())
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
			if dict_patient[expressions[index][0]].overall_survival != 'unknown':
				T_low.append(int(dict_patient[expressions[index][0]].overall_survival))
				if dict_patient[expressions[index][0]].last_vital_status == 'Dead':
					E_low.append(1)
				else:
					E_low.append(0)
		for index in range_high:
			if dict_patient[expressions[index][0]].overall_survival != 'unknown':
				T_high.append(int(dict_patient[expressions[index][0]].overall_survival))
				if dict_patient[expressions[index][0]].last_vital_status == 'Dead':
					E_high.append(1)
				else:
					E_high.append(0)
	elif curve_type == 1:
		for index in range_low:
			if dict_patient[expressions[index][0]].disease_free_survival != 'unknown':
				T_low.append(int(dict_patient[expressions[index][0]].disease_free_survival))
				if dict_patient[expressions[index][0]].min_new_tumor == 'unknown':
					E_low.append(0)
				else:
					E_low.append(1)
		for index in range_high:
			if dict_patient[expressions[index][0]].disease_free_survival != 'unknown':
				T_high.append(int(dict_patient[expressions[index][0]].disease_free_survival))
				if dict_patient[expressions[index][0]].min_new_tumor == 'unknown':
					E_high.append(0)
				else:
					E_high.append(1)


	result = logrank_test(T_high, T_low, E_high, E_low, alpha=0.95)
	return result



rows = []
headings = ['Rna_names''Partition_type','Curve_type','p_value','test_statistic','test_result','is_significant']
curve_name = ['Overall','Disease_Free']
partition_name = ['half','quarter']

i = 0
for index in range(len(rna_names)):
	for mode in range(0,4):
		result = km_test(index,mode)
		row = []
		row.append(rna_names)
		row.append(partition_name[int(mode%2)])
		row.append(curve_name[int(mode/2)])
		row.append(result.p_value)
		row.append(result.test_statistic)
		row.append(result.test_result)
		row.append(result.is_significant)
		rows.append(row)
		pprint(row)

with open('rna.csv', 'w') as f:
	f_csv = csv.writer(f)
	f_csv.writerow(headings)
	f_csv.writerows(rows)




















