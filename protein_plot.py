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

	def __init__(self,patient_id,last_vital_status,overall_survival,disease_free_survival,max_followup):
		self.patient_id = patient_id
		self.last_vital_status = last_vital_status
		self.overall_survival = overall_survival
		self.disease_free_survival = disease_free_survival
		self.protein_data = None
		self.max_followup = max_followup
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
		dict_patient[row['patient_id']] = Patient_data(row['patient_id'],row['last_vital_status'],row['overall_survival'],row['disease_free_survival'],row['max_followup'])

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

# Protein name

print('Available Proteins')
for index in range(len(protein_names)):
	print(index,protein_names[index])

p_index = int(input('Enter an index:'))
print('Protein Name:',protein_names[p_index])

# retrieve expressions by p_index

expressions = []
for key in dict_patient:
	a = (key,dict_patient[key].output_protein(p_index))
	if a[1] != None:
		expressions.append(a)

expressions.sort(key=itemgetter(1))


T_high = []
E_high = []
T_low = []
E_low = []



# partition

print('Available partition types:\n0\t1/2\n1\t1/4')
selection = int(input('Choose:'))
if selection == 0:
	range_low = range(int(len(expressions)/2))
	range_high = range(int(len(expressions)/2),len(expressions)) 
	partition_title = '1/2_Partition'
elif selection == 1:
	range_low = range(int(len(expressions)/4))
	range_high = range(3*int(len(expressions)/4),len(expressions))
	partition_title = '1/4_Partition'
else:
	raise NameError('Unexpected partition type')

pprint(range_low)
pprint(range_high)

# calculate p value


# plot

print('Available plot types:\n0\tOverall_Survival_Curve\n1\tDisease_Free_Survival_Curve')
selection = int(input('Choose:'))
if selection == 0:
	for index in range_low:
		check = False
		if dict_patient[expressions[index][0]].overall_survival != 'unknown':
			T_low.append(int(dict_patient[expressions[index][0]].overall_survival))
			check = True
		elif dict_patient[expressions[index][0]].max_followup != 'unknown':
			T_low.append(int(dict_patient[expressions[index][0]].max_followup))
			check = True
		if check == True:
			if dict_patient[expressions[index][0]].last_vital_status == 'Dead':
				E_low.append(1)
			else:
				E_low.append(0)
	for index in range_high:
		check = False
		if dict_patient[expressions[index][0]].overall_survival != 'unknown':
			T_high.append(int(dict_patient[expressions[index][0]].overall_survival))
			check = True
		elif dict_patient[expressions[index][0]].max_followup != 'unknown':
			T_high.append(int(dict_patient[expressions[index][0]].max_followup))
			check = True
		if check == True:
			if dict_patient[expressions[index][0]].last_vital_status == 'Dead':
				E_high.append(1)
			else:
				E_high.append(0)
	plot_title = 'Overall_Survival_Curve '
elif selection == 1:
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
	plot_title = 'Disease_Free_Survival_Curve '
else:
	raise NameError('Unexpected plot type')

print (T_high)
print (E_high)
print (T_low)
print (E_low)



results = logrank_test(T_high, T_low, E_high, E_low, alpha=0.95)
results.print_summary()

kmf = KaplanMeierFitter()

kmf.fit(T_high, event_observed=E_high, label='High')

kmf.survival_function_
kmf.median_
print(kmf)
ax = kmf.plot(show_censors=True)




kmf.fit(T_low, event_observed=E_low, label='Low')
kmf.plot(ax=ax,show_censors=True)

plt.title('Kaplan-Meier_'+plot_title+partition_title+'protein_names'+'Protein_Type:')
plt.show()

